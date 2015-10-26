#!/usr/bin/env Rscript

# Script to train multiple fRMA vectors in preparation for consistency
# evaluation

library(xlsx)
library(frma)
library(frmaTools)
library(stringr)
library(magrittr)
library(plyr)
library(affy)
library(preprocessCore)
library(ggplot2)
library(proto)
library(dplyr)

training.data.dir <- "Training Data"
datasets <- data.frame(Dataset=list.files(training.data.dir))
rownames(datasets) <- datasets$Dataset
datasets$Tissue <- factor(str_extract(datasets$Dataset, "\\b(PAX|BX)\\b"))

tsmsg <- function(...) {
  message(date(), ": ", ...)
}

## Some Scan Dates are marked as identical for multiple batches, which
## is bad. But the dates embedded in the file names for these batches
## are different, so we use those dates instead.
parse.date.from.filename <- function(fname) {
    res1 <- str_match(fname, "^(\\d\\d)(\\d\\d)(\\d\\d)")[,c(4,2,3)]
    res2 <- str_match(fname, "^20(\\d\\d)_(\\d\\d)_(\\d\\d)")[,-1]
    res1[is.na(res1)] <- res2[is.na(res1)]
    colnames(res1) <- c("year", "month", "day")
    res1[,"year"] %<>% str_c("20", .)
    as.Date(do.call(ISOdate, data.frame(res1)))
}

## Error: the following are not valid files:
##     Training Data/03 - TGCG ARADNRCANTX PAX Samples/10733.CEL
blacklist <- "Training Data/03 - TGCG ARADNRCANTX PAX Samples/10733.CEL"

## This reads in the xlsx file for each of the 7 datasets and combines
## them into one big table of all samples. The Batch column contains
## the partitioning of samples into unique combinations of Dataset,
## Scan Date, and Phenotype. Finally, we split based on Tissue type to
## get one table for biopsies (BX), and one for blood (PAX).
sample.tables <- ddply(datasets, .(Dataset), function(df) {
    df <- df[1,]
    rownames(df) <- NULL
    dset.dir <- file.path(training.data.dir, df$Dataset)
    x <- read.xlsx(list.files(dset.dir, pattern=glob2rx("*.xlsx"), full.names=TRUE)[1], 1) %>%
        setNames(c("Filename", "Phenotype", "ScanDate"))
    x$Filename <- as.character(x$Filename)
    missing.CEL <- !str_detect(x$Filename, "\\.CEL$")
    x$Filename[missing.CEL] <- str_c(x$Filename[missing.CEL], ".CEL")
    stopifnot(all(str_detect(x$Filename, "\\.CEL$")))
    parsed.date <- parse.date.from.filename(x$Filename)
    x$ScanDate[!is.na(parsed.date)] <- parsed.date[!is.na(parsed.date)]
    x %>% cbind(df) %>%
        transform(Filename=file.path(dset.dir, Filename),
                  Batch=droplevels(Tissue:Dataset:factor(ScanDate):Phenotype)) %>%
                      subset(! Filename %in% blacklist) %>%
                          subset(!duplicated(Filename))
}) %>%
    split(.$Tissue) %>%
        lapply(droplevels)

## fRMA requires equal-sized batches, so for each batch size from 3 to
## 15, compute how many batches have at least that many samples.
x <- sapply(3:15, function(i) sapply(sample.tables, . %$% Batch %>% table %>% as.vector %>% {sum(. >= i)}))
colnames(x) <- 3:15

## Based on the above and the recommendations in the frmaTools paper,
## I chose 5 as the optimal batch size. This could be optimized
## empirically, though.
arrays.per.batch <- 5

vectors <- lapply(names(sample.tables), function(ttype) {
    stab <- sample.tables[[ttype]]
    tsmsg("Reading full dataset for ", ttype)
    affy <- ReadAffy(filenames=stab$Filename, sampleNames=rownames(stab))
    tsmsg("Getting reference normalziation distribution from full dataset for ", ttype)
    normVec <- normalize.quantiles.determine.target(pm(bg.correct.rma(affy)))
    rm(affy); gc()
    ## Set the random seed for reproducibility.
    set.seed(1986)

    lapply(1:5, function(i) {
        on.exit(gc())
        tsmsg("Starting training run number ", i, " for ", ttype)
        tsmsg("Selecting batches for ", ttype)
        ## Keep only batches with enough samples
        big.enough <- stab$Batch %>% table %>% .[.>= arrays.per.batch] %>% names
        stab <- stab[stab$Batch %in% big.enough,] %>% droplevels

        ## Sample an equal number of arrays from each batch
        subtab <- ddply(stab, .(Batch), function(df) {
            df[sample(seq(nrow(df)), size=arrays.per.batch),]
        })

        tsmsg("Making fRMA vectors")
        ## Make fRMA vectors, using normVec from full dataset
        res <- makeVectorsAffyBatch(subtab$Filename, subtab$Batch, normVec=normVec)

        tsmsg("Finished training run number ", i, " for ", ttype)
        res
    }) %>% setNames(., str_c("V", seq_along(.)))
}) %>% setNames(names(sample.tables))

saveRDS(vectors, "consistency-vectors.RDS")
save.image("consistency.rda")
## Continues in consistency-evaluate.R
