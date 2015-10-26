#!/usr/bin/env Rscript

library(xlsx)
library(frma)
library(frmaTools)
library(stringr)
library(magrittr)
library(plyr)
library(affy)
library(preprocessCore)

training.data.dir <- "Training Data"
datasets <- data.frame(Dataset=list.files(training.data.dir))
rownames(datasets) <- datasets$Dataset
datasets$Tissue <- factor(str_extract(datasets$Dataset, "\\b(PAX|BX)\\b"))

tsmsg <- function(...) {
  message(date(), ": ", ...)
}

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
                      subset(! Filename %in% blacklist)
}) %>% split(.$Tissue) %>% lapply(droplevels)

annotation <- cleancdfname(affyio:::read.celfile.header(sample.tables[[1]]$Filename[1])$cdfName, FALSE)

esets <- list()

for (i in names(sample.tables)) {
    pkgname <- sprintf("DSalomon.%s.%sfrmavecs", i, annotation)
    message("Loading ", pkgname)
    require(pkgname, character.only=TRUE, quietly=TRUE)
    data(list = pkgname)

    message("Loading raw data for ", i)
    stab <- sample.tables[[i]]
    affy <- ReadAffy(filenames=stab$Filename, phenoData=stab)

    message("Running fRMA for ", i)
    esets[[i]] <- frma(affy, input.vecs=get(pkgname), verbose=TRUE)
    rm(affy)
    gc()
}
