#!/usr/bin/env Rscript

library(xlsx)
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

makeVectorsAffyBatch <- function (files, batch.id, background = "rma", normalize = "quantile",
    normVec = NULL, cdfname = NULL, file.dir = ".", verbose = TRUE)
{
    wd <- getwd()
    setwd(file.dir)
    object <- ReadAffy(filenames = files, cdfname = cdfname,
        verbose = verbose)
    setwd(wd)
    if (verbose)
        message("Data loaded \n")
    batch.size <- table(batch.id)[1]
    if (!all(table(batch.id) == batch.size))
        stop("Batches must be of the same size.")
    if (background == "rma") {
        object <- bg.correct.rma(object)
        if (verbose)
            message("Background Corrected \n")
        gc()
    }
    pms <- pm(object)
    pns <- probeNames(object)
    pmi <- unlist(pmindex(object))
    if (!all(sprintf("%i", pmi) == rownames(pms)))
        stop("Mismatch between pmindex and rownames of pms")
    rm(object)
    gc()
    if (normalize == "quantile") {
        if (is.null(normVec))
            normVec <- normalize.quantiles.determine.target(pms)
        pms <- normalize.quantiles.use.target(pms, normVec)
        names(normVec) <- as.character(pmi)
        if (verbose)
            message("Normalized \n")
    }
    pms <- log2(pms)
    gc()
    N <- 1:dim(pms)[1]
    S <- split(N, pns)
    nc <- ncol(pms)
    nr <- nrow(pms)
    resids <- matrix(ncol = nc, nrow = nr)
    probeVec <- vector(length = nr)
    if (verbose)
        message("Beginning Probe Effect Calculation ... \n")
    for (k in 1:length(S)) {
        fit <- rcModelPLM(pms[S[[k]], , drop = FALSE])
        resids[S[[k]], ] <- fit$Residuals
        probeVec[S[[k]]] <- fit$Estimates[(nc + 1):length(fit$Estimates)]
        if ((k%%1000) == 0) {
            message(paste("Finished probeset:", k, "\n"))
            gc()
        }
    }
    names(probeVec) <- as.character(pmi)
    if (verbose)
        message("Probe Effects Calculated \n")
    gc()
    tmp <- split(t(resids), batch.id)
    withinMean <- lapply(tmp, frmaTools:::getProbeMean, batch.size)
    withinVar <- lapply(tmp, frmaTools:::getProbeVar, batch.size)
    withinAvgVar <- rowMeans(matrix(unlist(withinVar), ncol = length(withinVar)))
    btwVar <- apply(matrix(unlist(withinMean), ncol = length(withinMean)),
        1, var)
    rm(tmp)
    rm(withinMean)
    rm(withinVar)
    names(withinAvgVar) <- names(btwVar) <- as.character(pmi)
    if (verbose)
        message("Probe Variances Calculated \n")
    gc()
    tmp <- split(resids, pns)
    psetMAD <- unlist(lapply(tmp, frmaTools:::getPsetMAD, nc, batch.id))
    names(psetMAD) <- names(tmp)
    rm(tmp)
    rm(resids)
    if (verbose)
        message("Probe Set SDs Calculated \n")
    gc()
    w <- 1/(withinAvgVar + btwVar)
    w[w == Inf] <- 1
    medianSE <- vector(length = length(psetMAD))
    if (verbose)
        message("Beginning Median SE Calculation ... \n")
    for (k in 1:length(S)) {
        fit <- frmaTools:::rwaFit2(pms[S[[k]], , drop = FALSE], w[S[[k]]],
            probeVec[S[[k]]], psetMAD[k])
        medianSE[k] <- median(fit$StdErrors)
        if ((k%%1000) == 0) {
            message(paste("Finished probeset:", k, "\n"))
            gc()
        }
    }
    names(medianSE) <- names(psetMAD)
    if (verbose)
        message("Median SEs Calculated \n")
    gc()
    rm(w)
    rm(pms)
    rm(pns)
    gc()
    if (is.null(cdfname)) {
        vers <- ""
    } else {
        vers <- as.character(packageVersion(cdfname))
    }
    ## vers <- ifelse(!is.null(cdfname), as.character(packageVersion(cdfname)),
    ##     "")
    return(list(normVec = normVec, probeVec = probeVec, probeVarWithin = withinAvgVar,
        probeVarBetween = btwVar, probesetSD = psetMAD, medianSE = medianSE,
        version = vers))
}

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

## For each tissue type, compute fRMA vectors.
vectors <- lapply(sample.tables, function(stab) {
    set.seed(1986)

    tsmsg("Reading full dataset")
    affy <- ReadAffy(filenames=stab$Filename, sampleNames=rownames(stab))
    tsmsg("Getting reference normalization distribution from full dataset")
    normVec <- normalize.quantiles.determine.target(pm(bg.correct.rma(affy)))
    rm(affy); gc()

    tsmsg("Selecting batches")
    ## Keep only arrays with enough samples
    big.enough <- stab$Batch %>% table %>% .[.>= arrays.per.batch] %>% names
    stab <- stab[stab$Batch %in% big.enough,] %>% droplevels

    ## Sample an equal number of arrays from each batch
    subtab <- ddply(stab, .(Batch), function(df) {
        df[sample(seq(nrow(df)), size=arrays.per.batch),]
    })

    tsmsg("Making fRMA vectors")
    ## Make fRMA vectors, using normVec from full dataset
    res <- makeVectorsAffyBatch(subtab$Filename, subtab$Batch, normVec=normVec)

    tsmsg("Finished.")
    res
})

## The code below here just takes the trained vectors and packages
## them up into installable R packages.
makePackageFromVectors <-
    function (vecs, version, maintainer, species, annotation,
              packageName, file.dir = ".",
              output.dir = ".", unlink = TRUE)
{
    platform <- gsub("cdf$", "", annotation)
    ## type <- match.arg(type, c("AffyBatch", "FeatureSet"))
    ## if (type == "AffyBatch")
    ##     platform <- gsub("cdf", "", annotation)
    ## if (type == "FeatureSet") {
    ##     platform <- annotation
    ##     require(oligo)
    ## }
    thispkg <- "frmaTools"
    desc <- packageDescription(thispkg)
    thispkgVers <- desc$Version
    symbolValues <- list(ARRAYTYPE = platform, VERSION = version,
        CREATOR = paste("package", thispkg, "version", thispkgVers),
        FRMATOOLSVERSION = thispkgVers, MAINTAINER = maintainer,
        SPECIES = species)
    createdPkg <- createPackage(packageName, destinationDir = output.dir,
        originDir = system.file("VectorPkg-template", package = "frmaTools"),
        symbolValues = symbolValues, unlink = unlink)
    ## if (type == "AffyBatch")
    ##     vecs <- makeVectorsAffyBatch(files, batch.id, background,
    ##         normalize, normVec, annotation, file.dir, verbose)
    ## if (type == "FeatureSet")
    ##     vecs <- makeVectorsFeatureSet(files, batch.id, annotation,
    ##         background, normalize, normVec, file.dir, verbose)
    assign(packageName, vecs)
    save(list = eval(packageName), file = file.path(createdPkg$pkgdir,
        "data", paste(packageName, ".rda", sep = "")), compress = TRUE)
}

annotation <- cleancdfname(affyio:::read.celfile.header(sample.tables[[1]]$Filename[1])$cdfName, FALSE)

dir.create("pkgs", FALSE, TRUE, mode="755")

for (i in names(vectors)) {
    vecs <- vectors[[i]]
    pkgname <- sprintf("DSalomon.%s.%sfrmavecs", i, annotation)
    message("Making ", pkgname)
    makePackageFromVectors(
        vecs,
        version="0.1",
        maintainer="Ryan C. Thompson <rcthomps@scripps.edu>",
        species="Homo_sapiens",
        annotation=annotation,
        packageName=pkgname,
        output.dir = "pkgs")
}

save.image("train-data.rda")
