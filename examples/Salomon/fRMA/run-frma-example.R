## Load packages
library(openxlsx)
library(affy)
library(frma)

## Load the fRMA data package and then load the data
library(DSalomon.PAX.hthgu133pluspmfrmavecs)
data(DSalomon.PAX.hthgu133pluspmfrmavecs)

## Get the sample table and list of CEL files
cel.dir <- "Training Data/01 - CTOT08 ABECASSIS PAX Samples"
## Read the first xlsx file in the directory, which is a spreadsheet
## containing all of the CEL file names
sample.table <- read.xlsx(list.files(cel.dir, pattern=glob2rx("*.xlsx"), full.names=TRUE)[1])
cel.files <- file.path(cel.dir, sample.table$Filename)

## Read the data from the CEL files
affy <- ReadAffy(filenames=cel.files, phenoData=sample.table)

## Apply fRMA
eset <- frma(affy, input.vecs=DSalomon.PAX.hthgu133pluspmfrmavecs)

## Extract expression matrix
expr <- exprs(eset)

## Extract sample table
sample.table <- pData(eset)
