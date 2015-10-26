#!/usr/bin/Rscript

source("common.R")

library(rtracklayer)
library(DiffBind)
library(plyr)
library(doMC)
registerDoMC()
options(cores=multicore:::detectCores())
library(snow)
library(stringr)
library(RColorBrewer)
library(xlsx)
library(edgeR)

tsmsg <- function(...) {
  message(date(), ": ", ...)
}

tsmsgf <- function(...) {
  tsmsg(sprintf(...))
}

## Additional args are passed to every call of addDataFrame
write.xlsx.multisheet <- function(data.frames, file, sheetNames=names(data.frames), ...) {
  if (is.null(sheetNames)) {
    sheetNames <- str_c("Sheet", seq_along(data.frames))
  }
  ## Ensure correct number of sheetNames
  stopifnot(length(sheetNames) == length(data.frames))
  ## Fill in missing names if needed
  sheetNames[is.na(sheetNames)] <- str_c("Sheet", seq_along(data.frames))[is.na(sheetNames)]
  wb <- createWorkbook()
  sheets <- llply(sheetNames, function(x) createSheet(wb, sheetName=x))
  mlply(cbind(sheet=sheets, x=data.frames), .fun=addDataFrame, .parallel=FALSE, ...)
  saveWorkbook(wb, file)
}

renice.self <- function(niceness=19)  {
  system2("renice", args=c("-n", as.character(niceness), Sys.getpid()))
}

makeNiceCluster <- function(...) {
  cl <- makeCluster(...)
  clusterCall(cl, renice.self)
  cl
}

select.nearest <- function(x, y) {
  y[nearest(x,y)]
}

kgxref <- get.ucsc.table("knownGene","kgXref")
get.kgxref <- function(kgID) {
  x <- merge(x=DataFrame(kgID=kgID), y=kgxref,
             all.x=TRUE, all.y=FALSE)
  x[names(x) != "kgID"]
}

read.htseq.counts <- function(f) {
  x <- read.table(f, header=FALSE, sep="\t")
  names(x) <- c("ID", "count")
  ## Discard the last 5 lines
  exception.rows <- (nrow(x)-4):nrow(x)
  attr(x, "exception_counts") <- x[exception.rows,]
  x <- x[-exception.rows,]
  x
}

get.transcript.lengths <- function(transcript.ids, exon.lengths) {
  exons.by.transcript <- split(exon.lengths, transcript.ids)
  sapply(exons.by.transcript, sum)
}

get.gene.lengths <- function(gene.ids, transcript.lengths, method="max") {
  if (is.character(method)) {
    method <- get(method)
  }
  transcripts.by.genes <- split(transcript.lengths, gene.ids)
  sapply(transcripts.by.genes, method)
}

str_matches <- function(string, pattern) {
  !is.na(str_match(string, pattern)[,1])
}

annotate.ensp.in.hsap.or.mmul <- function(ensp.ids) {
  ensembl=useMart("ensembl")
  datasets <- c("mmulatta_gene_ensembl", "hsapiens_gene_ensembl")
  x <- Reduce(rbind, llply(datasets, function(x) {
    ensembl <- useDataset(x, mart=ensembl)
    getBM(filters="ensembl_peptide_id",
          values=unique(ensp.ids),
          attributes=c("hgnc_symbol", "wikigene_name", "ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id", "description", "wikigene_description"),
          mart=ensembl,
          uniqueRows=TRUE)
  }))
  ## Convert empty string to NA in all columns
  x <- data.frame(llply(x, function(column) ifelse(column == "", NA, column)))

  ## Unify symbols & descriptions
  unified.symbol <- ifelse(is.na(x$hgnc_symbol), as.character(x$wikigene_name), as.character(x$hgnc_symbol))
  unified.desc <- ifelse(is.na(x$description), as.character(x$wikigene_description), as.character(x$description))
  x <- data.frame(symbol=unified.symbol, x[! names(x) %in% c("hgnc_symbol", "wikigene_name", "description", "wikigene_description")], description=unified.desc)

  ## Reorder rows so that more-preferred rows for the same input are on top
  x <- x[order(is.na(x$symbol), str_matches(x$symbol, "^LOC\\d+$"), is.na(x$description)),]

  ## Keep only the first row for each input
  x <- x[!duplicated(x$ensembl_peptide_id),]

  data.frame(x, row.names=x$ensembl_peptide_id)
}

## TODO: Replace by biomart

CE.name.annot <- {
  x <- setNames(nm=c("CE_ensembl_peptide_id", "symbol", "description"),
                read.table("annotation/CE.name", sep="\t", quote=""))
  x$ensembl_peptide_id <- str_replace(x$CE_ensembl_peptide_id, "^CE_", "")
  row.names(x) <- x$ensembl_peptide_id
  x
}

annotation.gff <- import("cuffmerge_results/merged.gff")

transcripts <- annotation.gff[annotation.gff$type == "transcript", ]
exons <- annotation.gff[annotation.gff$type == "exon", ]
transcript.lengths <- get.transcript.lengths(unlist(exons$Parent), width(exons))
gene.lengths <- get.gene.lengths(transcripts$geneID, transcript.lengths[transcripts$ID])

gene.ref.lists <- sapply(split(transcripts$nearest_ref, transcripts$geneID),
                         function(x) unique(x[!is.na(x)]))

## For genes with multiple ref IDs, join them with commas
gene.refs <- sapply(gene.ref.lists, function (x) {
  if (length(x) == 0) {
    NA
  } else {
    str_c(x, collapse=",")
  }
})

## For genes with multiple ref IDs, just pick the first one
gene.first.refs <- sapply(gene.ref.lists, function(x) {
  if (length(x) == 0) {
    NA
  } else {
    x[[1]]
  }
})

## Use delayed assignment so that the cluster won't be created until it is actually needed
delayedAssign("cl", makeNiceCluster(rep(c("salomon14", "salomon18", "salomon19"), 8)))

mart.raw.annot <- annotate.ensp.in.hsap.or.mmul(na.omit(gene.first.refs))
gene.annot <- data.frame(mart.raw.annot[gene.first.refs,], row.names=names(gene.first.refs))

## Fill in missing values that are available from the CE.name file
## provided by the cyno genome paper authors
suppl.gene.annot <- data.frame(CE.name.annot[gene.first.refs,], row.names=names(gene.first.refs))

suppl.ensp <- is.na(gene.annot$ensembl_peptide_id) & !is.na(suppl.gene.annot$ensembl_peptide_id)
suppl.symbols <- is.na(gene.annot$symbol) & !is.na(suppl.gene.annot$symbol)
gene.annot$ensembl_peptide_id <- ifelse(suppl.ensp,
                                        as.character(suppl.gene.annot$ensembl_peptide_id),
                                        as.character(gene.annot$ensembl_peptide_id))
gene.annot$symbol <- ifelse(suppl.symbols,
                            as.character(suppl.gene.annot$symbol),
                            as.character(gene.annot$symbol))
## Replace description whenever we replace the symbol to keep things consistent
gene.annot$description <- ifelse(suppl.symbols,
                            as.character(suppl.gene.annot$description),
                            as.character(gene.annot$description))
gene.annot$symbol[is.na(gene.annot$ensembl_peptide_id)] <- "[No annotation]"
gene.annot$symbol[is.na(gene.annot$symbol)] <- "[No symbol for ENSP]"

expdata <- {
  x <- read.xlsx("kenyon-cyno-expdata.xls", 1)
  x$Animal.ID <- str_trim(x$Animal.ID)
  x$Condition <- str_replace_all(x$Condition, "γ", "g")
  x$Sample.ID <-
    sprintf("%s_%s_%s",
            str_replace_all(x$Animal.ID, "[ -]", "_"),
            ifelse(x$Condition == "Control", "CTRL", "IFNg"),
            x$Passage)
  x$Sample.Name <-
    sprintf("%s_%s_%s_L%03i",
            str_replace_all(x$Animal.ID, "[ -]", "_"),
            ifelse(x$Condition == "Control", 1, 2),
            x$Index.Seq,
            x$Lane)
  x$FASTQ.Read1.File <- sprintf("seqprep_results/%s/%s_R1_001.fastq", x$Sample.Name, x$Sample.Name)
  x$FASTQ.Read2.File <- sprintf("seqprep_results/%s/%s_R2_001.fastq", x$Sample.Name, x$Sample.Name)
  x$BAM.File <- sprintf("tophat_results/%s/accepted_hits.bam", x$Sample.Name)
  x$GFF.File <- sprintf("cufflinks_quantification/%s/transcripts.gff", x$Sample.Name)
  x$Counts.File <- sprintf("htseq_counts/%s/counts.txt", x$Sample.Name)
  mcparallel(write.xlsx(x, "expdata.xlsx"))
  rownames(x) <- x$Sample.ID
  x
}

counts.vectors <- setNames(llply(expdata$Counts.File, .parallel=TRUE, function(f) {
  x <- read.htseq.counts(f)
  setNames(x$count, x$ID)
}), expdata$Sample.ID)

## Put the data into a matrix, making sure we account for the
## possibility that not all genes are listed in all samples.
all.geneIDs <- sort(unique(unlist(llply(counts.vectors, names))))

counts <- {
  x <- matrix(data=0, ncol=length(counts.vectors), nrow=length(all.geneIDs),
              dimnames=list(`geneID`=all.geneIDs, `sample`=names(counts.vectors)))
  for (i in expdata$Sample.ID) {
    cvec <- counts.vectors[[i]]
    x[names(cvec),i] <- cbind(cvec)
  }
  x
}

blocked.design <- model.matrix(~Animal.ID+Passage+Condition, data=expdata)
unblocked.design <- model.matrix(~Condition, data=expdata)
dge <- DGEList(counts=counts,
               group=expdata$Condition,
               genes=gene.annot[rownames(counts),])
dge <- calcNormFactors(dge)

blocked.dge <- estimateGLMCommonDisp(dge, blocked.design, verbose=TRUE)
blocked.dge <- estimateGLMTrendedDisp(blocked.dge, blocked.design)
blocked.dge <- estimateGLMTagwiseDisp(blocked.dge, blocked.design)
blocked.fit <- glmFit(blocked.dge, blocked.design)
blocked.lrt <- glmLRT(blocked.dge, blocked.fit, coef="ConditionIFNg activated")
blocked.tt <- topTags(blocked.lrt, n=nrow(counts))

unblocked.dge <- estimateGLMCommonDisp(dge, unblocked.design, verbose=TRUE)
unblocked.dge <- estimateGLMTrendedDisp(unblocked.dge, unblocked.design)
unblocked.dge <- estimateGLMTagwiseDisp(unblocked.dge, unblocked.design)
unblocked.fit <- glmFit(unblocked.dge, unblocked.design)
unblocked.lrt <- glmLRT(unblocked.dge, unblocked.fit, coef="ConditionIFNg activated")
unblocked.tt <- topTags(unblocked.lrt, n=nrow(counts))

write.csv(blocked.tt$table, "edgeR-genes-prelim")
write.xlsx.multisheet(list(blocked=blocked.tt$table,
                           unblocked=unblocked.tt$table),
                      "edgeR-genes.xlsx")
