#!/usr/bin/env Rscript

source("common.R")
library(BSgenome.Hsapiens.UCSC.hg19)

textConnectionFromLines <- function(lines, linesep="\n") {
  textConnection(str_c(sprintf("%s%s", lines, linesep), collapse=""))
}

## Takes a GRangesList
calculate.block.entropy <- function(grl, expr.column="tagExpression") {
  groupfac <- factor(rep(names(grl), elementLengths(grl)))
  exprs <- as.vector(elementMetadata(unlist(grl))[[expr.column]])
  total.exprs <- aggregate(exprs, by=list(groupfac), FUN=sum)$x
  qi <- exprs / rep(total.exprs, elementLengths(grl))
  qi.times.log <- qi * log2(qi)
  results <- -aggregate(qi.times.log, by=list(groupfac), FUN=sum)$x
  names(results) <- names(grl)
  results
}

## http://hoffmann.bioinf.uni-leipzig.de/LIFE/blockbuster.html
read.blockbuster.tag.output <- function(bbout) {
  x <- readLines(bbout)
  cluster.line.numbers <- which(str_sub(x, 1,1) == ">")
  cluster.lines <- str_sub(x[cluster.line.numbers], 2)
  cluster.table <- read.table(textConnectionFromLines(cluster.lines),
                              col.names=c("clusterID", "chrom", "clusterStart", "clusterEnd", "strand", "ClusterExpression", "tagCount", "blockCount"),
                              colClasses=c("character", "factor", "integer", "integer", "factor", "numeric", "integer", "integer"))
  tag.lines <- x[-cluster.line.numbers]
  cluster.tag.line.counts <- cluster.line.numbers[-1] - cluster.line.numbers[-length(cluster.line.numbers)] - 1
  cluster.tag.line.counts <- c(cluster.tag.line.counts, length(tag.lines) - sum(cluster.tag.line.counts))
  tag.table <- read.table(textConnectionFromLines(tag.lines),
                          col.names=c("tagChrom", "tagStart", "tagEnd", "tagID", "tagExpression", "tagStrand", "blockNb"),
                          colClasses=c("factor", "integer", "integer", "character", "numeric", "factor", "integer"))
  tag.table$clusterID <- rep(cluster.table$clusterID, cluster.tag.line.counts)
  tag.table$blockID <- sprintf("%s/B%s", tag.table$clusterID, tag.table$blockNb)
  tags <- table.to.granges(tag.table, seqnames.column="tagChrom", start.column="tagStart", end.column="tagEnd", strand.column="tagStrand", seqlengths="hg19")
  tags <- split(tags, as.vector(elementMetadata(tags)$blockID))
  return(tags)
}

read.blockbuster.output <- function(bbout, bbout.tags=NULL) {
  x <- readLines(bbout)
  cluster.line.numbers <- which(str_sub(x, 1,1) == ">")
  cluster.lines <- str_sub(x[cluster.line.numbers], 2)
  block.lines <- x[-cluster.line.numbers]
  cluster.table <- read.table(textConnectionFromLines(cluster.lines),
                              col.names=c("clusterID", "chrom", "clusterStart", "clusterEnd", "strand", "ClusterExpression", "tagCount", "blockCount"),
                              colClasses=c("character", "factor", "integer", "integer", "factor", "numeric", "integer", "integer"))
  clusters <- table.to.granges(cluster.table, seqnames.column="chrom", start.column="clusterStart", end.column="clusterEnd", seqlengths="hg19")
  block.table <- read.table(textConnectionFromLines(block.lines),
                            col.names=c("blockNb", "blockChrom", "blockStart", "blockEnd", "blockStrand", "blockExpression", "readCount"),
                            colClasses=c("integer", "factor", "integer", "integer", "factor", "numeric", "integer"))
  block.table$clusterID <- rep(cluster.table$clusterID, cluster.table$blockCount)
  block.table$blockID <- sprintf("%s/B%s", block.table$clusterID, block.table$blockNb)
  blocks. <- table.to.granges(block.table, seqnames.column="blockChrom", start.column="blockStart", end.column="blockEnd", strand.column="blockStrand", seqlengths="hg19")
  ## blocks. <- split(blocks., as.vector(elementMetadata(blocks.)$clusterID))
  retval <- list(clusters=clusters, blocks=blocks.)
  if (!is.null(bbout.tags)) {
    retval$tags <- read.blockbuster.tag.output(bbout.tags)
  }
  return(retval)
}


blockbuster.path <- "/home/ryan/bin/blockbuster"

run.blockbuster <- function(gr, ...) {
  temp.bed.file <- tempfile(fileext=".bed")
  temp.bbout.file <- tempfile(fileext=".bbout")
  temp.bbout.tag.file <- tempfile(fileext=".tag.bbout")
  tryCatch({
    export(sort(gr), temp.bed.file, format="BED")
    extra.args <- list(...)
    bbargs <- c(rbind(sprintf("-%s", names(extra.args)), extra.args), "-print", "1", temp.bed.file)
    bbargs.tags <- c(rbind(sprintf("-%s", names(extra.args)), extra.args), "-print", "2", temp.bed.file)
    system2(blockbuster.path, args=bbargs,
            stdout=temp.bbout.file)
    system2(blockbuster.path, args=bbargs.tags,
            stdout=temp.bbout.tag.file)
    x <- read.blockbuster.output(temp.bbout.file, temp.bbout.tag.file)
    x
  }, finally={
    unlink(c(temp.bed.file, temp.bbout.file, temp.bbout.tag.file))
  })
}

## Load the reads from the bam file
infile <- "./all-results.bam"
read.ranges <- {
  x <- readAlignedRanges(infile, include.reads=TRUE)
  ## Throw away unneeded columns
  elementMetadata(x) <- subset(elementMetadata(x), select=-c(id,qual,flag))
  read.multi.map.counts <- table(as.vector(elementMetadata(x)$seq))
  elementMetadata(x)$multimap <- Rle(as.vector(read.multi.map.counts[as.vector(elementMetadata(x)$seq)]))
  x
}

## Get needed annotations

## rRNA & tRNA (from the repeats table)
repeat.table <- get.ucsc.table("rmsk", "rmsk", genome="hg19")
repeat.ranges <- table.to.granges(repeat.table, seqnames.column="genoName", start.column="genoStart", end.column="genoEnd", seqlengths="hg19")
trna.ranges <- repeat.ranges[elementMetadata(repeat.ranges)$repClass == "tRNA"]
rrna.ranges <- repeat.ranges[elementMetadata(repeat.ranges)$repClass == "rRNA"]

## miRNA & snoRNA
small.rna.table <- subset(get.ucsc.table("wgRna", "wgRna", genome="hg19"), select=-c(thickStart,thickEnd,bin))
small.rna.ranges <- table.to.granges(small.rna.table, seqnames.column="chrom", start.column="chromStart", end.column="chromEnd", seqlengths="hg19")
miRNA.types <- "miRNA"
snoRNA.types <- c("CDBox", "HAcaBox")
miRNA.ranges <- small.rna.ranges[ elementMetadata(small.rna.ranges)$type %in% miRNA.types ]
snoRNA.ranges <- small.rna.ranges[ elementMetadata(small.rna.ranges)$type %in% snoRNA.types ]

## chrM
chrM.range <- GRanges(seqnames="chrM", IRanges(start=1,end=seqlengths(read.ranges)[["chrM"]]), seqlengths=seqlengths(read.ranges))

## For each sequence that maps once to an anotated rRNA, tRNA, miRNA,
## or chrM, remove *all* mappings for that sequence.
forbidden.ranges <-
  Reduce(append,
         llply(list(rrna.ranges,
                    trna.ranges,
                    ## miRNA.ranges,
                    chrM.range),
               function(x) { elementMetadata(x) <- NULL; x }))

generate.null.ranges <- function(y) {
  x <- GRanges(seqnames=names(seqlengths(y)), ranges=IRanges(1,0), strand="*", seqlengths=seqlengths(y))
  for (i in names(elementMetadata(y))) {
    elementMetadata(x)[[i]] <- as(NA, class(as.vector(elementMetadata(y)[[i]])))
  }
  x
}

select.nearest <- function(x, y) {
  y <- append(y, generate.null.ranges(y))
  y[nearest(x,y)]
}

annotate.by.granges <- function(peaks, gr, annot.columns) {
  for (i in names(elementMetadata(gr))) {
    elementMetadata(gr)[[i]] <- as.vector(elementMetadata(gr)[[i]])
  }
  nearest.ranges <- select.nearest(peaks, gr)
  nearest.distance <- distance(peaks, nearest.ranges, ignore.strand=TRUE)
  in.ranges <- Rle(nearest.distance == 0)
  annot.data <- elementMetadata(nearest.ranges)[annot.columns]
  if (!is.null(names(annot.columns))) {
    names(annot.data) <- names(annot.columns)
  }
  DataFrame(overlap=in.ranges, distance=nearest.distance, annot.data)
}

forbidden.seqs <- unique(as.vector(elementMetadata(read.ranges)$seq[read.ranges %in% forbidden.ranges]))
forbidden.indices <- !is.na(findOverlaps(read.ranges, forbidden.ranges, select="first", ignore.strand=TRUE))
forbidden.read.ranges <- read.ranges[forbidden.indices]
read.ranges <- read.ranges[!forbidden.indices]

## Read the counts table
read.counts <- {
  x <- read.csv("./data/R21_R82_read_expr_matrix_no_blanks",
                stringsAsFactors=FALSE, sep="\t", header=TRUE, row.names=1)
  row.names(x) <- str_trim(row.names(x))
  names(x) <- str_replace(names(x), "_collapsed_sorted_with_zeros$", "")
  x <- x[row.names(x) %in% elementMetadata(read.ranges)$seq,]
  x <- x[order(names(x))]
  x
}

## Create per-sample bed files with score = count / multimap
sample.read.ranges <-
  llply(names(read.counts),
        function (sample) {
          x <- read.ranges
          ## Score = sample count / multimap
          elementMetadata(x)$score <-
            as.vector(read.counts[as.vector(elementMetadata(x)$seq), sample] / elementMetadata(x)$multimap)
          ## Eliminate reads with zero score
          x <- x[elementMetadata(x)$score > 0]
          x
        }, .parallel=TRUE)
names(sample.read.ranges) <- names(read.counts)

## Run blockbuster on each sample
## x <- sample.read.ranges[[1]][1:500]
## y <- run.blockbuster(x)
## z <- calculate.block.entropy(y$tags)
blockbuster.results <- llply(sample.read.ranges, function(x) run.blockbuster(unstranded(x), minBlockHeight=5), .parallel=TRUE)

## calculate.block.entropy(blockbuster.results[[1]]$tags[1:5])

## Calculate entropy of blocks
block.entropy <- llply(blockbuster.results, function(x) calculate.block.entropy(x$tags), .parallel=TRUE)

for (i in names(blockbuster.results)) {
  elementMetadata(blockbuster.results[[i]]$blocks)$entropy <- block.entropy[[i]][as.character(elementMetadata(blockbuster.results[[i]]$blocks)$blockID)]
}

## Plot entropy vs block length
x <- blockbuster.results[[1]]$blocks
y <- cbind(as.data.frame(elementMetadata(x)), width=width(x))
ggplot(y, aes(x=width, y=entropy, color=..density..)) + stat_density2d(geom="tile", aes(fill=..density..), contour=FALSE) + scale_fill_gradient(low="blue", high="yellow") + scale_color_gradient(low="blue", high="yellow")

## Compute nearest miRNA/snoRNA to each block and add annotation
block.annot <- llply(blockbuster.results, function(br) {
  annotate.by.granges(br$blocks, small.rna.ranges, c(nearest.ncRNA="name", ncRNA.type="type"))
}, .parallel=TRUE)
for (i in names(blockbuster.results)) {
  elementMetadata(blockbuster.results[[i]]$blocks)[names(block.annot[[i]])] <- block.annot[[i]]
}

## write output
saveRDS(blockbuster.results, "blockbuster_results.RDS")
block.tables <- llply(blockbuster.results, function(br) as(granges.to.dataframe(br$blocks, ignore.strand=TRUE, include.width=TRUE), "data.frame"))
write.xlsx.multisheet(block.tables, "blockbuster_results.xlsx", row.names=FALSE)
