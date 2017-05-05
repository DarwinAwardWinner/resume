#!/usr/bin/env Rscript

default.align.opts <- list(match=1, mismatch=3,
                           gapOpening=5, gapExtension=2)

parse_arguments <- function() {
    suppressMessages({
        library(optparse)
        library(parallel)
    })
    option_list <-
        list(make_option(c("-c", "--min-call"), type="integer", default=10, metavar="10",
                         help="Minimum perfect overlap required to call the presence of the subject (paired only). Imperfect overlap will need to be longer, based on specified mismatch and gap penalties."),
             make_option(c("-l", "--min-length"), type="integer", default=36, metavar="36",
                         help="Minimum length allowed after trimming a read. Any reads shorter than this after trimming will be discarded."),
             make_option(c("-i", "--interleaved"), action="store_true", default=FALSE,
                         help="Specify this option if you have paired-end sequences interleaved in a single FASTQ file. The default is to read paired-end sequences from a matched pair of files, and this option is ignored if two fastq files are provided. When you use this option, skip the \"READ2.fastq\" argument."),
             make_option(c("-o", "--read1-orientation"), type="character", default="in", metavar="in/out",
                         help="Orientation of read1. Can be either \"in\" or \"out\" (paired only). Note that Illumina reads are \"in\"."),
             make_option(c("-q", "--read2-orientation"), type="character", default="in", metavar="in/out",
                         help="Orientation of read2. Can be either \"in\" or \"out\" (paired only)"),
             make_option(c("-j", "--jobs"), type="integer",
                         default=parallel:::detectCores(),
                         metavar=as.character(parallel:::detectCores()),
                         help="Number of jobs to run in parallel for alignment. This should be autodetected by default."),
             make_option(c("-y", "--yield-size"), type="integer",
                         default=100000,
                         metavar="100000",
                         help="Number of reads to process at a time. Setting this higher will read more data into memory at once and result in faster runtime. Setting this lower will require less memory."),
             make_option(c("-m", "--match-bonus"), type="double",
                         default=default.align.opts$match,
                         metavar=as.character(default.align.opts$match),
                         help="Score bonus for a matching nucleotide"),
             make_option(c("-p", "--mismatch-penalty"), type="double",
                         default=default.align.opts$mismatch,
                         metavar=as.character(default.align.opts$mismatch),
                         help="Score penalty for a mismatched nucleotide (specify as a positive number)"),
             make_option(c("-g", "--gap-open-penalty"), type="double",
                         default=default.align.opts$gapOpening,
                         metavar=as.character(default.align.opts$gapOpening),
                         help="Score penalty for opening a gap in the alignment (specifiy as a positive number)"),
             make_option(c("-e", "--gap-extension-penalty"), type="double",
                         default=default.align.opts$match,
                         metavar=as.character(default.align.opts$gapExtension),
                         help="Score penalty for extending an alignment gap by two nucleotides (specify as a positive number)"),
             make_option(c("-s", "--single-read-mode"), action="store_true", default=FALSE,
                         help="Tell DeLoxer to run in single-end mode instead of paired-end mode. In this mode, the only a single input fastq file is provided, and only a single output file is created. No classification is performed, only trimming.  When you use this option, skip the \"READ2.fastq\" argument, and specify the full file name for OUTPUT_NAME instead of just the base name."))
    option_parser <- OptionParser(option_list=option_list,
                                  usage="%prog [options] adapter.fasta READ1.fastq READ2.fastq OUTPUT_NAME")
    opt <- parse_args(option_parser, positional_arguments=TRUE)
    return(opt)
}

## Call this here to handle --help quickly, before we waste 10 seconds
## loading all the libraries
invisible(parse_arguments())

print.option.list <- function(opt=parse_arguments()) {
    args <- opt$args
    opts <- opt$options
    message("Options:")
    foreach (o=opts, n=names(opts)) %do% {
        if (n != "help")
            message(" ", n, ": ", o)
    }
    message("Args: ", paste("\"", args, "\"", sep="", collapse=", "))
}

unimplemented <- function() stop("UNIMPLEMENTED")

## Timestampped message
tsmsg <- function(...) {
    message("# ", date(), ": ", ...)
}

tsmsg("Starting deloxer and loading required packages")

suppressMessages({
    library(ShortRead)
    library(optparse)
    library(foreach)
    library(iterators)
    library(itertools)
    library(doMC)
    registerDoMC()
    mcoptions <- list(preschedule=TRUE, set.seed=FALSE)
})

## Merge l1 and l2 by names
merge.lists <- function(l1, l2) {
    new.names <- setdiff(names(l2), names(l1))
    l1[new.names] <- l2[new.names]
    l1
}

## Return an object sans names
strip.names <- function(x) {
    names(x) <- NULL
    x
}

## Define some missing type coercions
setAs(from="ShortRead", to="DNAStringSet", def=function(from) sread(from))
setAs(from="PhredQuality", to="FastqQuality", def=function(from) FastqQuality(BStringSet(from)))
setAs(from="SolexaQuality", to="SFastqQuality", def=function(from) SFastqQuality(BStringSet(from)))
setAs(from="QualityScaledXStringSet", to="ShortReadQ", def=function(from) {
    q <- quality(from)
    new.quality.class <- switch(class(q),
                                SolexaQuality="SFastqQuality",
                                PhredQuality="FastqQuality",
                                stop("Unknown quality type: ", class(q)))
    q <- as(q, new.quality.class)
    ShortReadQ(sread=as(from, "DNAStringSet"),
               quality=q,
               id=BStringSet(names(from)))
})
## Override the provided method to keep the sequence names
setAs(from="ShortReadQ", to="QualityScaledDNAStringSet",
      def=function (from, to = "QualityScaledDNAStringSet", strict = TRUE) {
          q <- quality(from)
          new.quality.class <- switch(class(q),
                                      SFastqQuality="SolexaQuality",
                                      FastqQuality="PhredQuality",
                                      "XStringQuality")
          q <- as(q, new.quality.class)
          x <- QualityScaledDNAStringSet(sread(from), q)
          names(x) <- as.character(id(from))
          x
      })

## Define functions for reading fastq into standard Biostrings object
## and writing it back out. The standard functions readFastq and
## writeFastq operate on ShortRead objects. These simply wrap them in
## conversion to/from QualityScaledDNAStringSet.
read.QualityScaledDNAStringSet <- function(filepath, format = "fastq", ...) {
    switch(format,
           fastq=as(readFastq(filepath, withIds=TRUE, ...), "QualityScaledDNAStringSet"),
           ## Default
           stop("Unknown quality-scaled sequence format: ", format))
}

write.QualityScaledDNAStringSet <- function (x, filepath, append = FALSE, format = "fastq") {
    if(length(x) > 0) {
        sr <- as(x, "ShortReadQ")
        switch(format,
               fastq={
                   if (!append)
                       unlink(filepath);
                   writeFastq(object=sr,
                              file=filepath, mode=ifelse(append, "a", "w"))
               },
               ## Default
               stop("Unknown quality-scaled sequence format: ", format))
    } else {
        ## Zero-length sequence; just truncate/touch the file
        sink(file=filepath, append=append)
        sink()
    }
}

discard.short.reads <- function(reads, min.length=1) {
    kept.reads <- reads[width(reads) >= min.length]
    return(kept.reads)
}

## Takes a set of interleaved reads (or anything else) and
## de-interleaves them
deinterleave.pairs <- function(reads) {
    stopifnot(length(reads) %% 2 == 0)
    mask <- seq(from=1, to=length(reads), by=2)
    return(list(read1=reads[mask], read2=reads[-mask]))
}

.delox.trimmed.ranges <- function(subj, reads, min.length=36,
                                  include.scores=TRUE,
                                  include.deleted.ranges=TRUE,
                                  align.opts=list()) {

    align.opts <- merge.lists(align.opts, default.align.opts)

    aln <- list(forward=pairwiseAlignment(pattern=reads,
                subject=subj,
                type="overlap",
                substitutionMatrix=nucleotideSubstitutionMatrix(match = align.opts$match, mismatch = -align.opts$mismatch),
                gapOpening=-align.opts$gapOpening, gapExtension=-align.opts$gapExtension),
                revcomp=pairwiseAlignment(pattern=reads,
                subject=reverseComplement(DNAString(subj)),
                type="overlap",
                substitutionMatrix=nucleotideSubstitutionMatrix(match = align.opts$match, mismatch = -align.opts$mismatch),
                gapOpening=-align.opts$gapOpening, gapExtension=-align.opts$gapExtension))

    aln.scores <- Map(score, aln)
    aln.pat <- Map(pattern, aln)
    aln.ranges <- Map(function(x) IRanges(start=start(x), end=end(x)), aln.pat)
    aln.threebands <- Map(function (x) threebands(IRanges(start=1, end=width(reads)),
                                                  start=start(x), end=end(x)),
                          aln.ranges)

    ## For each read, decide whether the forward or reverse alignment
    ## was better.
    revcomp.better <- aln.scores$forward < aln.scores$revcomp

    ## For each read, take the threebands for the better alignment.
    best.threebands <- aln.threebands$forward
    for (band in names(best.threebands)) {
        best.threebands[[band]][revcomp.better] <- aln.threebands$revcomp[[band]][revcomp.better]
    }

    ## Use the left band if it is longer than either min.length or
    ## length of right band.
    use.right.band <- width(best.threebands$left) < pmin(min.length, width(best.threebands$right))
    ranges <- best.threebands$left
    ranges[use.right.band] <- best.threebands$right[use.right.band]

    ## Record which ranges are shorter than min.length
    too.short <- width(ranges) < min.length
    ## ranges[too.short] <- IRanges(start=1,end=0)

    ## Record what was trimmed off of each read (NOT what was kept!)
    trim <- factor(ifelse(use.right.band, "left", "right"), levels=c("right", "left", "all", "none"))
    ## If it's too short, then we trim "all", i.e. discard the whole
    ## read.
    trim[too.short] <- "all"
    ## If the read is not shorter after trimming, then nothing was
    ## actually trimmed.
    trim[width(ranges) == width(reads)] <- "none"

    emeta <- list()

    emeta$trim <- trim

    if (include.deleted.ranges) {
        deleted.start <- ifelse(too.short, 1,
                                ifelse(use.right.band,
                                       start(best.threebands$left),
                                       start(best.threebands$middle)))
        deleted.end <- ifelse(too.short, width(reads),
                              ifelse(use.right.band,
                                     end(best.threebands$middle),
                                     end(best.threebands$right)))
        emeta$deleted.range <- IRanges(deleted.start, deleted.end)
    }

    if (include.scores) {
        ## If requested, take the best score out of each pair of forward
        ## and reverse scores.
        scores <- ifelse(revcomp.better, aln.scores$revcomp, aln.scores$forward)
        emeta$score <- scores
    }

    mcols(ranges) <- DataFrame(emeta)

    return(ranges)
}

## Always call delox on the underlying DNAStringSet object when called
## on something more complicated.
suppressMessages({
    invisible(setMethod(".delox.trimmed.ranges", signature=c(reads="ShortRead"),
                        function (subj, reads, min.length, include.scores, include.deleted.ranges, align.opts) {
                            callGeneric(subj, as(reads, "DNAStringSet"), min.length, include.scores, include.deleted.ranges, align.opts)
                        }))
    invisible(setMethod(".delox.trimmed.ranges", signature=c(reads="QualityScaledDNAStringSet"),
                        function (subj, reads, min.length, include.scores, include.deleted.ranges, align.opts) {
                            callGeneric(subj, as(reads, "DNAStringSet"), min.length, include.scores, include.deleted.ranges, align.opts)
                        }))
    invisible(setMethod(".delox.trimmed.ranges", signature=c(reads="QualityScaledXStringSet"),
                        function (subj, reads, min.length, include.scores, include.deleted.ranges, align.opts) {
                            callGeneric(subj, as(reads, "XStringSet"), min.length, include.scores, include.deleted.ranges, align.opts)
                        }))
})

delox.single <- function(subj, reads , min.length=36,
                         include.scores=TRUE, align.opts=list()) {
    tsmsg("Saving read names")
    saved.names <- BStringSet(names(reads))
    reads <- strip.names(reads)
    invisible(gc())

    tsmsg("Doing alignments")
    nchunks <- min(getDoParWorkers(), ceiling(length(reads)/1000))
    deloxed.ranges <- foreach(reads=isplitVector(reads, chunks=nchunks), .combine=c) %dopar% {
        .delox.trimmed.ranges(reads=reads, subj=subj, min.length=min.length,
                              include.scores=include.scores,
                              include.deleted.ranges=FALSE,
                              align.opts=align.opts)
    }
    ## maybe.chunkapply(.delox.trimmed.ranges,
    ##                  VECTOR.ARGS=list(reads=reads),
    ##                  SCALAR.ARGS=list(subj=subj, min.length=min.length,
    ##                    include.scores=include.scores,
    ##                    include.deleted.ranges=FALSE,
    ##                    align.opts=align.opts),
    ##                  min.chunk.size=1000,
    ##                  MERGE=c)

    tsmsg("Trimming reads")
    trimmed.reads <- narrow(reads, start(deloxed.ranges), end(deloxed.ranges))

    tsmsg("Restoring read names")
    names(trimmed.reads) <- as.character(saved.names)

    tsmsg("Adding metadata")
    emeta <- list()
    if (include.scores) {
        emeta$score <- mcols(deloxed.ranges)$score
    }
    if (length(emeta) > 0) {
        mcols(trimmed.reads) <- DataFrame(emeta)
    }

    return(discard.short.reads(trimmed.reads, min.length))
}

delox.paired <- function(subj, read1, read2,
                         min.call=10, min.length=36,
                         include.scores=TRUE, align.opts=list()) {
    align.opts <- merge.lists(align.opts, default.align.opts)

    tsmsg("Checking read counts")
    stopifnot(length(read1) == length(read2))

    tsmsg("Listing reads")
    original.reads <- list(read1=read1,
                           read2=read2)
    rm(read1, read2)

    tsmsg("Saving read names")
    read.names <- foreach(r=original.reads) %do% BStringSet(names(r))
    names(read.names) <- names(original.reads)
    original.reads <- Map(strip.names, original.reads)
    invisible(gc())

    tsmsg("Doing alignments")
    deloxed.ranges <- lapply(original.reads, function(x) {
        nchunks <- min(getDoParWorkers(), ceiling(length(x)/1000))
        foreach(reads=isplitVector(x, chunks=nchunks), .combine=c) %dopar% {
            .delox.trimmed.ranges(reads=reads, subj=subj, min.length=min.length,
                                  include.scores=include.scores,
                                  include.deleted.ranges=FALSE,
                                  align.opts=align.opts)
        }

        ## maybe.chunkapply(.delox.trimmed.ranges,
        ##                  VECTOR.ARGS=list(reads=strip.names(x)),
        ##                  SCALAR.ARGS=list(subj=subj,
        ##                    min.length=min.length,
        ##                    include.scores=TRUE,
        ##                    include.deleted.ranges=TRUE,
        ##                    align.opts=align.opts),
        ##                  MERGE=c,
        ##                  min.chunk.size=1000)
    })

    tsmsg("Extracting metadata")
    delox.meta <- lapply(deloxed.ranges, mcols)

    ## Decide whether enough was trimmed on the inside (right end) of
    ## either read to call it a mate-pair.
    tsmsg("Calculating inside trim score")
    inside.trim.score <- Reduce(pmax,
                                lapply(delox.meta,
                                       function(x) ifelse(x$trim == "right", x$score, 0)))

    ## Decide whether enough was trimmed on the outside (left end) of
    ## either read to call it a non-mate-pair.
    tsmsg("Calculating outside trim score")
    outside.trim.score <- Reduce(pmax,
                                 lapply(delox.meta,
                                        function(x) ifelse(x$trim == "left", x$score, 0)))

    tsmsg("Calling presence of subject")
    calls <- list(inside=inside.trim.score >= min.call * align.opts$match,
                  outside=outside.trim.score >= min.call * align.opts$match)

    tsmsg("Categorizing reads")
    category <- factor(rep(NA, length(original.reads$read1)), levels=c("mate", "non-mate", "negative", "unpaired", "discard"))
    category[calls$inside] <- "mate"
    category[calls$outside] <- "non-mate"
    ## If they're either both true or both false, then it's ambiguous
    category[calls$inside == calls$outside] <- "negative"
    ## All categories should be filled in now
    stopifnot(all(!is.na(category)))

    too.short <- lapply(deloxed.ranges, function(x) width(x) < min.length)
    ## If either read in a pair is too short, then its partner is no
    ## longer paired at all.
    one.too.short <- Reduce(`|`, too.short)
    category[one.too.short] <- "unpaired"
    ## If both reads in a pair are too short, then the entire pair is
    ## discarded. This is highly unlikely, since Cre-Lox should not
    ## appear in the middle of both sequences.
    both.too.short <- Reduce(`&`, too.short)
    category[both.too.short] <- "discard"

    tsmsg("Trimming reads and restoring read names")
    trimmed.reads <- lapply(names(original.reads), function(x) {
        trimmed <- narrow(original.reads[[x]],
                          start=start(deloxed.ranges[[x]]),
                          end=end(deloxed.ranges[[x]]))
        names(trimmed) <- as.character(read.names[[x]])
        trimmed
    })
    names(trimmed.reads) <- names(original.reads)

    tsmsg("Assembling metadata")
    foreach (r=names(trimmed.reads)) %do% {
        emeta <- list()
        emeta$category <- category
        emeta$category[too.short[[r]]] <- "discard"
        if (include.scores) {
            emeta$score <- delox.meta[[r]]$score
        }
        mcols(trimmed.reads[[r]]) <- DataFrame(emeta)
    }

    return(trimmed.reads)
}

## Wrapper for both single and paired as appropriate
delox <- function(subj, read1, read2=NULL,
                  min.call=10, min.length=36,
                  interleaved=FALSE,
                  read1.orientation=c("in", "out")[1],
                  read2.orientation=c("in", "out")[1],
                  align.opts=list()) {
    if (is.null(read2)) {
        if (interleaved) {
            x <- deinterleave.pairs(read1)
            read1 <- x$read1
            read2 <- x$read2
        } else {
            tsmsg("Doing single-read delox")
            return(delox.single(subj=subj, reads=read1, min.length=min.length, align.opts=align.opts))
        }
    }

    ## Make sure both reads are oriented "in" before calling
    tsmsg("Ensuring correct read orientation")
    if (tolower(read1.orientation) == "out") {
        read1 <- reverseComplement(read1)
    }
    if (!is.null(read2) && tolower(read2.orientation) == "out") {
        read2 <- reverseComplement(read2)
    }

    tsmsg("Doing paired-end delox")
    deloxed.reads <- delox.paired(subj, read1, read2,
                                  min.call=min.call, min.length=min.length,
                                  align.opts=align.opts)

    ## If reads started "out", put them back that way before returning
    tsmsg("Restoring original read orientation")
    if (tolower(read1.orientation) == "out") {
        deloxed.reads$read1 <- reverseComplement(deloxed.reads$read1)
    }
    if (tolower(read2.orientation) == "out") {
        deloxed.reads$read2 <- reverseComplement(deloxed.reads$read2)
    }

    return(deloxed.reads)
}

## ## Hack to work around a bug in BioConductor that prevents subsetting
## ## of named XStringSet objects. Apparently, since DeLoxer was first
## ## published, the BioConductor devs broke the XStringSet subsetting
## ## code so that it can no longer handle XStringSets with names. The
## ## code below strips the names from the XStringSet, then calls the old
## ## code to subset the nameless object while subsetting the names
## ## separately, then finally puts the names back on and returns the
## ## result.
## old.XStringSet.subset.method <- selectMethod("[", "XStringSet")
## invisible(setMethod("[", signature="XStringSet", definition=function(x, i, j, ..., drop=TRUE) {
##     ## Save the names into a seaprate variable
##     xnames <- names(x)
##     ## Do the old behavior, which works on unnamed objects
##     x <- old.XStringSet.subset.method(unname(x), i, j, ..., drop=drop)
##     ## Put the names back on and return
##     setNames(x, xnames[i])
## }))

save.deloxed.pairs.as.fastq <- function(read1, read2, output.base,
                                        mate.ext="matepaired",
                                        nonmate.ext="paired",
                                        negative.ext="negative",
                                        unpaired.ext="unpaired",
                                        append=FALSE) {

    extension <- c(mate=mate.ext,
                   `non-mate`=nonmate.ext,
                   negative=negative.ext,
                   unpaired=unpaired.ext)

    ## ## Make sure that read1 and read2 are a match for each other
    ## stopifnot(identical(as.character(mcols(read1)$category),
    ##                     as.character(mcols(read2)$category)))

    ## ## Discard the shorter read on "unpaired"
    ## read1.shorter <- width(read1) < width(read2)
    ## mcols(read1)$category[mcols(read1)$category == "unpaired" & read1.shorter] <- NA
    ## mcols(read2)$category[mcols(read2)$category == "unpaired" & !read1.shorter] <- NA

    filename.template <- "%s_read%s.%s.fastq"

    for (cat in names(extension)) {
        read1.for.category <- read1[mcols(read1)$category == cat]
        read1.file.for.category <- sprintf(filename.template, output.base, 1, extension[[cat]])
        tsmsg("Writing ", read1.file.for.category)
        write.QualityScaledDNAStringSet(read1.for.category,
                                        file=read1.file.for.category,
                                        append=append)

        read2.for.category <- read2[mcols(read2)$category == cat]
        read2.file.for.category <- sprintf(filename.template, output.base, 2, extension[[cat]])
        tsmsg("Writing ", read2.file.for.category)
        write.QualityScaledDNAStringSet(read2.for.category,
                                        file=read2.file.for.category,
                                        append=append)
    }

    return(TRUE)
}

get.category.counts <- function(deloxed.pairs) {
    r1cat <- mcols(deloxed.pairs$read1)$category
    r2cat <- mcols(deloxed.pairs$read2)$category
    x <- table(r1cat)[c("mate", "non-mate", "negative")]
    x["r1.single"] <- sum(r1cat == "unpaired")
    x["r2.single"] <- sum(r2cat == "unpaired")
    x["discard"] <- length(r1cat) - sum(x)
    x
}

mcparallel.quiet <- function(expr, ...) {
    parallel:::mcparallel(suppressMessages(expr), ...)
}

print.stats <- function(category.counts) {
    category.pct <- setNames(sprintf("%.3g%%", category.counts / sum(category.counts) * 100),
                             names(category.counts))
    x <- rbind(Counts=category.counts, Fractions=category.pct)
    names(dimnames(x)) <- c("Stat", "Category")
    print(x, quote=FALSE, justify="right")
}

main <- function() {
    opt <- parse_arguments()
    print.option.list(opt)
    args <- opt$args
    opts <- opt$options

    if (!(tolower(opts[["read1-orientation"]]) %in% c("in", "out") &&
          tolower(opts[["read2-orientation"]]) %in% c("in", "out") )) {
        stop("Valid orientations are \"in\" and \"out\"")
    }

    align.opts <- list(match = opts[["match-bonus"]],
                       mismatch = opts[["mismatch-penalty"]],
                       gapOpening = opts[["gap-open-penalty"]],
                       gapExtension = opts[["gap-extension-penalty"]])

    stopifnot(opts$`min-call` >= 1 &&
              opts$`min-length` >= 0 &&
              opts$`jobs` >= 0)

    ## Set jobs if requested
    if (opts$jobs > 0) {
        options(cores=opts$jobs)
    }
    tsmsg("Using ", getDoParWorkers(), " cores.")

    paired <- !opts[["single-read-mode"]]
    interleaved <- opts[["interleaved"]]

    if (!paired && interleaved) {
        stop("ERROR: You cannot specify both --interleaved and --single-read-mode")
    } else if (!paired) {
        if (length(args) != 3) {
            stop("DeLoxer in single-read mode requires exactly 3 arguments")
        }
        subject.file <- args[[1]]
        read1.file <- args[[2]]
        read2.file <- NULL
        output.file <- args[[3]]
    } else if (interleaved) {
        if (length(args) != 3) {
            stop("DeLoxer interleaved input mode requires exactly 3 arguments")
        }
        subject.file <- args[[1]]
        read1.file <- args[[2]]
        read2.file <- NULL
        output.basename <- args[[3]]
    } else {
        if (length(args) != 4) {
            stop("DeLoxer requires exactly 4 arguments")
        }
        subject.file <- args[[1]]
        read1.file <- args[[2]]
        read2.file <- args[[3]]
        output.basename <- args[[4]]
    }

    subj <- readDNAStringSet(subject.file, format="fasta", nrec=1)[[1]]

    yieldSize <- opts[["yield-size"]]
    if (paired) {
        tsmsg("Deloxing and classifying paired sequences")
        read1.stream <- FastqStreamer(read1.file, n=yieldSize)
        read2.stream <- if (!interleaved) FastqStreamer(read2.file, n=yieldSize)
        process.chunk <- function(fq1, fq2, append) {
            if (length(fq1) < 1)
                return(TRUE)
            if (interleaved) {
                stopifnot(is.null(fq2))
                deint <- deinterleave.pairs(fq1)
                fq1 <- deint[[1]]
                fq2 <- deint[[2]]
            } else {
                if (length(fq1) != length(fq2))
                    stop("Both input files must have equal numbers of reads")
            }
            read1 <- as(fq1, "QualityScaledDNAStringSet")
            read2 <- as(fq2, "QualityScaledDNAStringSet")
            deloxed.pairs <-
                delox(subj, read1, read2,
                      min.call=opts[["min-call"]],
                      interleaved=interleaved,
                      read1.orientation=opts[["read1-orientation"]],
                      read2.orientation=opts[["read2-orientation"]],
                      align.opts=align.opts)
            save.deloxed.pairs.as.fastq(deloxed.pairs$read1, deloxed.pairs$read2, output.basename, append=append)

            ret <- get.category.counts(deloxed.pairs)
            return(ret)
        }
        fq1 <- yield(read1.stream)
        fq2 <- if (!interleaved) yield(read2.stream)
        if (length(fq1) == 0)
            warning("No reads were read from the input file.")
        proc <- mcparallel.quiet(process.chunk(fq1, fq2, append=FALSE))
        reads.processed <- length(fq1) / ifelse(interleaved, 2, 1)
        category.stats <-
            category.counts <- NULL
        while (length(fq1 <- yield(read1.stream))) {
            if (!interleaved)
                fq2 <- yield(read2.stream)
            prev.result <- mccollect(proc)[[1]]
            if (is(prev.result, "try-error")) {
                tsmsg("Encountered error in deloxing subprocess:")
                stop(attr(prev.result, "condition"))
            }
            if (is.null(category.counts)) {
                category.counts <- prev.result
            } else {
                category.counts <- category.counts + prev.result
            }
            tsmsg("Category stats after processing ", reads.processed, " reads:")
            ## category.pct <- setNames(sprintf("%.3g%%", category.counts / sum(category.counts) * 100),
            ##                          names(category.counts))
            print.stats(category.counts)
            proc <- mcparallel.quiet(process.chunk(fq1, fq2, append=TRUE))
            reads.processed <- reads.processed + length(fq1) / ifelse(interleaved, 2, 1)
        }
        close(read1.stream)
        if (!interleaved) close(read2.stream)
        prev.result <- mccollect(proc)[[1]]
        if (is.null(category.counts)) {
            category.counts <- prev.result
        } else {
            category.counts <- category.counts + prev.result
        }
        if (is(prev.result, "try-error")) {
            tsmsg("Encountered error in deloxing subprocess:")
            stop(attr(prev.result, "condition"))
            stop("Encountered error in deloxing")
        }
        tsmsg("Final category stats after processing ", reads.processed, " reads:")
        print.stats(category.counts)
    } else {
        tsmsg("Deloxing single sequences")
        read1.stream <- FastqStreamer(read1.file, n=yieldSize)
        process.chunk <- function(fq, append) {
            if (length(fq) < 1)
                return(TRUE)
            reads <- as(fq, "QualityScaledDNAStringSet")
            deloxed.reads <-
                delox(subj, reads, NULL,
                      min.call=opts[["min-call"]],
                      interleaved=interleaved,
                      read1.orientation=opts[["read1-orientation"]],
                      read2.orientation=opts[["read2-orientation"]],
                      align.opts=align.opts)
            write.QualityScaledDNAStringSet(deloxed.reads, output.file, append=append)
            return(TRUE)
        }
        ## First chunk is processed with append=FALSE to start the file
        fq <- yield(read1.stream)
        if (length(fq) == 0)
            warning("No reads were read from the input file.")
        proc <- mcparallel.quiet(suppressMessages(process.chunk(fq, append=FALSE)))
        reads.processed <- length(fq)
        while (length(fq <- yield(read1.stream))) {
            prev.result <- mccollect(proc)[[1]]
            if (is(prev.result, "try-error")) {
                tsmsg("Encountered error in deloxing subprocess:")
                stop(attr(prev.result, "condition"))
                stop("Encountered error in deloxing")
            }
            tsmsg("Processed ", reads.processed, " reads")
            proc <- mcparallel.quiet(suppressMessages(process.chunk(fq, append=TRUE)))
            reads.processed <- reads.processed + length(fq)
        }
        close(read1.stream)
        prev.result <- mccollect(proc)[[1]]
        if (is(prev.result, "try-error")) {
            tsmsg("Encountered error in deloxing subprocess:")
            stop(attr(prev.result, "condition"))
            stop("Encountered error in deloxing")
        }
        tsmsg("Processed ", reads.processed, " reads")
    }
    tsmsg("Finished successful run")
}

main()
