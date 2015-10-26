#!/usr/bin/env Rscript

# Script to check reproducibility of fRMA training process

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
load("consistency.rda")

## Select a random subset of 20 arrays from each tissue (Would rather
## use entire dataset but ENOMEM)
set.seed(1986)
norm.exprs <- lapply(names(vectors), function(ttype) {
    stab <- sample.tables[[ttype]] %>% sample_n(20)
    tsmsg("Reading 20 random arrays for ", ttype)
    affy <- ReadAffy(filenames=stab$Filename, sampleNames=rownames(stab))
    tsmsg("Normalizing with RMA for comparison")
    eset.rma <- rma(affy)
    rma.exprs <- eset.rma %>% exprs %>% as.vector
    tsmsg("Normalizing with 5 trains fRMA vector sets")
    esets.frma <- lapply(vectors[[ttype]], . %>% frma(affy, input.vecs=.))
    frma.exprs <- sapply(esets.frma, . %>% exprs %>% as.vector)
    data.frame(RMA=rma.exprs, fRMA=frma.exprs)
}) %>% setNames(names(vectors))

## Save because the above takes a while
save.image("consistency.rda")

dir.create("fRMA_consistency_results", FALSE)

## Compute M/A data for all pairwise comparisons
ma.data <- lapply(norm.exprs, function(normexprs) {
    normexprs %>% names %>% combn(2, simplify=FALSE) %>% ldply(. %>% {
        col1 <- .[1]
        col2 <- .[2]
        x1 <- normexprs[[col1]]
        x2 <- normexprs[[col2]]
        data.frame(Comparison = str_c(col1,".vs.",col2),
                   x1=x1, x2=x2,
                   M = x2 - x1,
                   A = (x1+x2)/2)
    })
})

for (ttype in names(norm.exprs)) {
    madata <- ma.data[[ttype]]

    pdf(str_c("fRMA_consistency_results/MA Plots for ", ttype, ".pdf"))
    ## MA plot for every pair of normalizations
    ddply(madata, .(Comparison), . %$% {
        smoothScatter(x=A, y=M, nbin=512, main=Comparison[1], xlab="A", ylab="M")
    })
    dev.off()

    ## M boxplots & violin plots
    pdf(str_c("fRMA_consistency_results/M Boxplots for ", ttype, ".pdf"),
        width=8, height=12)
    p <- ggplot(madata) + aes(x=Comparison, y=M) +
        scale_x_discrete(limits = rev(levels(madata$Comparison))) +
            coord_flip()
    print(p + geom_boxplot(notch=TRUE, outlier.shape = NA) +
          ggtitle(str_c("Boxplots of M value distributions for ", ttype)))
    print(p + geom_violin(scale="width") +
          ggtitle(str_c("Violin plots of M value distributions for ", ttype)))
    dev.off()
}
