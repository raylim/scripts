#!/usr/bin/env Rscript
# convert a bam file to Rdata

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("Rsamtools"));

arguments <- commandArgs(trailingOnly = T)

if (length(arguments) < 1) {
    cat("Need bam to convert");
    stop();
}

for (bam in arguments) {
    rdata <- sub("\\.bam$", ".Rdata", bam, perl = T);
    bam <- scanBam(bam);
    save(bam, file = rdata)
}


