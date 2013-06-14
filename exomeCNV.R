#!/usr/bin/env Rscript
# generates amplicon coverage report (hg19)

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("doMC"));
suppressPackageStartupMessages(library("multicore"));
suppressPackageStartupMessages(library("ExomeCNV"));

registerDoMC(8)
options(warn = -1, error = traceback)

optList <- list(
                make_option("--outDir", default = NULL, type = "character", action = "store", help ="Output directory (required)"))

parser <- OptionParser(usage = "%prog [options] tumorGatkCoverageFile normalGatkCoverageFile", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 2) {
    cat("Need GATK coverage file (sample_interval_summary) for tumor and normal\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outDir)) {
    cat("Need output dir\n\n");
    print_help(parser);
    stop();
} else {
    tumorFile <- arguments$args[1];
    normalFile <- arguments$args[2];
}

normal <- read.coverage.gatk(normalFile)
tumor <- read.coverage.gatk(tumorFile)

chrs <- levels(normal$chr)

# analysis overview
#1. Calculate log coverage ratio between case and control
logR <- calculate.logR(normal, tumor)
eCNV <- foreach(i = 1:length(chrs), .combine = rbind) %dopar% {
    idx = (normal$chr == chrs[i])
    classify.eCNV(normal = normal[idx,], tumor = tumor[idx,], logR = logR[idx], min.spec = 0.9999, min.sens = 0.9999, option="spec", c = 0.5, l = 70)
}

cnv = multi.CNV.analyze(normal, tumor, logR=logR, all.cnv.ls=list(eCNV), coverage.cutoff=5, min.spec=0.99, min.sens=0.99, option="auc", c=0.5)
#2. Call CNV/LOH for each exon individually
#3. Combine exonic CNV/LOH into segments using Circular Binary Segmentation (CBS)

