#!/usr/bin/env Rscript
# extract allele read counts for titan from bam file using specified vcf

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("TitanCNA"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
        make_option("--out", default = NULL, type = "character", action = "store", help ="Output file (required)"),
        make_option("--vcf", default = NULL, type = "character", action = "store", help ="positions in vcf format (required)"))

parser <- OptionParser(usage = "%prog [options] [bam file]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need input bam\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$out)) {
    cat("Need output file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$vcf)) {
    cat("Need vcf file\n\n")
    print_help(parser);
    stop();
}

fn <- arguments$args[1]
if (!file.exists(fn)) {
    cat("bam file missing\n\n")
    print_help(parser);
    stop();
}
ifn <- paste(fn, 'bai', sep = '.')
if (!file.exists(ifn)) {
    cat("bam index missing\n\n")
    print_help(parser);
    stop();
}

extractAlleleReadCounts(fn, ifn, opt$vcf, opt$out)
