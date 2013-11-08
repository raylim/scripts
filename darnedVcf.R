#!/usr/bin/env Rscript
# Read a variant table and check darned database

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("VariantAnnotation"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--genome", default = 'hg19', help = "genome build [default %default]"),
                make_option("--darnedDb", default = NULL, help = "darned database"),
                make_option("--outFile", default = stdout(), help = "vcf output file [default %default]"))

parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$darnedDb)) {
    cat("Need DARNED RNA-editing database\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    #cat("Reading from stdin...\n");
    tmp <- tempfile()
    f <- file('stdin')
    open(f)
    x <- readLines(f)
    close(f)
    write(x, file = tmp)
    vcf <- readVcf(tmp, genome = opt$genome)
} else {
    fn <- arguments$args[1];
    vcf <- readVcf(fn, genome = opt$genome)
}


darnedDb <- read.table(opt$darnedDb, sep = '\t', header = T, as.is = T, quote = '', comment.char = '')

al <- sapply(rowData(vcf)$ALT, length)
X <- data.frame(rowId = 1:nrow(vcf), seq = as.character(seqnames(rowData(vcf))), pos = as.numeric(start(rowData(vcf))))
X <- as.data.frame(apply(X, 2, rep, times = al))
X$pos <- as.integer(as.character(X$pos))
X$alt <- as.character(unlist(rowData(vcf)$ALT))

XX <- merge(X, darnedDb, by.x = c('seq', 'pos', 'alt'), by.y = c('chrom', 'coordinate', 'inchr'), all.x = T)

