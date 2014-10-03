#!/usr/bin/env Rscript
# run bicseq on tumor normal samples

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("BICseq"));

options(error = quote(dump.frames("testdump", TRUE)))

optList <- list(
                make_option("--genome", default = 'hg19', help = "genome build [default %default]"),
                make_option("--includeY", default = F, action = 'store_true', help = "include Y chromosome [default %default]"),
                make_option("--lambda", default = 10, help = "lambda [default %default]"),
                make_option("--bin", default = 100, help = "bin size [default %default]"),
                make_option("--winSize", default = 200, help = "window size [default %default]"),
                make_option("--quant", default = 0.95, help = "probability of the read count quantile [default %default]"),
                make_option("--mult", default = 1, help = "a genomic position s is  considered as an outlier if it has more than mult*qunatile number of aligned reads, where quantile is the quant^th quantile of the read counts calculated from the genomic window of s [default %default]"),
                make_option("--outFile", default = NULL, help = "output file [default %default]"))

parser <- OptionParser(usage = "%prog [tumor.bam] [normal.bam]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    cat("Need bam files\n");
    print_help(parser);
    stop();
}

tumorBam <- arguments$args[1]
normalBam <- arguments$args[2]

if (opt$genome == "hg19") {
    seqnames <- c(1:22, "X")
} else if (opt$genome == "mm10") {
    seqnames <- c(1:19, "X")
} else {
    cat("unsupported genome\n");
    stop();
}

if (opt$includeY) {
    seqnames <- c(seqnames, "Y")
}

bicseq <- BICseq(sample = tumorBam, reference = normalBam, seqNames = seqnames)
segs <- getBICseg(object = bicseq, bin = opt$bin, lambda = opt$lambda, winSize = opt$winSize, quant = opt$quant, mult = opt$mult)



