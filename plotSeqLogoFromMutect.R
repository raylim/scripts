#!/usr/bin/env Rscript
# plot berry logos using mutect input

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("berrylogo"));

options(warn = -1, error = traceback)

optList <- list(
                make_option("--out", default = NULL, type = "character", action = "store", help = "Output .png (required)"))

parser <- OptionParser(usage = "%prog [options]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 1) {
    cat("Need mutect results table\n");
    print_help(parser);
    stop();
} else if (is.null(opt$out)) {
    cat("Need output file\n\n");
    print_help(parser);
    stop();
} else {
    mutectFile <- arguments$args[1];
}

d <- read.table(mutectFile, header=T, as.is = T);
d <- subset(d, judgement == "KEEP" & dbsnp_site == "NOVEL")

context <- t(sapply(d$context, function(x) unlist(strsplit(x, ''))))
contextTable <- do.call('cbind', apply(context, 2, function(x) table(factor(x))))
contextTable[,4] <- table(d$alt_allele)
pwm <- apply(contextTable, 2, function(x) x / sum(x))

png(opt$out, height = 800, width = 800)
berrylogo(pwm)
dev.off()
