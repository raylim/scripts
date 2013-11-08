#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("VennDiagram"));
suppressPackageStartupMessages(library("gsalib"));
suppressPackageStartupMessages(library("hwriter"));


options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outDir", default = ".", help = "Output dir"))

parser <- OptionParser(usage = "%prog [options] [hs_metrics.txt]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

arguments <- commandArgs(T);

if (length(arguments$args) < 1) {
    cat("Need input hs metrics file\n");
    print_help(parser);
    stop();
} else {
    f <- arguments$args[1];
}

grp <- gsa.read.gatkreport(f)


venn.plot <- venn.diagam(

