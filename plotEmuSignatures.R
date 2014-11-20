#!/usr/bin/env Rscript
# plots emu signature

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("VariantAnnotation"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("org.Hs.eg.db"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outPrefix", default = NULL, help = "Output prefix (required)"),
                make_option("--samples", default = NULL, help = "samples file; labels axes if specified"),
                make_option("--inPrefix", default = NULL, help = "input prefix (required)"))

parser <- OptionParser(usage = "%prog [options]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n");
    print_help(parser);
    stop();
} else if (is.null(opt$inPrefix)) {
    cat("Need input prefix\n");
    print_help(parser);
    stop();
} else {
    files <- arguments$args;
}

glob <- paste(opt$inPrefix, '_*_ml_spectra.txt', sep = '')
spectraFiles <- Sys.glob(glob)

glob <- paste(opt$inPrefix, '_*_map_activities.txt', sep = '')
activityFiles <- Sys.glob(glob)

glob <- paste(opt$inPrefix, '_*_assigned.txt', sep = '')
assignedFiles <- Sys.glob(glob)
