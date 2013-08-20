#!/usr/bin/env Rscript
# segment copy number data and generate plot

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("DNAcopy"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--prefix", default = NULL, type = "character", action = "store", help ="Output prefix (required)"))

parser <- OptionParser(usage = "%prog [options] copynumFile", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 1) {
    cat("Need copy number file\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$prefix)) {
    cat("Need output prefix\n\n");
    print_help(parser);
    stop();
} else {
    cnFile <- arguments$args[1];
}

cat("Reading", cnFile, " ... ")
cn <- read.table(cnFile, header = T)
cat("done\n")
CNA.object <- CNA(genomdat = cn$adjusted_log_ratio, chrom = cn$chrom, maploc = cn$chr_start, data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
cat("Segmenting ... ")
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
cat("Done\n")

fn <- paste(opt$prefix, '.segment.txt', sep = '')
cat("Writing to", fn, " ... ")
write.table(segs$output[ ,-1], file = fn , row.names=F, quote=F, sep="\t")
cat("done\n")

fn <- paste(opt$prefix, '.cnv_full.png', sep = '')
cat("Plotting to", fn, " ... ")
png(filename = fn, res = 70, width = 1000,
    height = 800, pointsize = 16, type = "cairo-png")
plot(segs)
null <- dev.off()
cat("done\n")

fn <- paste(opt$prefix, '.cnv_by_chr.png', sep = '')
cat("Plotting to", fn, " ... ")
png(filename = fn, res = 70, width = 2000,
    height = 2400, pointsize = 16, type = "cairo-png")
plot(segs, plot.type = "samplebychrom")
null <- dev.off()
cat("done\n")

cat("Finished\n")
