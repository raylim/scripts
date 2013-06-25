#!/usr/bin/env Rscript
# generates amplicon coverage report (hg19)

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("ExomeCNV"));

options(warn = -1, error = traceback)

optList <- list(
                make_option("--outDir", default = NULL, type = "character", action = "store", help ="Output directory (required)"),
                make_option("--lohMethod", default = "two.sample.fisher", type = "character", action = "store", help ="LoH calling method (default: %default)"),
                make_option("--alpha", default = 0.05, action = "store", help ="LoH alpha (default: %default)"))

parser <- OptionParser(usage = "%prog [options] tumorBAF normalBAF", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) != 2) {
    cat("Need BAF files for tumor and normal\n\n");
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

normalName <- sub('\\..*', '', sub('.*/', '', normalFile))
tumorName <- sub('\\..*', '', sub('.*/', '', tumorFile))

cat("Reading", normalName, ":", normalFile, "\n")
normal = read.delim(normalFile, header = T)

cat("Reading", tumorName, ":", tumorFile, "\n")
tumor = read.delim(tumorFile, header = T)

cat("Analyzing LOH using", opt$lohMethod, "\n")
cat("alpha:", opt$alpha, "\n")
eLOH = LOH.analyze(normal, tumor, alpha = opt$alpha, method = opt$lohMethod)

cat("Merging segments\n")
loh = multi.LOH.analyze(normal, tumor, all.loh.ls=list(eLOH), test.alpha=0.001, method="variance.f", sdundo=c(0,0), alpha=c(0.05,0.01))

prefix <- paste(opt$outDir, "/", tumorName, "_", normalName, sep = "")
cat("Writing output (prefix: ", prefix, ")\n", sep = "")
write.loh.output(eLOH, prefix)

fn <- paste(prefix, ".loh.png", sep = "")
cat("Plotting to", fn, "\n")
png(filename = fn, res = 70, width = 2000,
    height = 1200, pointsize = 16, type = "cairo-png")
do.plot.loh(loh, normal, tumor, "two.sample.fisher", plot.style="baf")
dev.off()
