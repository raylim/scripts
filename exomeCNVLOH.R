#!/usr/bin/env Rscript
# generates amplicon coverage report (hg19)

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("ExomeCNV"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outPrefix", default = NULL, type = "character", action = "store", help ="Output prefix (required)"),
                make_option("--tumor", default = NULL, type = "character", action = "store", help ="tumor BAF file (required)"),
                make_option("--normal", default = NULL, type = "character", action = "store", help ="normal BAF file (not required)"),
                make_option("--lohMethod", default = "two.sample.fisher", type = "character", action = "store", help ="LoH calling method (default: %default)"),
                make_option("--cbsLohMethod", default = "variance.f", type = "character", action = "store", help ="CBS LoH calling method (default: %default)"),
                make_option("--alpha", default = 0.05, action = "store", help ="LoH alpha (default: %default)"))

parser <- OptionParser(usage = "%prog [options] ", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$tumor)) {
    cat("Need tumor BAF file\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n\n");
    print_help(parser);
    stop();
}

tumorName <- sub('\\..*', '', sub('.*/', '', opt$tumor))

if (!is.null(opt$normal)) {
    normalName <- sub('\\..*', '', sub('.*/', '', opt$normal))
    cat("Reading", normalName, ":", opt$normal, "\n")
    normal <- read.delim(opt$normal, header = T)
} else {
    normal <- NULL
    opt$lohMethod <- "only.tumor"
    opt$cbsLohMethod <- "only.tumor"
}

cat("Reading", tumorName, ":", opt$tumor, "\n")
tumor <- read.delim(opt$tumor, header = T)


if (!any(grepl('chr', normal$chr))) {
    if (!is.null(opt$normal)) {
        normal$chr <- sub('^', 'chr', normal$chr)
    }
    tumor$chr <- sub('^', 'chr', tumor$chr)
}


if (!is.null(opt$normal)) {
    x <- apply(normal, 1, function(x) any(is.na(x))) | apply(tumor, 1, function(x) any(is.na(x)))
    normal <- normal[!x, ]
}
tumor <- tumor[!x, ]


cat("Analyzing LOH using", opt$lohMethod, "\n")
cat("alpha:", opt$alpha, "\n")
eLOH = LOH.analyze(tumor = normal, tumor = tumor, alpha = opt$alpha, method = opt$lohMethod)

cat("Merging segments\n")
loh = multi.LOH.analyze(normal = normal, tumor = tumor, all.loh.ls = list(eLOH), test.alpha = 0.001, method = opt$cbsLohMethod, sdundo = c(0,0), alpha = c(0.05,0.01))

prefix <- paste(opt$outPrefix, sep = "")
cat("Writing output (prefix: ", opt$outPrefix, ")\n", sep = "")
write.loh.output(loh, opt$outPrefix)

fn <- paste(opt$outPrefix, ".loh.png", sep = "")
cat("Plotting to", fn, "\n")
png(filename = fn, res = 70, width = 2000,
    height = 1200, pointsize = 16, type = "cairo-png")
do.plot.loh(loh, normal = normal, tumor = tumor, method = opt$cbsLohMethod, plot.style = "baf")
dev.off()
