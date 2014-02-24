#!/usr/bin/env Rscript
# applies absCNseq

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("absCNseq"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outPrefix", default = NULL, type = "character", action = "store", help ="Output prefix (required)"))

parser <- OptionParser(usage = "%prog [options] varscan.seg.file snv.table", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(argumetns$args) != 2) {
    cat("Need varscan seg file and SNV table\n\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n\n");
    print_help(parser);
    stop();
}

segFile <- arguments$args[1]
snvFile <- arguments$args[2]

tumorName <- sub('\\..*', '', sub('.*/', '', sub('_.*', '', segFile)))

segData <- read.table(segFile, sep = '\t', header = T)
segRle <- rle(segData$Segmented)
segData <- transform(segData, segId = as.factor(rep(1:length(segRle$values), segRle$lengths)))
segData <- transform(segData, length = End - Start)


chrom <- tapply(segData$Chrom, segData$segId, function(x) x[1])
start <- tapply(segData$Start, segData$segId, function(x) x[1])
end <- tapply(segData$End, segData$segId, function(x) x[length(x)])
effSegLen <- tapply(segData$length, segData$segId, sum)
normRatio <- tapply(segData$Segmented, segData$segId, function(x) x[1])

absSegData <- data.frame(chrom = chrom, loc.start = start, loc.end = end, eff.seg.len = effSegLen, normalized.ratio = normRatio)

snvData <- read.table(snvFile, sep = '\t', header = T, comment.char = '')
af <- paste(tumorName
absSnvData <- with(snvData, data.frame(chrom = X.CHROM, position = POS, tumor_var_freq = snvData[,af])

