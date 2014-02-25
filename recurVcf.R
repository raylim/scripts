#!/usr/bin/env Rscript
# output positions that occur in more than one vcf file

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("plyr"));
suppressPackageStartupMessages(library("VariantAnnotation"));

#options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))
#options(error = recover)
options(error = quote(dump.frames("testdump", TRUE)))

optList <- list(
                make_option("--genome", default = 'hg19', help = "genome build [default %default]"),
                make_option("--outFile", default = NULL, help = "output file [default %default]"))

parser <- OptionParser(usage = "%prog vcf.files", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) <= 1) {
    cat("Need vcf files\n");
    print_help(parser);
    stop();
}

files <- arguments$args;

vcfs <- list()

for (f in files) {
    vcfs <- append(vcfs, readVcf(f, genome = opt$genome))
}

all <- do.call('rbind', lapply(vcfs, function(x) as.data.frame(subset(rowData(vcf), FILTER == "PASS"))))
all <- all[, c("seqnames", "start", "end")]
cnt <- ddply(all, .(seqnames, start, end), nrow)

cnt <- subset(cnt, V1 > 1)

write.table(cnt, file = opt$outFile, sep = '\t', quote = F)
