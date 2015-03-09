#!/usr/bin/env Rscript
# parse gene list
# usage: ./addGeneListAnnotationToVcf.pl 
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("rtracklayer"));
suppressPackageStartupMessages(library("VariantAnnotation"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
        make_option("--genome", default = 'hg19', help = "genome build [default %default]"),
        make_option("--geneBed", default = NULL, type = "character", action = "store", help ="input gene bed (required)"),
        make_option("--name", default = 'genelist', type = "character", action = "store", help ="annotation name (default %default)"),
        make_option("--outFile", default = NULL, type = "character", action = "store", help ="output file (required)"))

parser <- OptionParser(usage = "%prog [options] [vcf file]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need input vcf file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$geneBed)) {
    cat("Need input gene bed\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output file name\n\n")
    print_help(parser);
    stop();
}

fn <- arguments$args[1];
outfn <- opt$outFile
null <- suppressWarnings(file.remove(outfn))
out <- file(outfn, open = 'a')

cat('Reading target bed ... ')
genes <- import(opt$geneBed)
cat('done\n')

cat('Indexing vcf ... ')
temp <- tempfile()
zipped <- bgzip(fn, temp)
idx <- indexTabix(temp, "vcf")
cat('done\n')

tab <- TabixFile(zipped, idx, yieldSize = 2000)
open(tab)

cat('Processing vcf by chunk\n')
i <- 1
while(nrow(vcf <- readVcf(tab, genome = opt$genome))) {
    # replace header
    newInfo <- DataFrame(Number = 0, Type = "Flag", Description = paste(opt$name, ": variant is in gene list", sep = ''), row.names = opt$name)
    ol <- findOverlaps(rowData(vcf), genes, select = 'first')
    info(vcf)[,opt$name] <- !is.na(ol)

    #fix sample genotype order
    if ("GT" %in% names(geno(vcf))) {
        x <- which(names(geno(vcf)) == "GT")
        ord <- c(x, (1:length(geno(vcf)))[-x])
        geno(vcf) <- geno(vcf)[ord]
    }

    cat("Appending vcf chunk to", opt$outFile, "... ")
    writeVcf(vcf, out)
    cat("done\n")
}
close(tab)
close(out)
