#!/usr/bin/env Rscript
# plots control freec 

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("VariantAnnotation"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("org.Hs.eg.db"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--out", default = NULL, help = "Output file (pdf)"),
                make_option("--ref", default = 'hg19', help = "ref genome (default = %default)"),
                make_option("--geneList", help = "gene list"));

parser <- OptionParser(usage = "%prog [options] [list of annotated vcf files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need input vcf files\n");
    print_help(parser);
    stop();
} else {
    files <- arguments$args;
}

effFormat <- c("impact", "class", "codon change", "aa change", "aa length", "gene", "transcript biotype", "gene coding", "transcript id", "exon rank", "genotype number")

for (vcfFile in files) {
    vcf <- readVcf(vcfFile, opt$ref)
    effLists <- info(vcf)$EFF
    effDfs <- list()
    for (i in 1:length(effLists)) {
        effList <- effLists[[i]]
        if (length(effList) > 0) {
            effF <- sub('.*\\(', '', sub('\\)', '', effList))
            effM <- do.call('rbind', lapply(strsplit(effF, '\\|'), function(x) x[1:length(effFormat)]))
            effM[effM == ""] <- NA
            colnames(effM) <- effFormat
            effType <- tolower(gsub('_', ' ', sub('\\(.*', '', effList)))
            effDfs[[i]]<- data.frame(chr = seqnames(rowData(vcf))[i], pos = start(rowData(vcf))[i], id = names(rowData(vcf))[i], eff.type = effType, effM)
        }
    }
    effDf <- do.call('rbind', effDfs)

}



