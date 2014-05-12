#!/usr/bin/env Rscript
# cluster vcf and output dendrogram

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("VariantAnnotation"));

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
} else if (length(arguments$args) < 1) {
    cat("Need vcf files\n");
    print_help(parser);
    stop();
}

vcfFile <- arguments$args[1]


vcf <- readVcf(vcfFile, opt$genome)
gt <- geno(vcf)$GT
X <- matrix(0, nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))
X[gt == "0/0"] <- 0
X[gt == "0/1"] <- 1
X[gt == "1/1"] <- 2
X[!gt %in% c("0/0", "0/1", "1/1")] <- NA
plot(hclust(dist(t(X), method = 'manhattan')))

gt <- matrix(as.integer(factor(gt)), nrow = nrow(gt), ncol = ncol(gt), dimnames = list(rownames(gt), colnames(gt)))

png(opt$outFile, height = 900, width = 1000)
null <- plot(hclust(dist(t(gt)), method = 'ward'))
dev.off()

