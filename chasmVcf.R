#!/usr/bin/env Rscript
# Read a variant table and extract uniprot accession ids

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("biomaRt"));
suppressPackageStartupMessages(library("VariantAnnotation"));
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene));
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(org.Hs.eg.db))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--genome", default = 'hg19', help = "genome build [default %default]"),
                make_option("--chasmDir", default = NULL, help = "CHASM dir"),
                make_option("--classifier", default = 'Breast', help = "CHASM classifier [default %default]"),
                make_option("--outFile", default = stdout(), help = "vcf output file [default %default]"))

parser <- OptionParser(usage = "%prog vcf.file");
arguments <- parse_args(parser, positional_arguments = T, option_list = optList);
opt <- arguments$options;

if (is.null(opt$chasmDir)) {
    cat("Need CHASM dir\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) != 1) {
    cat("Reading from stdin...\n");
    vcf <- readVcf(stdin(), genome = opt$genome)
} else {
    fn <- arguments$args[1];
    vcf <- readVcf(fn, genome = opt$genome)
}

X <- cbind(seq = as.character(seqnames(rowData(vcf))), start = start(rowData(vcf)), end = start(rowData(vcf)) + width(rowData(vcf)), strand = '+', ref = as.character(rowData(vcf)$REF), alt = as.character(unlist(rowData(vcf)$ALT)))

setwd(opt$chasmDir)
tmp <- tempfile()
write.table(X, file = tmp, quote = F, sep = '\t', row.names = F, col.names = F)
cmd <- paste("./RunChasm ", opt$classifier, ' -g' ,sep = '')
system(cmd)
results <- read.table(file = paste(tmp, '.output', sep = ''), sep = '\t', header = T, as.is = T)

hinfoprime <- apply(as.data.frame(info(header(vcf))), 2, as.character)
rownames(hinfoprime) <- rownames(info(header(vcf)))
hinfoprime <- rbind(hinfoprime, chasm_mut = c("A", "String", "CHASM mutation"))
hinfoprime <- rbind(hinfoprime, chasm_pval = c("A", "Float", "CHASM p-value"))
hinfoprime <- rbind(hinfoprime, chasm_score = c("A", "Float", "CHASM score"))
hinfoprime <- rbind(hinfoprime, chasm_fdr = c("A", "Float", "CHASM B-H FDR"))
hinfoprime <- DataFrame(hinfoprime)
exptData(vcf)$header <- VCFHeader(samples = header(vcf)@samples, meta = meta(header(vcf)), info = hinfoprime, geno = geno(header(vcf)))

infoprime <- info(vcf)
infoprime[as.integer(results$MutationID),"chasm_mut"] <- as.character(results$Mutation)
infoprime[as.integer(results$MutationID),"chasm_score"] <- results$CHASM
infoprime[as.integer(results$MutationID),"chasm_pval"] <- results$PValue
infoprime[as.integer(results$MutationID),"chasm_fdr"] <- results$BHFDR
info(vcf) <- infoprime

writeVcf(vcf, opt$outFile)

