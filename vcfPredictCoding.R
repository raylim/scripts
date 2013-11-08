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
                make_option("--outFile", default = NULL, help = "vcf output file [default %default]"),
                make_option("--ensemblTxdb", default = NULL, help = "Ensembl TxDb SQLite"),
                make_option("--ref", default = NULL, help = "Reference fasta file")
                )
parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) != 1) {
    cat("Need vcf file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$ref)) {
    cat("Need reference fasta file\n");
    print_help(parser);
    stop();
}


#fn <- 'vcf/AdCC10T_AdCC10N.mutect.dp_ft.dbsnp.nsfp.chasm.vcf'
#opt$ref <- '/home/limr/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta'
#opt$ensemblTxdb <- '~/share/reference/ensembl_biomart.sqlite'
#opt$fathmmDir <- '~/share/usr/fathmm/'
#opt$genome <- 'hg19'
fn <- arguments$args[1];

cat('Reading vcf ... ')
vcf <- readVcf(fn, genome = opt$genome)
passIds <- which(rowData(vcf)$FILTER == "PASS")
if (length(passIds) == 0) {
    predCod <- NULL
    save(predCod, file = opt$outFile)
} else {
    vcfPass <- vcf[passIds, ]
    cat('done\n')

    if (is.null(opt$ensemblTxdb)) {
        txdb <- makeTranscriptDbFromBiomart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
    } else {
        txdb <- loadDb(opt$ensemblTxdb)
    }

    #saveDb(txdb, file = 'ensembl_biomart.sqlite')

    #ref = FaFile('/home/limr/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta')
    cat("Predicting coding from reference ...")
    ref = FaFile(opt$ref)
    predCod <- predictCoding(vcfPass, txdb, ref)
    cat(" done\n")

    save(predCod, file = opt$outFile)
}
