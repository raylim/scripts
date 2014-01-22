#!/usr/bin/env Rscript
# Read a variant table and extract uniprot accession ids

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("biomaRt"));
suppressPackageStartupMessages(library("VariantAnnotation"));
suppressPackageStartupMessages(library("AnnotationDbi"));

options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))

optList <- list(
                make_option("--genome", default = 'hg19', help = "genome build [default %default]"),
                make_option("--ensemblTxdb", default = NULL, help = "Ensembl TxDb SQLite"),
                make_option("--predCoding", default = NULL, help = "predicted coding R data file"),
                make_option("--outFile", default = NULL, help = "Rdata output file [default %default]")
                )
parser <- OptionParser(usage = "%prog predCoding.Rdata", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) != 1) {
    cat("Need predicted coding Rdata file\n");
    print_help(parser);
    stop();
}


#fn <- 'vcf/AdCC10T_AdCC10N.mutect.dp_ft.dbsnp.nsfp.chasm.vcf'
#opt$ref <- '/home/limr/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta'
#opt$ensemblTxdb <- '~/share/reference/ensembl_biomart.sqlite'
#opt$fathmmDir <- '~/share/usr/fathmm/'
#opt$genome <- 'hg19'
fn <- arguments$args[1];

cat("Loading predicted coding R data file ... ")
load(fn)

if (sum(predCod$CONSEQUENCE == "nonsynonymous") == 0) {
    fathmmInput <- NULL
    ids <- NULL
    save(fathmmInput, ids, file = opt$outFile)
} else {
    predCod <- subset(predCod, CONSEQUENCE == "nonsynonymous")
    cat(" done\n")

    if (is.null(opt$ensemblTxdb)) {
        txdb <- makeTranscriptDbFromBiomart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
    } else {
        txdb <- loadDb(opt$ensemblTxdb)
    }

    cat('Connecting to ensembl ... ')
    ensembl = useMart("ensembl") #, host = 'localhost', port = 9000)
    ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
    cat('done\n')

    #saveDb(txdb, file = 'ensembl_biomart.sqlite')

    #txByGene <- transcriptsBy(txdb, 'gene')


    x <- transcripts(txdb, vals = list(tx_id = predCod$TXID), columns = c('tx_id', 'tx_name'))
    enstIds <- x$tx_name
    names(enstIds) <- x$tx_id
    aa = cbind(queryId = predCod$QUERYID, aa = paste(as.character(predCod$REFAA), lapply(predCod$PROTEINLOC, function(x) x[1]), as.character(predCod$VARAA), sep = ''))
    rownames(aa) <- predCod$TXID

    cat("Looking up ensembl peptide IDs ... ")
    ids <- getBM(filters = 'ensembl_transcript_id', attributes = c('ensembl_transcript_id', 'ensembl_peptide_id'), values = enstIds, mart = ensembl)
    rownames(ids) <- names(enstIds)[match(ids$ensembl_transcript_id, enstIds)]
    cat("done\n")

    #ids <- cbind(enst = enstIds, esnp = enspIds[names(enstIds), "ensembl_peptide_id"])
    ids <- cbind(aa, ids[rownames(aa), ])

    #X <- cbind(aa, enspIds[names(aa), ])

    fathmmInput <- subset(ids, ensembl_peptide_id != "", select = c('ensembl_peptide_id', 'aa'))

    save(fathmmInput, ids, file = opt$outFile)
}

