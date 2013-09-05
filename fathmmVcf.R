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
                make_option("--fathmmDir", default = NULL, help = "fathmm dir"),
                make_option("--fathmmAlg", default = 'Cancer', help = "fathmm algorithm [default %default]"),
                make_option("--fathmmOnt", default = 'DO', help = "fathmm ontology [default %default]"),
                make_option("--ensemblTxdb", default = NULL, help = "Ensembl TxDb SQLite"),
                make_option("--ref", default = NULL, help = "Reference fasta file"),
                make_option("--outFile", default = NULL, help = "vcf output file [default %default]"),
                make_option("--python", default = 'python', help = "python executable [default %default]")
                )
parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$fathmmDir)) {
    cat("Need fathmm dir\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
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


fn <- arguments$args[1];

cat('Reading vcf ... ')
vcf <- readVcf(fn, genome = opt$genome)
cat('done\n')

if (is.null(opt$ensemblTxdb)) {
    txdb <- makeTranscriptDbFromBiomart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
} else {
    txdb <- loadDb(opt$ensemblTxdb)
}

cat('Connecting to ensembl ... ')
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
cat('done\n')

#saveDb(txdb, file = 'ensembl_biomart.sqlite')

#txByGene <- transcriptsBy(txdb, 'gene')

#vcf <- vcf[rowData(vcf)$FILTER == "PASS", ]
#ref = FaFile('/home/limr/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta')
cat("Predicting coding from reference ...")
ref = FaFile(opt$ref)
predCod <- predictCoding(vcf, txdb, ref)
predCod <- subset(predCod, CONSEQUENCE == "nonsynonymous")
cat(" done\n")


x <- transcripts(txdb, vals = list(tx_id = predCod$TXID), columns = c('tx_id', 'tx_name'))
enstIds <- x$tx_name
names(enstIds) <- x$tx_id
aa = cbind(queryId = predCod$QUERYID, aa = paste(as.character(predCod$REFAA), predCod$PROTEINLOC, as.character(predCod$VARAA), sep = ''))
rownames(aa) <- predCod$TXID

cat("Looking up ensembl peptide IDs ... ")
ids <- getBM(filters = 'ensembl_transcript_id', attributes = c('ensembl_transcript_id', 'ensembl_peptide_id'), values = enstIds, mart = ensembl)
rownames(ids) <- names(enstIds)[match(ids$ensembl_transcript_id, enstIds)]
cat("done\n")

#ids <- cbind(enst = enstIds, esnp = enspIds[names(enstIds), "ensembl_peptide_id"])
ids <- cbind(aa, ids[rownames(aa), ])

#X <- cbind(aa, enspIds[names(aa), ])
cat("Calling fathmm: ")
oldwd <- getwd()
setwd(paste(opt$fathmmDir, '/cgi-bin', sep = ''))
tmp1 <- tempfile()
tmp2 <- tempfile()
write.table(subset(ids, ensembl_peptide_id != "", select = c('ensembl_peptide_id', 'aa')), file = tmp1, quote = F, sep = ' ', row.names = F, col.names = F)
cmd <- paste(opt$python, 'fathmm.py -w', opt$fathmmAlg, '-p', opt$fathmmOnt, tmp1, tmp2)
#cmd <- paste('python fathmm.py -w Cancer', tmp1, tmp2)
system(cmd)
results <- read.table(tmp2, sep = '\t', header = T, comment.char = '', row.names = 1)
results <- merge(ids, results, by.x = c('aa', 'ensembl_peptide_id'), by.y = c('Substitution', 'Protein.ID'))

hinfoprime <- apply(as.data.frame(info(header(vcf))), 2, as.character)
rownames(hinfoprime) <- rownames(info(header(vcf)))
hinfoprime <- rbind(hinfoprime, fathmm = c("A", "String", "fathmm prediction"))
hinfoprime <- rbind(hinfoprime, fathmm_score = c("A", "Float", "fathmm score"))
hinfoprime <- DataFrame(hinfoprime, row.names = rownames(hinfoprime))
hlist <- header(exptData(vcf)$header)
hlist$INFO <- hinfoprime
exptData(vcf)$header <- new("VCFHeader", samples = header(vcf)@samples, header = hlist)

if (nrow(results) > 0) {
    cat("Merging fathmm results ... ")
    infoprime <- info(vcf)
    infoprime[as.integer(results$queryId),"fathmm"] <- as.character(results$Prediction)
    infoprime[as.integer(results$queryId),"fathmm_score"] <- results$Score
    info(vcf) <- infoprime
    cat("done\n")
} else {
    cat("No results from fathmm\n")
}

cat("Writing vcf to", opt$outFile, "... ")
setwd(oldwd)
writeVcf(vcf, opt$outFile)
cat("done\n")



