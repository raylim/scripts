#!/usr/bin/env Rscript
# Read a vcf file and append fathmm results

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("biomaRt"));
suppressPackageStartupMessages(library("VariantAnnotation"));
suppressPackageStartupMessages(library("data.table"));
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene));
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(org.Hs.eg.db))

options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))

optList <- list(
        make_option("--genome", default = 'hg19', help = "genome build [default %default]"),
        make_option("--fathmmDir", default = NULL, help = "fathmm dir"),
        make_option("--fathmmAlg", default = 'Cancer', help = "fathmm algorithm [default %default]"),
        make_option("--fathmmOnt", default = 'DO', help = "fathmm ontology [default %default]"),
        make_option("--ensemblTxdb", default = NULL, help = "Ensembl TxDb SQLite"),
        make_option("--ref", default = NULL, help = "Reference fasta file"),
        make_option("--python", default = 'python', help = "python executable [default %default]"),
        make_option("--outFile", default = NULL, help = "vcf output file [default %default]")
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
} else if (is.null(opt$ref)) {
    cat("Need reference fasta file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) != 1) {
    cat("Need vcf file\n");
    print_help(parser);
    stop();
}

fn <- arguments$args[1];
outfn <- opt$outFile
out <- file(outfn, open = 'a')

cat('Loading transcriptdb ... ')
if (is.null(opt$ensemblTxdb)) {
    txdb <- makeTranscriptDbFromBiomart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
} else {
    txdb <- loadDb(opt$ensemblTxdb)
}
cat('done\n')

ref = FaFile(opt$ref)

cat('Connecting to ensembl ... ')
ensembl = useMart("ensembl") #, host = 'localhost', port = 9000)
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
cat('done\n')

#fn <- 'vcf/AdCC10T_AdCC10N.mutect.dp_ft.dbsnp.nsfp.chasm.vcf'
#opt$ref <- '/home/limr/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta'
#opt$fathmmDir <- '~/share/usr/fathmm/'
#opt$genome <- 'hg19'
fn <- arguments$args[1];

cat('Reading vcf header ... ')
# create new header
vcfHeader <- scanVcfHeader(fn)
hinfoprime <- apply(as.data.frame(info(vcfHeader)), 2, as.character)
rownames(hinfoprime) <- rownames(info(vcfHeader))
hinfoprime <- rbind(hinfoprime, fathmm_query = c("A", "String", "fathmm query"))
hinfoprime <- rbind(hinfoprime, fathmm_pred = c("A", "String", "fathmm prediction"))
hinfoprime <- rbind(hinfoprime, fathmm_score = c("A", "Float", "fathmm score"))
hinfoprime <- DataFrame(hinfoprime, row.names = rownames(hinfoprime))
hlist <- header(vcfHeader)
hlist$INFO <- hinfoprime
newVcfHeader <- new("VCFHeader", samples = vcfHeader@samples, header = hlist)
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

    exptData(vcf)$header <- newVcfHeader

    cat(paste('Chunk', i, "\n"))
    passIds <- which(rowData(vcf)$FILTER == "PASS")
    if (length(passIds) == 0) {
        cat("No unfiltered variants\n")
        next()
    }

    cat("Predicting coding from reference...\n")
    predCod <- predictCoding(vcf[passIds, ], txdb, ref)
    cat(" done\n")

    if (sum(predCod$CONSEQUENCE == "nonsynonymous") == 0) {
        cat("No non-syn variants\n")
        next()
    }

    predCod <- subset(predCod, CONSEQUENCE == "nonsynonymous")

    # retrieve transcript ids
    x <- transcripts(txdb, vals = list(tx_id = predCod$TXID), columns = c('tx_id', 'tx_name'))
    enstIds <- x$tx_name
    names(enstIds) <- x$tx_id
    aa = cbind(queryId = predCod$QUERYID, aa = paste(as.character(predCod$REFAA), lapply(predCod$PROTEINLOC, function(x) x[1]), as.character(predCod$VARAA), sep = ''))
    rownames(aa) <- predCod$TXID

    cat("Looking up ensembl peptide IDs ... ")
    ids <- getBM(filters = 'ensembl_transcript_id', attributes = c('ensembl_transcript_id', 'ensembl_peptide_id'), values = enstIds, mart = ensembl)
    rownames(ids) <- names(enstIds)[match(ids$ensembl_transcript_id, enstIds)]
    cat("done\n")

    ids <- cbind(aa, ids[rownames(aa), ])
    fathmmInput <- subset(ids, ensembl_peptide_id != "", select = c('ensembl_peptide_id', 'aa'))

    cat("Calling fathmm: ")
    oldwd <- getwd()
    tmp1 <- tempfile()
    tmp2 <- tempfile()
    setwd(paste(opt$fathmmDir, '/cgi-bin', sep = ''))
    cmd <- paste(opt$python, 'fathmm.py -w', opt$fathmmAlg, '-p', opt$fathmmOnt, tmp1, tmp2)
    write.table(subset(ids, ensembl_peptide_id != "", select = c('ensembl_peptide_id', 'aa')), file = tmp1, quote = F, sep = ' ', row.names = F, col.names = F)
    #cmd <- paste('python fathmm.py -w Cancer', tmp1, tmp2)
    system(cmd)
    cat("\ndone\n")

    cat("Reading results ... ")
    results <- read.table(tmp2, sep = '\t', header = T, comment.char = '', row.names = 1, quote = '')
    cat("done\n")
    results <- merge(ids, results, by.x = c('aa', 'ensembl_peptide_id'), by.y = c('Substitution', 'Protein.ID'))

    split.results <- split(results, factor(results$queryId))
    cat("Selecting minimum scores ... ")
    results <- rbindlist(lapply(split.results, function(x) x[which.min(x$Score), ]))
    cat("done\n")

    if (!is.null(results) && nrow(results) > 0) {
        cat("Merging fathmm results ... ")
        infodprime <- info(vcf[passIds, ])
        infodprime[as.integer(as.character(results$queryId)),"fathmm_query"] <- with(results, paste(ensembl_peptide_id, aa, sep = "_"))
        infodprime[as.integer(as.character(results$queryId)),"fathmm_pred"] <- as.character(results$Prediction)
        infodprime[as.integer(as.character(results$queryId)),"fathmm_score"] <- results$Score
        info(vcf)[passIds, ] <- infodprime
        cat("done\n")
    } else {
        cat("No results from fathmm\n")
    }

    #fix sample genotype order
    x <- which(names(geno(vcf)) == "GT")
    ord <- c(x, (1:length(geno(vcf)))[-x])
    geno(vcf) <- geno(vcf)[ord]

    cat("Appending vcf chunk to", opt$outFile, "... ")
    setwd(oldwd)
    writeVcf(vcf, out)
    cat("done\n")

    i <- i + 1
}

close(tab)
close(out)


