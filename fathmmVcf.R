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
        make_option("--python", default = 'python', help = "python executable [default %default]"),
        make_option("--outFile", default = NULL, help = "vcf output file [default %default]")
        )
parser <- OptionParser(usage = "%prog vcf.file fathmm.input.Rdata", option_list = optList);
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
} else if (length(arguments$args) != 2) {
    cat("Need vcf file and fathmm input Rdata\n");
    print_help(parser);
    stop();
}


#fn <- 'vcf/AdCC10T_AdCC10N.mutect.dp_ft.dbsnp.nsfp.chasm.vcf'
#opt$ref <- '/home/limr/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta'
#opt$fathmmDir <- '~/share/usr/fathmm/'
#opt$genome <- 'hg19'
fn <- arguments$args[1];
fn2 <- arguments$args[2];

cat('Reading vcf ... ')
vcf <- readVcf(fn, genome = opt$genome)
passIds <- which(rowData(vcf)$FILTER == "PASS")
cat('done\n')

load(fn2)

oldwd <- getwd()
if (!is.null(ids)) {
    cat("Calling fathmm: ")

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

    hinfoprime <- apply(as.data.frame(info(header(vcf))), 2, as.character)
    rownames(hinfoprime) <- rownames(info(header(vcf)))
    hinfoprime <- rbind(hinfoprime, fathmm_query = c("A", "String", "fathmm query"))
    hinfoprime <- rbind(hinfoprime, fathmm_pred = c("A", "String", "fathmm prediction"))
    hinfoprime <- rbind(hinfoprime, fathmm_score = c("A", "Float", "fathmm score"))
    hinfoprime <- DataFrame(hinfoprime, row.names = rownames(hinfoprime))
    hlist <- header(exptData(vcf)$header)
    hlist$INFO <- hinfoprime
    exptData(vcf)$header <- new("VCFHeader", samples = header(vcf)@samples, header = hlist)

    infoprime <- info(vcf)
    infoprime[, "fathmm_query"] <- rep(NA, nrow(vcf))
    infoprime[, "fathmm_pred"] <- rep(NA, nrow(vcf))
    infoprime[, "fathmm_score"] <- rep(NA, nrow(vcf))
    info(vcf) <- infoprime

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
}
#fix sample genotype order
x <- which(names(geno(vcf)) == "GT")
ord <- c(x, (1:length(geno(vcf)))[-x])
geno(vcf) <- geno(vcf)[ord]

cat("Writing vcf to", opt$outFile, "... ")
setwd(oldwd)
writeVcf(vcf, opt$outFile)
cat("done\n")



