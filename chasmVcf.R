#!/usr/bin/env Rscript
# Read a variant table and extract uniprot accession ids

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("VariantAnnotation"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--genome", default = 'hg19', help = "genome build [default %default]"),
                make_option("--chasmDir", default = NULL, help = "CHASM dir"),
                make_option("--pythonDir", default = NULL, help = "python directory"),
                make_option("--classifier", default = 'Breast', help = "CHASM classifier [default %default]"),
                make_option("--outFile", default = stdout(), help = "vcf output file [default %default]"))

parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$chasmDir)) {
    cat("Need CHASM dir\n");
    print_help(parser);
    stop();
} else if (is.null(opt$pythonDir)) {
    cat("Need python dir\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    #cat("Reading from stdin...\n");
    tmp <- tempfile()
    f <- file('stdin')
    open(f)
    x <- readLines(f)
    close(f)
    write(x, file = tmp)
    vcf <- readVcf(tmp, genome = opt$genome)
} else {
    fn <- arguments$args[1];
    vcf <- readVcf(fn, genome = opt$genome)
}

al <- sapply(rowData(vcf)$ALT, length)
X <- cbind(id = 1:nrow(vcf), seq = as.character(seqnames(rowData(vcf))), zero = start(rowData(vcf)) - 1, one = start(rowData(vcf)) + 1, strand = rep('+', nrow(vcf)), ref = as.character(rowData(vcf)$REF), filt = rowData(vcf)$FILTER)
XX <- apply(X, 2, rep, times = al)
X <- cbind(XX, alt = as.character(unlist(rowData(vcf)$ALT)))
X <- as.data.frame(X, stringsAsFactors = F)
X <- subset(X, sapply(alt, nchar) == 1 & sapply(ref, nchar) == 1 & filt == "PASS")
X <- X[,-which(colnames(X) == 'filt')]
X$seq <- sub('^', 'chr', X$seq)

oldwd <- getwd()
setwd(opt$chasmDir)
tmp <- tempfile()
write.table(X, file = tmp, quote = F, sep = '\t', row.names = F, col.names = F)
cmd <- paste("CHASMDIR=", opt$chasmDir, " PATH=", opt$pythonDir, "/bin ./RunChasm ", opt$classifier, ' ', tmp, ' -g' ,sep = '')
#cat(cmd)
system(cmd, ignore.stdout = T)
results <- read.table(file = paste(tmp, opt$classifier, '.output', sep = ''), sep = '\t', header = T, as.is = T)

hinfoprime <- apply(as.data.frame(info(header(vcf))), 2, as.character)
rownames(hinfoprime) <- rownames(info(header(vcf)))
hinfoprime <- rbind(hinfoprime, chasm_mut = c("A", "String", "CHASM mutation"))
hinfoprime <- rbind(hinfoprime, chasm_pval = c("A", "Float", "CHASM p-value"))
hinfoprime <- rbind(hinfoprime, chasm_score = c("A", "Float", "CHASM score"))
hinfoprime <- rbind(hinfoprime, chasm_fdr = c("A", "Float", "CHASM B-H FDR"))
hinfoprime <- DataFrame(hinfoprime, row.names = rownames(hinfoprime))
hlist <- header(exptData(vcf)$header)
hlist$INFO <- hinfoprime
exptData(vcf)$header <- new("VCFHeader", samples = header(vcf)@samples, header = hlist)

if (nrow(results) > 1) {
    infoprime <- info(vcf)
    infoprime[as.integer(as.character(results$MutationID)),"chasm_mut"] <- as.character(results$Mutation)
    infoprime[as.integer(as.character(results$MutationID)),"chasm_score"] <- results$CHASM
    infoprime[as.integer(as.character(results$MutationID)),"chasm_pval"] <- results$PValue
    if (nrow(results) > 10) {
        infoprime[as.integer(as.character(results$MutationID)),"chasm_fdr"] <- results$BHFDR
    }
    info(vcf) <- infoprime
}

setwd(oldwd)
writeVcf(vcf, opt$outFile)
