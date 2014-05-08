#!/usr/bin/env Rscript
# rbinds together tab-delimited tables and outputs to STDOUT

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("org.Hs.eg.db"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--outDir", default = NULL, help = "Output dir"),
                make_option("--includeChrY", action = 'store_true', default = F, help = "include Y chromosome"),
                make_option("--knownVariants", default = NULL, help = "known variants file"),
                make_option("--txdb", default = NULL, help = "txdb"))

parser <- OptionParser(usage = "%prog [options] [list of ratio.txt files]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need input controlFreeC CNV files\n");
    print_help(parser);
    stop();
}

files <- arguments$args;


grs <- list()
tables <- list()
samples <- c()
for (f in files) {
    s <- sub('\\..*', '', f)
    samples <- c(samples, s)
    d <- read.table(file = f, sep = '\t', header = F, as.is = T, comment.char = '');
    tables[[s]] <- d
    grs[[s]] <- GRanges(seqnames = d[, 1], ranges = IRanges(d[, 2], end = d[, 3]), copynum = d[,4])
}

posns <- do.call('rbind', lapply(tables, function(x) x[,1:3]))
colnames(posns) <- c("chr", "start", "end")
rownames(posns) <- NULL

gr <- GRanges(seqnames = posns$chr, ranges = IRanges(posns$start, end = posns$end))
gr <- disjoin(gr)
gr <- gr[width(gr) > 1]
if (!opt$includeChrY) {
    gr <- gr[seqnames(gr) != "Y"]
}
x <- as.vector(seqnames(gr))
if (any(x == "X")) {
    x[x == "X"] <- 23
}
if (any(x == "Y")) {
    x[x == "Y"] <- 24
}
if (any(x == "MT")) {
    x[x == "MT"] <- 25
}
x <- as.integer(x)
oo <- order(x, start(gr))
gr <- gr[oo, ]
for (s in samples) {
    mcols(gr)[, s] <- rep(2, length(gr))
}

for (s in samples) {
    overlaps <- findOverlaps(grs[[s]], gr)
    mcols(gr[subjectHits(overlaps), ])[[s]] <- mcols(grs[[s]][queryHits(overlaps), ])$copynum
}

X <- as.matrix(mcols(gr))
X[X > 3] <- 3


if (!is.null(opt$txdb)) {
    txdb <- loadDb(opt$txdb)
} else {
    txdb <- makeTranscriptDbFromBiomart('ensembl', 'hsapiens_gene_ensembl')
}

txs <- transcriptsBy(txdb, by = "gene")
overlaps <- findOverlaps(gr, txs)
ensids <- names(txs[subjectHits(overlaps)])
x <- tapply(ensids, queryHits(overlaps), function(x) {
    egids <- sapply(mget(x, org.Hs.egENSEMBL2EG, ifnotfound = NA), function(x) x[1]);
    if (sum(!is.na(egids) > 0)) {
        genes <- sapply(mget(egids[!is.na(egids)], org.Hs.egSYMBOL), function(x) x[1]);
        paste(genes, collapse = "|");
    } else {
        ""
    }
})
mcols(gr)$Genes <- rep("", length(gr))
mcols(gr)[unique(queryHits(overlaps)), 'Genes'] <- x

genes <- unique(unlist(sapply(mcols(gr)$Genes, strsplit, '\\|')))
geneCN <- matrix(2, nrow = length(genes), ncol = length(samples), dimnames = list(Gene = genes, Sample = samples))
for (i in 1:length(gr)) {
    gs <- unlist(strsplit(mcols(gr)[i, "Genes"], '\\|'))
    x <- as.integer(as.data.frame(mcols(gr))[i, samples])
    if (sum(x != 2) > 1) {
        for (g in gs) {
            geneCN[g, samples][x != 2] <- x[x != 2]
        }
    }
}

fn <- paste(opt$outDir, "/gene_copynum_recurrent.txt", sep = "")
x <- rowSums(geneCN != 2) > 1
geneRecurrentCNV <- geneCN[x, ]
write.table(geneRecurrentCNV, file = fn, sep = "\t", quote = F)

fn <- paste(opt$outDir, "/gene_copynum_recurrent_gain.txt", sep = "")
x <- rowSums(geneCN > 2) > 1
geneRecurrentGainCNV <- geneCN[x, ]
write.table(geneRecurrentGainCNV, file = fn, sep = "\t", quote = F)

fn <- paste(opt$outDir, "/gene_copynum_recurrent_loss.txt", sep = "")
x <- rowSums(geneCN < 2) > 1
geneRecurrentLossCNV <- geneCN[x, ]
write.table(geneRecurrentLossCNV, file = fn, sep = "\t", quote = F)

if (!is.null(opt$knownVariants)) {
    dgv <- read.table(opt$knownVariants, header = T, quote = '', comment.char = '', as.is = T, sep = '\t')
        grs[[s]] <- GRanges(seqnames = d[, 1], ranges = IRanges(d[, 2], end = d[, 3]), copynum = d[,4])
    dgv.gr <- GRanges(seqnames = dgv$chr, ranges = IRanges(dgv$start, end = dgv$end), id = dgv$variantaccession, observedgains = dgv$observedgains, observedlosses = dgv$observedlosses)
    small.dgv.gr <- sort(dgv.gr[width(dgv.gr) <= 500000])
    ols <- findOverlaps(gr, small.dgv.gr, minoverlap = 25000L)
    rs <- ranges(ols, ranges(gr), ranges(small.dgv.gr))
    olf <- width(rs) / width(gr)[queryHits(ols)]
    x <- ols[olf > 0.5]
    mcols(gr)$ID <- rep(".", length(gr))
    ids <- mcols(dgv.gr[subjectHits(x)])$id
    ids <- tapply(ids, queryHits(x), paste, collapse = '|')
    small.gr.index <- which(width(gr) <= 500000)
    mcols(gr)$ID[intersect(as.integer(names(ids)), small.gr.index)] <- ids
    mcols(gr) <- cbind(mcols(gr), mcols(small.dgv.gr[subjectHits(x)]))
}

fn <- paste(opt$outDir, "/annotated_cnv.txt", sep = "")
write.table(as.data.frame(gr), file = fn, sep = "\t", row.names = F, quote = F)
