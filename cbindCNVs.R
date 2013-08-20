#!/usr/bin/env Rscript
# rbinds together tab-delimited tables and outputs to STDOUT

suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("org.Hs.eg.db"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

arguments <- commandArgs(T);

arguments <- paste(samples, '.bam_CNVs', sep = '')

grs <- list()
tables <- list()
samples <- c()
for (arg in arguments) {
    s <- sub('\\..*', '', arg)
    samples <- c(samples, s)
    d <- read.table(file = arg, sep = '\t', header = F, as.is = T, comment.char = '');
    tables[[s]] <- d
    grs[[s]] <- GRanges(seqnames = d[, 1], ranges = IRanges(d[, 2], end = d[, 3]), copynum = d[,4])
}

posns <- do.call('rbind', lapply(tables, function(x) x[,1:3]))
colnames(posns) <- c("chr", "start", "end")
rownames(posns) <- NULL

gr <- GRanges(seqnames = posns$chr, ranges = IRanges(posns$start, end = posns$end))
gr <- disjoin(gr)
gr <- gr[width(gr) > 1]
x <- as.vector(seqnames(gr))
x[x == "X"] <- 23
x[x == "Y"] <- 24
x[x == "MT"] <- 25
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
#x <- round(width(gr) / 50000)
#Z <- apply(X, 2, rep, times = x)
#sn <- rep(as.vector(seqnames(gr)), times = x)
#chrstart <- seqnames(gr)

#txdb <- makeTranscriptDbFromUCSC('hg19', tablename = 'knownGene')
txdb <- makeTranscriptDbFromBiomart('ensembl', 'hsapiens_gene_ensembl')
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

d <- as.data.frame(gr)[,-9]
x <- rowSums(d[,c("AdCC10T", "AdCC1T", "AdCC2T", "AdCCPCT")] != 2) > 1
write.table(d[x,], file = "adcc_cnv_genes_sans9.txt", sep = "\t", row.names = F, quote = F)


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

geneCN <- geneCN[, -4]

x <- rowSums(geneCN != 2) > 1
geneRecurrentCNV <- geneCN[x, ]
write.table(geneRecurrentCNV, file = "gene_recurrent_cnv.txt", sep = "\t", quote = F)

x <- rowSums(geneCN > 2) > 1
geneRecurrentGainCNV <- geneCN[x, ]
write.table(geneRecurrentGainCNV, file = "gene_recurrent_gain_cnv.txt", sep = "\t", quote = F)

x <- rowSums(geneCN < 2) > 1
geneRecurrentLossCNV <- geneCN[x, ]
write.table(geneRecurrentLossCNV, file = "gene_recurrent_loss_cnv.txt", sep = "\t", quote = F)


arguments <- paste(samples, '.bam_ratio.txt', sep = '')

grs <- list()
tables <- list()
samples <- c()
for (arg in arguments) {
    s <- sub('\\..*', '', arg)
    samples <- c(samples, s)
    d <- read.table(file = arg, sep = '\t', header = T, as.is = T, comment.char = '');
    tables[[s]] <- d
    grs[[s]] <- GRanges(seqnames = d[, 1], ranges = IRanges(d[, 2], width = 50000), medianratio = d[,4], copynum = d[,5])
}

posns <- do.call('rbind', lapply(tables, function(x) x[,1:3]))
colnames(posns) <- c("chr", "start" )
rownames(posns) <- NULL

gr <- GRanges(seqnames = posns$chr, ranges = IRanges(posns$start, width = 50000))
gr <- disjoin(gr)
gr <- gr[width(gr) > 1]
gr <- gr[seqnames(gr) != "Y"]
x <- as.vector(seqnames(gr))
x[x == "X"] <- 23
#x[x == "MT"] <- 25
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

centromerePosns <- subset(read.table('../centromere_posns.txt', sep = '\t'), V8 == "centromere", as.is = T)
centromerePosns <- centromerePosns[, 2:4]
centromerePosns[, 1] <- sub('chr', '', centromerePosns[, 1])
centromerePosns[centromerePosns[,1] == "X", 1] <- 23
centromerePosns <- centromerePosns[-which(centromerePosns[,1] == "Y"), ]
oo <- order(as.integer(centromerePosns[,1]))
centromerePosns <- centromerePosns[oo, ]
centromerePosns[centromerePosns[,1] == 23, 1] <- "X"
centromereGR <- GRanges(seqnames = centromerePosns[, 1], ranges = IRanges(centromerePosns[, 2], end = centromerePosns[, 3]))

cmPos <- nearest(centromereGR, gr)
cmPos <- (cmPos - 1) / (length(gr) - 1)

chrPos <- (start(seqnames(gr)) - 1) / (length(gr) - 1)

X <- as.matrix(mcols(gr))
X[X > 3] <- 3
X <- X[,-4]
#x <- round(width(gr) / 50000)
#Z <- apply(X, 2, rep, times = x)
#sn <- rep(as.vector(seqnames(gr)), times = x)
#chrstart <- seqnames(gr)


cols <- c('black', 'blue', 'white', 'red')
png("adcc_copynum.png", height = 1000, width = 30000)
image(X, col = cols, axes = F)
axis(1, at = chrPos, labels = as.character(runValue(seqnames(gr))))
abline(v = chrPos, col = 'grey')
abline(v = cmPos, lty = 2, col = 'grey')
axis(2, at = 0:(ncol(X) - 1) / (ncol(X) - 1), labels = colnames(X))
box()
#legend('top', legend = c("deletion", "loss", "neutral", "gain"), fill = cols, horiz = T)
null <- dev.off()

