#suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene));
suppressPackageStartupMessages(library(GenomicFeatures));
suppressPackageStartupMessages(library(org.Hs.eg.db));


setwd('~/projects/salivary_rnaseq')
txdb <- loadDb('ensembl_biomart.sqlite')
txByGenes <- transcriptsBy(txdb, "gene")

fusions <- read.table('fusions.out', sep = '\t', quote = '', comment.char = '')[, 1:10]
colnames(fusions) <- c("Chromosomes", "Pos1", "Pos2", "Orientation", "NumSpan", "NumSplit", "NumSplitSpan", "NumContradict", "Num5pBases", "Num3pBases")

fusions <- transform(fusions, Strand1 = ifelse(grepl('f.', Orientation), "+", "-"))
fusions <- transform(fusions, Strand2 = ifelse(grepl('.f', Orientation), "+", "-"))
fusions <- transform(fusions, Chr1 = sub('-.*', '', Chromosomes))
fusions <- transform(fusions, Chr2 = sub('.*-', '', Chromosomes))

gr1 <- GRanges(seqnames = fusions$Chr1, ranges = IRanges(start = fusions$Pos1, width = 1), strand = fusions$Strand1)
gr2 <- GRanges(seqnames = fusions$Chr2, ranges = IRanges(start = fusions$Pos2, width = 1), strand = fusions$Strand2)

tx1 <- transcriptsByOverlaps(txdb, gr1)
olaps <- findOverlaps(tx1, gr1)
X <- tx1[queryHits(olaps)]
X$id <- subjectHits(olaps)
X$pos <- start(gr1)[subjectHits(olaps)]
tx1 <- X
#tx1$id[ <- subjectHits(findOverlaps(tx1, gr1))
egids <- mget(tx1$tx_name, org.Hs.egENSEMBLTRANS2EG, ifnotfound = NA)
egids <- sapply(egids, function (x) x[1])
tx1$egid <- egids


tx1 <- subset(tx1, !is.na(egid))

chrloc <- mget(tx1$egid, org.Hs.egCHRLOC)
chrloc <- sapply(chrloc, function (x) x[1])
strand <- ifelse(chrloc > 0, '+', '-')

tx1$tx_strand <- strand
tx1 <- subset(tx1, !is.na(strand))

tx1$genepos <- ifelse(tx1$tx_strand == as.character(strand(tx1)), "5prime", "3prime")

tx2 <- transcriptsByOverlaps(txdb, gr2)
olaps <- findOverlaps(tx2, gr2)
X <- tx2[queryHits(olaps)]
X$id <- subjectHits(olaps)
X$pos <- start(gr2)[subjectHits(olaps)]
tx2 <- X
egids <- mget(tx2$tx_name, org.Hs.egENSEMBLTRANS2EG, ifnotfound = NA)
egids <- sapply(egids, function (x) x[1])
tx2$egid <- egids

tx2 <- subset(tx2, !is.na(egid))

chrloc <- mget(tx2$egid, org.Hs.egCHRLOC)
chrloc <- sapply(chrloc, function (x) x[1])
strand <- ifelse(chrloc > 0, '+', '-')

tx2$tx_strand <- strand
tx2 <- subset(tx2, !is.na(strand))

tx2$genepos <- ifelse(tx2$tx_strand == as.character(strand(tx2)), "5prime", "3prime")

XX <- as.data.frame(tx1)
YY <- as.data.frame(tx2)
Data <- merge(XX, YY, by = "id", suffixes = c(".1", ".2"))
Data <- subset(Data, genepos.1 != genepos.2)
Data <- transform(Data, upstreamChrom = ifelse(genepos.1 == "5prime", seqnames.1, seqnames.2))
Data <- transform(Data, upstreamPos = ifelse(genepos.1 == "5prime", pos.1, pos.2))
Data <- transform(Data, downstreamChrom = ifelse(genepos.1 != "5prime", seqnames.1, seqnames.2))
Data <- transform(Data, downstreamPos = ifelse(genepos.1 != "5prime", pos.1, pos.2))

oncofuseCoords <- Data[, c("upstreamChrom", "upstreamPos", "downstreamChrom", "downstreamPos")]
if (any(duplicated(oncofuseCoords))) {
    oncofuseCoords <- oncofuseCoords[!duplicated(oncofuseCoords), ]
}
