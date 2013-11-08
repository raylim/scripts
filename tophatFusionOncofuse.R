suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene));
suppressPackageStartupMessages(library(GenomicFeatures));
suppressPackageStartupMessages(library(org.Hs.eg.db));


txdb <- loadDb('~/share/reference/ensembl_biomart.sqlite')
txByGenes <- transcriptsBy(txdb, "gene")

fusions <- read.table('fusions.out', sep = '\t', quote = '', comment.char = '')[, 1:10]
colnames(fusions) <- c("Chromosomes", "Pos1", "Pos2", "Orientation", "NumSpan", "NumSplit", "NumSplitSpan", "NumContradict", "Num5pBases", "Num3pBases")

fusions <- transform(fusions, Strand1 = ifelse(grepl('f.', Orientation), "+", "-"))
fusions <- transform(fusions, Strand2 = ifelse(grepl('.f', Orientation), "+", "-"))
fusions <- transform(fusions, Chr1 = sub('-.*', '', Chromosomes))
fusions <- transform(fusions, Chr2 = sub('.*-', '', Chromosomes))

gr1 <- GRanges(seqnames = fusions$Chr1, ranges = IRanges(start = fusions$Pos1, width = 1), strand = fusions$Strand1)
gr2 <- GRanges(seqnames = fusions$Chr2, ranges = IRanges(start = fusions$Pos2, width = 1), strand = fusions$Strand2)

tx <- transcriptsByOverlaps(txdb, gr1)

