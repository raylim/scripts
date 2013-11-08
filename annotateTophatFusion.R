
#!/usr/bin/env Rscript
# Read a vcf file and append fathmm results

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene));
suppressPackageStartupMessages(library(org.Hs.eg.db))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

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

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txByGenes <- transcriptsBy(txdb, "gene")

fusions <- read.table('fusions.out', sep = '\t', quote = '', comment.char = '')[, 1:10]
colnames(fusions) <- c("Chromosomes", "Pos1", "Pos2", "Orientation", "NumSpan", "NumSplit", "NumSplitSpan", "NumContradict", "Num5pBases", "Num3pBases")
fusions <- transform(fusions, Chr1 = sub('-.*', '', Chromosomes))
fusions <- transform(fusions, Chr2 = sub('.*-', '', Chromosomes))
fusions$Chr1 <- sub('^', 'chr', fusions$Chr1)
fusions$Chr2 <- sub('^', 'chr', fusions$Chr2)
fusions <- transform(fusions, Strand1 = ifelse(grepl('f.', Orientation), "+", "-"))
fusions <- transform(fusions, Strand2 = ifelse(grepl('.f', Orientation), "+", "-"))

gr1 <- GRanges(seqnames = fusions$Chr1, ranges = IRanges(start = fusions$Pos1, width = 1), strand = fusions$Strand1)
gr2 <- GRanges(seqnames = fusions$Chr2, ranges = IRanges(start = fusions$Pos2, width = 1), strand = fusions$Strand2)

tx <- transcriptsByOverlaps(txdb, gr1)
transform(fusions, chr1 = s
