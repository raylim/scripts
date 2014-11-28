#!/usr/bin/env Rscript
# extract allele read counts for titan from bam file using specified vcf

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("TitanCNA"));
suppressPackageStartupMessages(library("VariantAnnotatoin"));

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
        make_option("--out", default = NULL, type = "character", action = "store", help ="Output file (required)"),
        make_option("--ref", default = 'hg19', type = "character", action = "store", help ="reference (default = %default)"),
        make_option("--vcf", default = NULL, type = "character", action = "store", help ="positions in vcf format (required)"))

parser <- OptionParser(usage = "%prog [options] [bam file]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need input bam\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$out)) {
    cat("Need output file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$vcf)) {
    cat("Need vcf file\n\n")
    print_help(parser);
    stop();
}

fn <- arguments$args[1]
if (!file.exists(fn)) {
    cat("bam file missing\n\n")
    print_help(parser);
    stop();
}
ifn <- paste(fn, 'bai', sep = '.')
if (!file.exists(ifn)) {
    cat("bam index missing\n\n")
    print_help(parser);
    stop();
}
extractAlleleReadCounts <- function (bamFile, bamIndex, positions, ref, outputFilename = NULL,
    pileupParam = PileupParam())
{
    vcfPosns <- readVcf(positions, ref)
    which <- GRanges(as.character(vcfPosns$CHROM), IRanges(vcfPosns$POS,
        width = 1))
    sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE),
        which = which)
    tumbamObj <- BamFile(bamFile, index = bamIndex)
    counts <- pileup(tumbamObj, scanBamParam = sbp, pileupParam = pileupParam)
    countsMerge <- xtabs(count ~ which_label + nucleotide, counts)
    label <- do.call(rbind, strsplit(rownames(countsMerge), ":"))
    posn <- do.call(rbind, strsplit(label[, 2], "-"))
    countsMerge <- cbind(data.frame(chr = label[, 1]), position = posn[,
        1], countsMerge[, 1:7])
    countMat <- data.frame(chr = vcfPosns$CHROM, position = as.numeric(vcfPosns$POS),
        ref = vcfPosns$REF, refCount = 0, Nref = vcfPosns$ALT,
        NrefCount = 0, stringsAsFactors = FALSE)
    countMat <- merge(countMat, countsMerge, by = c("chr", "position"),
        sort = FALSE, stringsAsFactors = FALSE)
    NT <- c("A", "T", "C", "G")
    for (n in 1:length(NT)) {
        indRef <- countMat$ref == NT[n]
        countMat[indRef, "refCount"] <- countMat[indRef, NT[n]]
        countMat[indRef, "NrefCount"] <- rowSums(countMat[indRef,
            NT[-n]])
    }
    countMat$chr <- gsub("chr", "", countMat$chr)
    countMat <- countMat[countMat$chr %in% c(as.character(1:22),
        "X", "Y"), ]
    countMat <- countMat[, 1:6]
    if (!is.null(outputFilename)) {
        message("extractAlleleReadCounts: writing to ", outputFilename)
        write.table(countMat, file = outputFilename, row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep = "\t")
    }
    return(countMat)
}

extractAlleleReadCounts(fn, ifn, opt$vcf, opt$ref, opt$out)
