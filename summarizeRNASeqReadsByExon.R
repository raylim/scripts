#!/usr/bin/env Rscript
# Description: This script is used to generate exon raw counts and expression values 
# Authors: Fong Chun Chan <fongchunchan@gmail.com>
library("GenomicFeatures")
library("Rsamtools")
library('optparse')

optionList <- list(
	make_option(c('-a', '--addChr'), action='store_true', default = FALSE, help = 'Set the flag to add chr as a prefix to each seqlevel [%default]')
	)
posArgs <- c('txdbFile', 'bamFile', 'outFile')
parser <- OptionParser(usage = paste('%prog [options]', paste(posArgs, collapse=' ')),  option_list=optionList)
arguments <- parse_args(parser, positional_arguments = TRUE)

if (length(arguments$args) != length(posArgs)) {
	print_help(parser)
	print(arguments$args)
	stop('Incorrect number of required positional arguments')
} else {
	cmdArgs <- arguments$args
	for (i in 1:length(cmdArgs)){
		assign(posArgs[i], cmdArgs[i])
	}
	opt <- arguments$options
}

#For debugging
if (FALSE){
	opt <- list('addChr' = TRUE, 'geneListFile' = NULL)
	txdbFile <- '~/ensg69.biomart.13012013.sqlite'
	bamFile <- '~/share/data/DLBCL/WTSS/bam/HS0653.bam'
	outFile <- 'tmp.txt'
}
#txdb <- makeTranscriptDbFromBiomart( biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl' )
#saveFeatures(txdb, '~/ensg69.biomart.13012013.sqlite')
#txdb <- makeTranscriptDbFromUCSC(genome = 'hg19', tablename = 'ensGene')
#
cat("Loading", txdbFile, " ... ")
txdb <- loadDb(txdbFile)
allExons <- exons(txdb, columns = c('gene_id', 'exon_id', 'exon_name'))

if ( opt$addChr ){
	print('Prefixing chr to the chromosome names ...')
	newSeqNames <- paste('chr', seqlevels(allExons), sep = '')
	names(newSeqNames) <- seqlevels(allExons)
	allExons <- renameSeqlevels( allExons, newSeqNames )
}
cat("Finished\n")

cat("Reading", bamFile, " ... ")
si <- seqinfo(BamFile(bamFile));
gr <- GRanges(seqnames(si), IRanges(100, seqlengths(si)-100));
scf <- scanBamFlag( isDuplicate = FALSE ) # remove duplicate reads
reads <- readBamGappedAlignments( bamFile, param = ScanBamParam( which = gr, flag = scf ) ); # grab reads in specific region cat("Finished\n")
#reads <- GRanges(seqnames = rname(reads), ranges = IRanges(start = start(reads), end = end(reads)), strand = rep('*', length(reads)));
cat('Finished\n')

print('Count raw exon read counts ...')
#countsForExons <- countOverlaps(allExons, reads);
summarizedExpt <- summarizeOverlaps(allExons, reads)
countsForExons <- as.numeric( assays(summarizedExpt)$counts )
names(countsForExons) <- rownames(summarizedExpt)
print('... Done')

print('Generating expression values ...')
numBases <- width(allExons)
numKBases <- numBases / 1000
millionsMapped <- sum(countsForExons) / 10^6
rpm <- countsForExons / millionsMapped
rpkm <- rpm / numKBases
print('... Done')

print('Retrieving annotation data ...')
annotDf <- values(allExons)
print('...Done')

exonsReadDf <- data.frame( geneID = sapply(annotDf[, 'gene_id'], '[[', 1), exonID = annotDf[, 'exon_id'], exonName = annotDf[, 'exon_name'], exonCount = countsForExons, exonRPM = rpm, exonRPKM = rpkm, stringsAsFactors = FALSE )

print(paste('Writing data to', outFile))
write.table(exonsReadDf, file = outFile, sep = '\t', quote = F, row.names=F)
print('...Done')
