#!/usr/bin/env Rscript
# This script calculates RPKM values.  
# Author: Fong Chun Chan <fongchunchan@gmail.com>
library('optparse')
library("GenomicFeatures")

optionList <- list()
posArgs <- c('txdbFile', 'summarized_read_file', 'rpkm_file')
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

#Debugging#
if (FALSE){
	txdbFile <- '~/hg19_ensGene.06022012.sqlite'
	summarized_read_file <- '/share/lustre/gascoyne/DLBCL/WTSS/summarized_reads/HS0639.gsnap.summarized_reads.txt'
	rpkm_file <- '/share/lustre/gascoyne/DLBCL/WTSS/expression/rpkm/HS0639.gsnap.rpkm.txt'
}

cat("Loading", txdbFile, " ... ")
txdb <- loadFeatures(txdbFile)
txByGene <- transcriptsBy(txdb, 'gene');
cat("Finished\n")

print('Loading summarized reads file')
sumReadsDf <- read.table(summarized_read_file, sep='\t', as.is=T, header=T)
genes <- sumReadsDf[, 'gene_id']
numOfAlignedReads = sum(sumReadsDf$summarized_reads)/1e+6

rpkmDf <- data.frame(matrix(NA, length(genes), 2, dimnames=list(NULL, c('gene_id', 'rpkm'))))
for(i in 1:length(genes)){
	gene <- genes[i]
	print(paste('Calculating RPKM for gene -', gene))

	geneLength <- max(end(ranges(txByGene[[gene]]))) - min(start(ranges(txByGene[[gene]])))
	geneLengthInKB <- geneLength/1000
	counts <- sumReadsDf[sumReadsDf$gene_id==gene, 'summarized_reads']
	rpm <- counts/numOfAlignedReads
	rpkm <- rpm/geneLengthInKB
	rpkmDf[i, ] <- c(gene, rpkm)
}
print(paste('Writing RPKM results to', rpkm_file))
write.table(rpkmDf, rpkm_file, sep='\t', quote=F, row.names=F)
