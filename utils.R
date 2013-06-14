#!/usr/bin/env Rscript
# Description: This script is used for storing general functions shared amongst scripts
# Authors: Fong Chun Chan <fongchunchan@gmail.com>

# Description: This function is used to get the raw reads over genomic features
# Inputs: 
#	1) features : A GRanges object
#	2) reads: A GappedAlignments object generated from readBamGappedAlignments
#	3) countMethod: Two alternative counting methods can be specified. i) countOverlaps allows for reads to be counted by than once, ii) summarizeOverlaps counts reads only once 
# Outputs:
getCounts <- function( features, reads, countMethod = 'countOverlaps' ){
	print(paste('Using the', countMethod, 'counting method'))
	if ( countMethod == 'countOverlaps' ){
		counts <- countOverlaps( features, reads )
	} else if ( countMethod == 'summarizeOverlaps' ){
		summarizedExpt <- summarizeOverlaps( features, reads )
		counts <- as.numeric( assays(summarizedExpt)$counts )
		names(counts) <- rownames(summarizedExpt)
	}
	return(counts)
}

# Description: This function is used to get RPKM values over genomic features
# Inputs: 
#	1) features : A GRanges object
#	2) featureCounts: A numeric vector containing the number of raw counts aligning to the feature
# Outputs:
getExprs <- function( features, featureCounts ){
	numBases <- sum(width(features))
	numKBases <- numBases / 1000

	millionsMapped <- sum(featureCounts) / 10^6

	rpm <- featureCounts / millionsMapped

	rpkm <- rpm / numKBases

	return( list( 'rpm' = rpm, 'rpkm' = rpkm) )
}
