
segData <- read.table(fn, sep = '\t', header = T)
segRle <- rle(segData$Segmented)
segData <- transform(segData, segId = as.factor(rep(1:length(segRle$values), segRle$lengths)))
segData <- transform(segData, length = End - Start)


chrom <- tapply(segData$Chrom, segData$segId, function(x) x[1])
start <- tapply(segData$Start, segData$segId, function(x) x[1])
end <- tapply(segData$End, segData$segId, function(x) x[length(x)])
effSegLen <- tapply(segData$length, segData$segId, sum)
normRatio <- tapply(segData$Segmented, segData$segId, function(x) x[1])

Data <- data.frame(chrom = chrom, loc.start = start, loc.end = end, eff.seg.len = effSegLen, normalized.ratio = normRatio)


