disease <- "breastcancer"
platform <- "SNP_6.0"
#seg.dat.fn
#tumour <- ""
#normal <- ""
#results.dir <- ""

args <- (commandArgs(TRUE))

if (length(args) == 0) {
	print("No arguments supplied.\n")
    exit(1)
} else {
	for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
    }
}

if (tumour == "") {
	print("Need tumour sample name\n")
    exit(1)
}
if (results.dir == "") {
	print("Need results directory\n")
    exit(1)
}
if (out.file == "") {
	print("Need output file\n")
    exit(1)
}

library(ABSOLUTE)
RunAbsolute(seg.dat.fn = seg.data.fn, out.file = out.file,
    sigma.p=0, max.sigma.h=0.02,
    min.ploidy=0.95, max.ploidy=10, primary.disease=disease,
    platform=platform, sample.name=tumour,
    results.dir=results.dir,
    max.as.seg.count=1500, copy_num_type="allelic",
    max.neg.genome=0, max.non.clonal=0,
    verbose=TRUE)
