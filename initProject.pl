#!/usr/bin/env perl

my $MAKEFILE = <<ENDL;
export REF = hg19
DUP_TYPE = markdup

#TARGETS_FILE = intervals.bed
EXOME = true

# gatk options
HARD_FILTER_SNPS = true

QSUB_PRIORITY = -500

include modules/Makefile
ENDL

system "git clone https://github.com/raylim/scripts.git";
system "git clone https://github.com/raylim/modules.git";

open OUT, ">Makefile";
print OUT $MAKEFILE;

