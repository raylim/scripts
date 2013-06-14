#!/bin/bash

# $Id: $
find -L /projects/analysis /projects/seq_analysis /archive/analysis* /archive/solexa* -maxdepth 2 -name 'HS*' -or -name 'A*' > $1


