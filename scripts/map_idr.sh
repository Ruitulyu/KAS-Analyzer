#!/bin/bash
  
# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-Analyzer idr [ -h ] [ -1 rep1 peaklist ] [ -2 rep2 peaklist ] [ -t file type ] [ -r rank ] [ -o output file ]"
exampleHelp="Example:
     narrowPeak (macs2):
     nohup KAS-Analyzer idr -1 KAS-seq.narrow_peakslist.rep1.bed -2 KAS-seq.narrow_peakslist.rep2.bed -t narrowPeak -r signal.value -o KAS-seq.narrow_peakslist &
     broadPeak (epic2):
     nohup KAS-Analyzer idr -1 KAS-seq.broad_peakslist.rep1.bed -2 KAS-seq.broad_peakslist.rep2.bed -t broadPeak -r signal.value -o KAS-seq.broad_peakslist &"
rep1peakslistHelp="-1 [rep1 peaklist]: please provide the list of KAS-seq peaks for the first replicate. REQUIRED"
rep1peakslistHelp="-2 [rep2 peaklist]: please provide the list of KAS-seq peaks for the second replicate. REQUIRED"
filetypeHelp="-t: please specify the output file type, narrowPeak or broadPeak. DEFAULT: broadPeak"
rankHelp="-r [rank]: please specify the type of signals to rank the peaks in your provided lists. signal.value, p.value or q.value. DEFAULT: signal.value"
outputHelp="-o [output file]: please specify the file where the output should be written."
helpHelp="-h: print this help and exit.
Note: The 'KAS-Analyzer idr' shell script mainly invoke the idr (Irreproducible Discovery Rate) framework is a uniﬁed approach to measure the reproducibility of ﬁndings identiﬁed from KAS-seq replicate experiments and provide highly stable thresholds based on reproducibility, please refer to https://github.com/nboley/idr for more information."
