#!/bin/bash
# 'KAS-pipe2 normalize' was developed by Ruitu Lyu on 12-10-2021.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-pipe2 normalize [ -h/--help ] [ -k KAS-seq ] [ -r ratios ]"
exampleHelp="Example: nohup KAS-pipe2 normalize -k KAS-seq_data.txt -r ratios.txt &"
KASseqHelp="-k [KAS-seq_data.txt]: please input the text file containing the bedGraph files generated from 'KAS-pipe2 (sp)KAS-seq'. REQUIRED.
Example:
KAS-seq_WT.rep1.bg
KAS-seq_WT.rep2.bg
KAS-seq_KO.rep1.bg 
KAS-seq_KO.rep2.bg   ---KAS-seq_data.txt"
ratiosHelp="-r [ratios.txt]: please input the text file containing ratios that used to normalize KAS-seq data, which can be calculated based on mapped reads number or SpikeIn reads. The order and number of ratios should be the consistent with KAS-seq bedGraph files. REQUIRED.
Example:
1.10
1.20
1.30
1.23                 ---ratios.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 normalize' shell script is applied to normalize spKAS-seq or KAS-seq data."

# print help function.
printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$ratiosHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 normalize' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'hk:r:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        k) KASseq=$OPTARG ;;
        r) ratios=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing the bedGraph files generated from 'KAS-pipe2 KAS-seq'. -k [KAS-seq.txt]"
   echo ""
   exit -1
fi

if test -z $ratios ;then
   echo ""	 
   echo "Please input the text file containing ratios that used to normalize KAS-seq data. -r [ratios.txt]"
   echo ""
   exit -1
fi

# Test if the number of ratios is consistent with the number of samples.
number_of_samples=$(awk 'END {print NR}' $KASseq )
number_of_ratios=$(awk 'END {print NR}' $ratios )

if [[ ${number_of_samples} != ${number_of_ratios} ]] ;then
   echo ""
   echo "Error: the number of ratios isn't consistent with the number of samples."
   echo ""
   exit
fi

# Normalize bedGraph files
for ((i=1; i<=${number_of_samples}; i++))
do
    sample_selected=$(sed -n ''$i'p' $KASseq)
    ratio_selected=$(sed -n ''$i'p' $ratios)
    KASseq_basename=$(basename ${sample_selected} .bg) 
    echo "Normalizing $sample_selected ..."
    echo ""
    awk -v ratio="$ratio_selected" '{printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$4*ratio)}' $sample_selected > ${KASseq_basename}.nor.bg 
    echo "done."
    echo ""
done

echo "'KAS-pipe2 normalize' run successfully!"

