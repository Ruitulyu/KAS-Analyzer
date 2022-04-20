#!/bin/bash
# 'KAS-pipe2 FRiP' was developed by Ruitu Lyu on 12-11-2021.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 FRiP [ -h/--help ] [ -o prefix ] [ -p peaks ] [ -l labels ] [ -k KAS-seq ]"
exampleHelp="Example: nohup KAS-pipe2 FRiP -o KAS-seq_FRiP -p peaks_files.txt -l labels.txt -k KAS-seq.txt &"
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 FRiP' output files. REQUIRED."
peaksHelp="-p [peaks]: please input the text file containing the peaks files. REQUIRED.
Example:
KAS-seq_WT_rep1_peaks.bed
KAS-seq_WT_rep2_peaks.bed
KAS-seq_KO_rep1_peaks.bed
KAS-seq_KO_rep2_peaks.bed      ---peaks_files.txt"
labelsHelp="-l [labels]: please input the text file containing the labels of KAS-seq or spKAS-seq data that show in fraction of reads in peaks (FRiP) plot. Default: basename of KAS-seq files.
Example:
WT_rep1
WT.rep2
KO.rep1
KO.rep2                        ---labels.txt"
KASseqHelp="-k [KAS-seq]: please input the text file containing bed files (uniquely mapped reads used for 'KAS-pipe2 peakcalling'), which are used to calcuate fraction of reads in peaks (FRiP) score. The order and number of (sp)KAS-seq data should be the consistent with the labels file. REQUIRED.
Example:
KAS-seq_WT_rep1.bed
KAS-seq_WT_rep2.bed
KAS-seq_KO_rep1.bed
KAS-seq_KO_rep2.bed            ---KAS-seq.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 FRiP' shell script is applied to calculate and plot fraction of reads in peaks (FRiP) scores."

printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$peaksHelp"
    echo -e ""
    echo -e "$labelsHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 FRiP' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
   printHelpAndExit
fi

# get the value of options.
while getopts 'ho:p:l:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        o) prefix=$OPTARG ;;
        p) peaks=$OPTARG ;;
        l) labels=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $prefix ;then
   echo ""
   echo "Please input the prefix (basename) of 'KAS-pipe2 FRiP' output files. -o [prefix]"
   echo ""
   exit -1
fi

if test -z $peaks ;then
   echo ""
   echo "Please input the text file containing the peaks files. -p [peaks]"
   echo ""
   exit -1
fi

if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing bed files (uniquely mapped reads used for 'KAS-pipe2 peakcalling'). -k [KAS-seq]"
   echo ""
   exit -1
fi

# setup the default parameters
# test if the number of labels is same to the number of samples.
number_of_samples=$( awk 'END {print NR}' $KASseq )
number_of_peaks=$( awk 'END {print NR}' $peaks )

if [[ $number_of_peaks != $number_of_samples ]] ;then
   echo ""
   echo "Error:the number of peak files isn't consistent with the number of samples."
   echo ""
   exit -1
fi

if test -n "$labels" ;then
number_of_labels=$( awk 'END {print NR}' $labels )
   if [[ $number_of_labels != $number_of_samples ]] ;then 
      echo "" 
      echo "Error: the number of labels isn't consistent with the number of samples." 
      echo ""
      exit -1
   fi
fi

# setup the default parameter of $labels.
if test -z $labels ;then
for ((i=1; i<=${number_of_samples}; i++))
do
sample_selected=$( sed -n ''$i'p' $KASseq )
label_basename=$( basename ${sample_selected} .bed )
echo $label_basename >> ${prefix}.labels_basename.txt
done
labels="${prefix}.labels_basename.txt"
fi

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# output the header of 'KAS-pipe2 FRiP' output files.
echo -e "Types\tlabels\tPercentage" > ${prefix}_FRiP.txt

# calculate the number of reads on peaks(FRiP) for every sample.
for ((j=1; j<=${number_of_samples}; j++))
do
sample_selected=$( sed -n ''$j'p' $KASseq )
peak_selected=$( sed -n ''$j'p' $peaks )
label_selected=$( sed -n ''$j'p' $labels )
echo "Calculate the number of reads on peaks for $sample_selected ..."
readsnum_on_peaks=$( intersectBed -a $sample_selected -b $peak_selected -wa | sort -u | wc -l )
readsnum=$( wc -l $sample_selected | awk '{print $1}' )

Inside_reads_percentage=$( awk -v x=${readsnum_on_peaks} -v y=${readsnum} 'BEGIN{printf "%.2f\n",x*100/y}' )
Outside_reads_percentage=$( awk -v x=${readsnum_on_peaks} -v y=${readsnum} 'BEGIN{printf "%.2f\n",100-x*100/y}' )

echo -e "Outside\t${label_selected}\t${Outside_reads_percentage}" >> ${prefix}_FRiP.txt
echo -e "Inside\t${label_selected}\t${Inside_reads_percentage}" >> ${prefix}_FRiP.txt
echo "done."
echo ""
done

rm -f ${prefix}.labels_basename.txt

echo "Plot fraction of reads in peaks (FRiP) scores of KAS-seq data."
Rscript --vanilla ${SH_SCRIPT_DIR}/../R/FRiP_barplot.R ${prefix}_FRiP.txt
echo "done."
echo ""

mv KAS-seq_fraction_of_reads_in_peaks.png ${prefix}_FRiP.png
mv KAS-seq_fraction_of_reads_in_peaks.svg ${prefix}_FRiP.svg
# rm -f ${prefix}_FRiP.txt

echo "'KAS-pipe2 FRiP' run successfully!"
