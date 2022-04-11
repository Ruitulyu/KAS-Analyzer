#!/bin/bash
# 'KAS-pipe2 complexity' was developed by Ruitu Lyu on 1-10-2022.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 complexity [ -h/--help ] [ -o prefix ] [ -l labels ] [ -k KAS-seq ]"
exampleHelp="Example: nohup KAS-pipe2 complexity -o KAS-seq_complexity -l labels.txt -k KAS-seq.txt &"
prefixHelp="-o [KAS-seq_complexity]: please input the prefix (basename) of 'KAS-pipe2 complexity' output files. REQUIRED."
labelsHelp="-l [labels.txt]: please input the text file containing the labels of KAS-seq or spKAS-seq data complexity metric. Default: basename of KAS-seq files.
Example:
WT.rep1
WT.rep2
WT.rep3
WT.rep4              ---labels.txt"
KASseqHelp="-k [KAS-seq.txt]: please input the text file containing the bam files without removing redundant reads, which are used to calculate the complexity metric including PCR Bottlenecking Coefficients (PBC) and Non-Redundant Fraction (NRF). The order and number of KAS-seq data should be the consistent with the labels file. REQUIRED.
Example:
KAS-seq.rep1.bam
KAS-seq.rep2.bam
KAS-seq.rep3.bam
KAS-seq.rep4.bam     ---KAS-seq.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 complexity' shell script is applied to calculate the library complexity metric of (sp)KAS-seq data including, the PCR Bottlenecking Coefficient and Non-Redundant Fraction (NRF) for KAS-seq. Please refer to https://www.encodeproject.org/data-standards/terms/ for more details about the library complexity metric."

printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$labelsHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 complexity' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

while getopts 'ho:l:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
	o) prefix=$OPTARG ;;
	l) labels=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing the bam files without removing redundant reads generated from 'KAS-pipe2 KAS-seq'. -k [KAS-seq.txt]"
   echo ""
   exit -1
fi

if test -z $prefix ;then
   echo ""
   echo "Please input the prefix (basename) of 'KAS-pipe2 complexity' output files. -o [prefix]."
   echo ""
   exit -1
fi

# setup the labels of spKAS-seq samples.
# get the number of samples.
number_of_samples=$( awk 'END {print NR}' $KASseq )

if test -e "$labels" ;then	
number_of_labels=$( awk 'END {print NR}' $labels ) 
   if [[ $number_of_labels != $number_of_samples ]] ;then
      echo "" 
      echo "Error: the number of labels isn't consistent with the number of samples." 
      echo ""
      exit -1
   fi
fi

if test -z $labels ;then
   for ((i=1; i<=${number_of_samples}; i++))
   do
   sample_selected=$( sed -n ''$i'p' $KASseq )
   label_basename=$( basename ${sample_selected} .bam )
   echo $label_basename >> ${prefix}.labels_basename.txt
   done
   labels="${prefix}.labels_basename.txt"
fi

echo -e "Samples\tPBC\tBottlenecking_level\tNRF\tComplexity" > ${prefix}_library_complexity_metric.txt

for ((i=1; i<=${number_of_samples}; i++))
do
sample_selected=$( sed -n ''$i'p' $KASseq )
label_selected=$( sed -n ''$i'p' $labels )
echo "Calculate library complexity metric for $sample_selected."
echo ""
samtools view -q 10 -b $sample_selected | bamToBed -i | awk '{printf("%s\t%d\t%d\n",$1,$2,$3)}' | sort | uniq -c > ${sample_selected}.unique.number.bed
M1=$( awk '$1==1 {print $0}' ${sample_selected}.unique.number.bed | wc -l )
M_Distinct=$( cat ${sample_selected}.unique.number.bed | wc -l )
M_All=$( awk 'BEGIN {sum=0} {sum+=$1} END {print sum}' ${sample_selected}.unique.number.bed )
PBC=$( awk -v x=${M1} -v y=${M_Distinct} 'BEGIN{printf "%.2f\n",x/y}' )
PBC_percent=$( awk -v x=${M1} -v y=${M_Distinct} 'BEGIN{printf "%d\n",x*100/y}' )
        
     if [ $PBC_percent -ge 70 ] ;then
	 Bottlenecking="None"
     elif [ $PBC_percent -ge 50 ] && [ $PBC_percent -lt 70 ] ;then	  
         Bottlenecking="Moderate"  
     elif [ $PBC_percent -ge 0 ] && [ $PBC_percent -lt 50 ] ;then	  
         Bottlenecking="Severe"
     fi	  

NRF=$( awk -v x=${M_Distinct} -v y=${M_All} 'BEGIN{printf "%.2f\n",x/y}' )
NRF_percent=$( awk -v x=${M_Distinct} -v y=${M_All} 'BEGIN{printf "%d\n",x*100/y}' )
       
     if [ $NRF_percent -ge 70 ] ;then
	 Complexity="Ideal"
     elif [ $NRF_percent -ge 50 ] && [ $NRF_percent -lt 70 ] ;then
         Complexity="Acceptable"
     elif [ $NRF_percent -ge 0 ] && [ $NRF_percent -lt 50 ] ;then
         Complexity="Concerning"
     fi       
echo -e "${label_selected}\t${PBC}\t${Bottlenecking}\t${NRF}\t${Complexity}" >> ${prefix}_library_complexity_metric.txt

# remove intermediate files.
rm -f ${sample_selected}.unique.number.bed

echo "done."
echo ""
done

rm -f ${prefix}.labels_basename.txt

echo "'KAS-pipe2 complexity' run successfully!"
