#!/bin/bash
# 'KAS-Analyzer RNaseH' was developed by Ruitu Lyu on 12-15-2022.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-Analyzer RNaseH [ -h/--help ] [ -p threads ] [ -o prefix ] [ -r WT_R-loop.bed ] [ -f fold change ] [ -c WT_R-loop_density.txt ] [ -t RNaseH_R-loop_density.txt ]"
exampleHelp="Example: nohup KAS-Analyzer RNaseH -p 10 -o RNaseH_sensitive_R-loops -r WT_R-loop.bed -f 2 -c WT_R-loop_density.txt -t RNaseH_R-loop_density.txt &"
threadsHelp="-p [threads]: please specify the number of threads used for R-loops identification. DEFAULT: 1."
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-Analyzer RNaseH' output files. REQUIRED."
RloopsHelp="-r [WT_R-loops.bed]: please specify the R-loop list defined using spKAS-seq data in cells without RNase H treatment. REQUIRED."
foldchangeHelp="-f [fold change]: please specify the fold change threshold of R-loop density difference between WT and RNase H samples used for the identification of R-loops sensitive to RNase H treatment. DEFAULT: 2."
WTRloopHelp="-c [WT_R-loop_density.txt]: please input the text file containing WT R-loop density files. REQUIRED.
Example:
WT_R-loop.rep1.bigWig
WT_R-loop.rep2.bigWig           ---WT_R-loop_density.txt"
RNaseHRloopHelp="-t [RNaseH_R-loop_density.txt]: please input the text file containing RNase H R-loop density files. REQUIRED.
Example:
RNaseH_R-loop.rep1.bigWig
RNaseH_R-loop.rep2.bigWig       ---RNaseH_R-loop_density.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-Analyzer RNaseH' shell script is applied to identify R-loops sensitive to RNase H treatment."

printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$threadsHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$RloopsHelp"
    echo -e ""
    echo -e "$foldchangeHelp"
    echo -e ""
    echo -e "$WTRloopHelp"
    echo -e ""
    echo -e "$RNaseHRloopHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer RNaseH' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

while getopts 'hp:o:r:f:c:t:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        p) threads=$OPTARG ;;
        o) prefix=$OPTARG ;;
        r) Rloops=$OPTARG ;;
        f) foldchange=$OPTARG ;;
        c) WT=$OPTARG ;;
        t) RNaseH=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# check deeptools were installed in your system.

if ! type deeptools > /dev/null 2>&1 ;then
   echo "deeptools was not installed or not export to the \$PATH'"
   echo ""
   echo "Install deeptools with 'conda install -c bioconda deeptools' or refer the official website of 'deeptools'."
   echo ""
   exit 1
fi

# Required options.
if test -z $prefix ;then
   echo ""      
   echo "Please input the prefix (basename) of 'KAS-Analyzer RNaseH' output files. REQUIRED: -o prefix"
   echo ""
   exit 1
fi

if test -z $Rloops ;then
   echo ""
   echo "Please specify the R-loop list defined by spKAS-seq data in cells without RNase H treatment. REQUIRED: -r WT_R-loops.bed"
   echo ""
   exit 1
fi

if test -z $WT ;then
   echo ""
   echo "Please input the text file containing WT R-loop density files. REQUIRED: -c WT_R-loop_density.txt"
   echo ""
   exit 1
fi

if test -z $RNaseH ;then
   echo ""
   echo "Please input the text file containing RNase H R-loop density files. REQUIRED: -t RNaseH_R-loop_density.txt"
   echo ""
   exit 1
fi

# setup parameters of default options.
if test -z $threads ;then
   threads=1
fi

if test -z $foldchange ;then
   foldchange=2
fi

# get the absolute path of 'KAS-Analyzer R-loop' shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# get the list of KASseq bigWig files and labels.
WTRloop_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $WT)
RNaseHRloop_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $RNaseH)

echo "Calculate the averaged WT R-loop density on ${Rloops} ..."
echo ""

Rloops_basename=$( basename ${Rloops} .bed )
multiBigwigSummary BED-file --bwfiles $WTRloop_list --BED $Rloops -p $threads -out ${prefix}_WT_on_${Rloops_basename}.npz --outRawCounts ${prefix}_WT_on_${Rloops_basename}.tab > /dev/null 2>&1
sed "s/nan/0/g" ${prefix}_WT_on_${Rloops_basename}.tab | sed "1d" > ${prefix}_WT_on_${Rloops_basename}.bed

awk '{printf("%s\t%d\t%d\n",$1,$2,$3)}' $Rloops > ${Rloops_basename}.regions
awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/(NF-3) }' ${prefix}_WT_on_${Rloops_basename}.bed > ${prefix}_WT_on_${Rloops_basename}.average

echo "done."
echo ""

echo "Calculate the averaged RNaseH R-loop density on ${Rloops} ..."
echo ""

Rloops_basename=$( basename ${Rloops} .bed )
multiBigwigSummary BED-file --bwfiles $RNaseHRloop_list --BED $Rloops -p $threads -out ${prefix}_RNaseH_on_${Rloops_basename}.npz --outRawCounts ${prefix}_RNaseH_on_${Rloops_basename}.tab > /dev/null 2>&1
sed "s/nan/0/g" ${prefix}_WT_on_${Rloops_basename}.tab | sed "1d" > ${prefix}_WT_on_${Rloops_basename}.bed

awk '{printf("%s\t%d\t%d\n",$1,$2,$3)}' $Rloops > ${Rloops_basename}.regions
awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/(NF-3) }' ${prefix}_RNaseH_on_${Rloops_basename}.bed > ${prefix}_RNaseH_on_${Rloops_basename}.average

echo "done."
echo ""

echo "Output the R-loop list with averaged WT and RNaseH R-loop density and filter R-loops sensitive to RNase H treatment."
echo ""

paste ${Rloops_basename}.regions ${prefix}_WT_on_${Rloops_basename}.average ${prefix}_RNaseH_on_${Rloops_basename}.average > ${prefix}_RNaseH_vs_WT_on_${Rloops_basename}.bed
awk -v x=$foldchange '($4+0.01)/($5+0.01)>=x {printf("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n",$1,$2,$3,$4,$5,log(($5+0.01)/($4+0.01))/log(2))}' ${prefix}_RNaseH_vs_WT_on_${Rloops_basename}.bed > ${prefix}_R-loops_sensitive_to_RNaseH.bed

echo "done."
echo ""

echo "Clean up the intermediate files."
echo ""

rm -f ${prefix}_WT_on_${Rloops_basename}.tab
rm -f ${prefix}_WT_on_${Rloops_basename}.npz
rm -f ${prefix}_WT_on_${Rloops_basename}.bed
rm -f ${prefix}_RNaseH_on_${Rloops_basename}.tab
rm -f ${prefix}_RNaseH_on_${Rloops_basename}.npz
rm -f ${prefix}_RNaseH_on_${Rloops_basename}.bed
rm -f ${Rloops_basename}.regions
rm -f ${prefix}_WT_on_${Rloops_basename}.average
rm -f ${prefix}_RNaseH_on_${Rloops_basename}.average
rm -f ${prefix}_RNaseH_vs_WT_on_${Rloops_basename}.bed

echo "done."
echo ""

echo "'KAS-Analyzer RNaseH' run successfully!"
