#!/bin/bash
# 'KAS-Analyzer fingerprint' was developed by Ruitu Lyu on 12-11-2021.

# Stop on error
set -e

# arguments                                                     
usageHelp="Usage: KAS-Analyzer fingerprint [ -h/--help ] [ -t threads ] [ -s assembly id ] [ -o prefix ] [ -l labels ] [ -k KAS-seq ] "
exampleHelp="Example: nohup KAS-Analyzer fingerprint -t 10 -s hg19 -o KAS-seq_fingerprint -l labels.txt -k KAS-seq_data.txt &"
threadsHelp="-t [threads]: please input the number of threads used for generating KAS-seq fingerprint plot. Default: 1."
assemblyidHelp="-s [assembly id]: please specify the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. REQUIRED."
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-Analyzer fingerprint' output files. REQUIRED"
labelsHelp="-l [labels] please input the text file containing the labels of KAS-seq or spKAS-seq data that show in fingerprint plot. REQUIRED.
Example:
KAS-seq.rep1
KAS-seq.rep2
Input.rep1
Input.rep2                       ---labels.txt"
KASseqHelp="-k [KAS-seq]: please input the text file containing the bam files of (sp)KAS-seq data. REQUIRED.
Example:
KAS-seq.rep1.bam
KAS-seq.rep2.bam
KAS-seq_Input.rep1.bam
KAS-seq_Input.rep2.bam           ---KAS-seq.txt"
helpHelp="-h/--help: print this help and exit.
Note: the 'KAS-Analyzer fingerprint' shell script is applied to generate the fingerprint plot of (sp)KAS-seq data. For the more details about fingerprint plot, please refer to https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html."

printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$threadsHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
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

# if no parameters was provided, 'KAS-Analyzer fingerprint' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'ht:s:o:l:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        t) threads=$OPTARG ;;
        s) assemblyid=$OPTARG ;;
        o) prefix=$OPTARG ;;
        l) labels=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# check if deeptools was installed.
if ! type deeptools > /dev/null 2>&1 ;then
   echo "deeptools was not installed or not export to the \$PATH'"
   echo ""
   echo "Install deeptools with 'conda install -c bioconda deeptools' or refer the official website of 'deeptools'."
   echo ""
   exit 1
fi

# Required options.
if test -z $assemblyid ;then
   echo ""      
   echo "Please specify the reference genome assembly id of your KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -s [assembly id]"
   echo ""
   printHelpAndExit 0
fi

# supported assembly id.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
   echo ""
   echo "Error: unsupported assembly id: $assemblyid "
   echo ""
   exit 1
fi

if test -z $prefix ;then
   echo ""
   echo "Please input the prefix (basename) of 'KAS-Analyzer fingerprint' output files. -o [prefix]"
   echo ""
   printHelpAndExit 0
fi

if test -z $labels ;then
   echo ""
   echo "Please input the text file containing the labels of KAS-seq or spKAS-seq data that show in fingerprint plot. -l [labels.txt]"
   echo ""
   printHelpAndExit 0
fi

if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing the bigWig files of KAS-seq or spKAS-seq data that used to generate fingerprint plot. -k [KAS-seq]"
   echo ""
   printHelpAndExit 0
fi

# setup the default parameters.
if test -z $threads ;then
   threads=1
fi

# test if the number of samples is same to the number of labels.
number_of_samples=$( awk 'END {print NR}' $KASseq )
number_of_labels=$( awk 'END {print NR}' $labels )

if [[ $number_of_labels != $number_of_samples ]] ;then
   echo ""
   echo "Error: the number of labels isn't consistent with the number of samples."
   echo ""
   printHelpAndExit 0
fi

# get the list of KAS-seq samples and labels.
KASseq_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $KASseq)
labels_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $labels)

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# generate the index for bam files.
echo "Generate new index of bam files."
for i in $KASseq_list; do samtools index $i ; done;
echo "done."
echo ""

# plot the fingerprint.
echo "Generating fingerprint plot of svg format ..."
# plotFingerprint -b $KASseq_list --labels $labels_list --minMappingQuality 10 --skipZeros --ignoreDuplicates --region 1 --numberOfSamples 500000 -T "${prefix}_fingerprint_plot" --plotFile ${prefix}_fingerprint_plot.svg --plotFileFormat svg --outRawCounts ${prefix}_fingerprint_plot.tab --numberOfProcessors $threads --blackListFileName ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed
# rm -f ${prefix}_fingerprint_plot.tab
echo "done."
echo ""

echo "Generating fingerprint plot of png format."
plotFingerprint -b $KASseq_list --labels $labels_list --minMappingQuality 10 --skipZeros --ignoreDuplicates --region 1 --numberOfSamples 500000 -T "${prefix}_fingerprint_plot" --plotFile ${prefix}_fingerprint_plot.png --plotFileFormat png --outRawCounts ${prefix}_fingerprint_plot.tab --numberOfProcessors $threads --blackListFileName ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed
# rm -f ${prefix}_fingerprint_plot.tab
echo "done."
echo ""

echo "'KAS-Analyzer fingerprint' run successfully!"
