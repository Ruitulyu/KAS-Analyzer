#!/bin/bash
# 'KAS-Analyzer bedGraphToBigWig' was developed by Ruitu Lyu on 12-10-2021.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-Analyzer ToBigWig [ -h/--help ] [ -k KAS-seq ] [ -s assembly id ]"
exampleHelp="Example: nohup KAS-Analyzer ToBigWig -k KAS-seq_data.txt -s hg19 &"
KASseqHelp="-k [KAS-seq]: please input the text file containing the bedGraph files generated from 'KAS-Analyzer KAS-seq'. REQUIRED.
Example:
KAS-seq_WT.rep1.nor.bg
KAS-seq_WT.rep2.nor.bg
KAS-seq_KO.rep1.nor.bg 
KAS-seq_KO.rep2.nor.bg    ---KAS-seq_data.txt"
assemblyidHelp="-s [assembly id]: please input the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the normalized KAS-seq bedGraph files. REQUIRED."
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-Analyzer ToBigWig' shell script is applied to convert (sp)KAS-seq bedGraph files to bigWig files."

printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer bedGraphToBigWig' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
   printHelpAndExit
fi

# get the value of options.
while getopts 'hk:s:' opt; do
    case $opt in
        h) printHelpAndExit ;;
        k) KASseq=$OPTARG ;;
        s) assemblyid=$OPTARG ;;
        ?) printHelpAndExit ;;
    esac
done

# check bedtools and bedGraphToBigWig were installed in your system.
if ! type bedtools > /dev/null 2>&1 ;then
   echo "bedtools was not installed or not export to the \$PATH'"
   echo ""
   echo "Install bedtools with 'conda install -c bioconda bedtools' or refer the official website of 'bedtools'."
   echo ""
   exit 1
fi

if ! type bedGraphToBigWig > /dev/null 2>&1 ;then
   echo "bedGraphToBigWig was not installed or not export to the \$PATH'"
   echo ""
   echo "Install bedGraphToBigWig with 'conda install -c bioconda ucsc-bedgraphtobigwig' or refer the official website of UCSC command lines."
   echo ""
   exit 1
fi

# Required options.
if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing the normalized bedGraph files. -k [KAS-seq]"
   echo ""
   exit -1
fi

# supported assembly id.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
        echo ""
        echo "Error: unsupported assembly id: $assemblyid "
        echo ""
        exit 1
fi


# get the number of samples.
number_of_samples=$( awk 'END {print NR}' $KASseq )

# get the path of 'KAS-Analyzer bedGraphToBigWig' shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)


# transfer spKAS-seq or KAS-seq bedGraph to bigWig files
for ((i=1; i<=${number_of_samples}; i++))
do
sample_selected=$(sed -n ''$i'p' $KASseq)
KASseq_basename=$(basename ${sample_selected} .bg)
echo "Transfering ${KASseq_basename}.bg into ${KASseq_basename}.bigWig ..."
echo ""
sed '/^chrM/d' $sample_selected | sortBed -i - > ${KASseq_basename}.sort.bg
bedGraphToBigWig ${KASseq_basename}.sort.bg ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes ${KASseq_basename}.bigWig
rm -f ${KASseq_basename}.sort.bg
echo "done."
echo ""
done

echo "'KAS-Analyzer ToBigWig' run successfully!"
