#!/bin/bash
# 'KAS-Analyzer readsnum' was developped by Ruitu Lyu on 1-18-2022.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-Analyzer readsnum [ -h/--help ] [ -o prefix ] [ -f format ] "
exampleHelp="Example: nohup KAS-Analyzer readsnum -o KAS-seq_reads_num -f fastq.gz &"
prefixHelp="-o [prefix]: please specify the prefix (basename) of 'KAS-Analyzer readsnum' output files. REQUIRED."
formatHelp="-f [format]: please specify the format of raw reads data. e.g. fastq, fq, fastq.gz, fasta, fa or fa.gz. REQUIRED."
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-Analyzer readsnum' shell script is applied to calculate the reads number of raw sequencing files."

printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$formatHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer readsnum' will print the help.
if [[ $# == 1 ]] || [[ "$1" == "--help" ]] || [[ $1 == "-help" ]];then
    printHelpAndExit
fi

while getopts 'ho:f:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        o) prefix=$OPTARG ;;
        f) format=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $prefix ;then
   echo ""
   echo "Please input the prefix (basename) of 'KAS-Analyzer readsnum' output files. -o [prefix]"
   echo ""
   exit 1
fi

if test -z $format ;then
   echo ""	
   echo "Please specify the format of raw reads data. e.g. fastq, fq, fastq.gz, fasta, fa or fa.gz. -f [format]"
   echo ""
   exit 1
fi

# Test the supported format of raw data files.
if [[ $format != "fastq" ]] && [[ $format != "fq" ]] && [[ $format != "fastq.gz" ]] && [[ $format != "fasta" ]] && [[ $format != "fa" ]] && [[ $format != "fa.gz" ]]  ;then
   echo ""	
   echo "Error: unsupported formats of raw data files: $format."
   echo ""
   exit 1
fi

if [[ $format == "fastq.gz" ]]; then
   echo -e "samples\treads_num" > ${prefix}_reads_num.txt
   for file in ./*.gz;
   do
   echo "Calculating the reads number of $file ..."
   rawdata_prefix=$(basename $file .fastq.gz)
   reads_num=$(zcat $file | wc -l | awk '{print $1/4}')

   echo -e "$rawdata_prefix\t$reads_num" >> ${prefix}_reads_num.txt
   echo "done."
   echo ""
   done

elif [[ $format == "fastq" ]]; then
   echo -e "samples\treads_num" > ${prefix}_reads_num.txt
   for file in ./*.fastq;
   do
   echo "Calculating the reads number of $file ..."
   rawdata_prefix=$(basename $file .fastq)
   reads_num=$(wc -l $file | awk '{print $1/4}')
   echo -e "$rawdata_prefix\t$reads_num" >> ${prefix}_reads_num.txt
   echo "done."
   echo ""
   done

elif [[ $format == "fq" ]]; then
   echo -e "samples\treads_num" > ${prefix}_reads_num.txt
   for file in ./*.fq;
   do
   echo "Calculating the reads number of $file ..."
   rawdata_prefix=$(basename $file .fq)
   reads_num=$(wc -l $file | awk '{print $1/4}')
   echo -e "$rawdata_prefix\t$reads_num" >> ${prefix}_reads_num.txt
   echo "done."
   echo ""
   done


elif [ "${rawdata_format}" == "fa.gz" ]; then
   echo -e "samples\treads_num" > ${prefix}_reads_num.txt
   for file in ./*.fa.gz;
   do
   echo "Calculating the reads number of $file ..."
   rawdata_prefix=$(basename $file .fa.gz)
   reads_num=$(zcat $file | wc -l | awk '{print $1/4}')
   echo -e "$rawdata_prefix\t$reads_num" >> ${prefix}_reads_num.txt
   echo "done."
   echo ""
   done

elif [ "${rawdata_format}" == "fa" ]; then
   echo -e "samples\treads_num" > ${prefix}_reads_num.txt
   for file in ./*.fa;
   do
   echo "Calculating the reads number of $file ..."
   rawdata_prefix=$(basename $file .fa)
   reads_num=$(wc -l $file | awk '{print $1/2}')
   echo -e "$rawdata_prefix\t$reads_num" >> ${prefix}_reads_num.txt
   done

elif [ "${rawdata_format}" == "fasta" ]; then
   echo -e "samples\treads_num" > ${prefix}_reads_num.txt
   for file in ./*.fasta;
   do
   echo "Calculating the reads number of $file ..."
   rawdata_prefix=$(basename $file .fasta)
   reads_num=$(wc -l $file | awk '{print $1/2}')
   echo -e "$rawdata_prefix\t$reads_num" >> ${prefix}_reads_num.txt
   done

fi

echo "'KAS-Analyzer readsnum' run successfully!"
