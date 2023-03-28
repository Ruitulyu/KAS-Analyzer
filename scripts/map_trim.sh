#!/bin/bash

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-Analyzer trim [ -h ] [ -a adapter ] [ -t threads ] [ -f ] [ -q quality ] [ -l length ] [ -1 read1 ] [ -2 read2 ]"
exampleHelp="Example:
       Single-end:
       nohup KAS-Analyzer trim -a illumina -t 10 -1 KAS-seq.fastq.gz &
       Paired-end:
       nohup KAS-Analyzer trim -a illumina -t 10 -1 KAS-seq.R1.fastq.gz -2 KAS-seq.R2.fastq.gz &"
adapterHelp="-a [adapter types]: adapter sequence to be trimmed. e.g. illumina, nextera or small_rna. Hint: most of the NGS data used Illumina adapter. If not specified explicitly. KAS-Analyzer trim will auto-detect."
threadsHelp="-t [threads]: number of threads to be used for trimming. DEFAULT: 1."
fastqcHelp="-f: instruct 'KAS-Analyzer trim' to check quality control before trimming. DEFAULT: off."
qualityHelp="-q [quality]: trim low-quality ends from reads in addition to adapter removal. Default Phred score(ASCII+33): 20."
lengthHelp="-l [length]: discard reads that became shorter than length INT bp because of either quality or adapter trimming. DEFAULT: 30."
read1Help="-1 [read1]: please input single-end KAS-seq raw fastq file or read 1 of paired-end KAS-seq raw fastq files. REQUIRED."
read2Help="-2 [read2]: please input read2 of paired-end KAS-seq raw fastq files."
helpHelp="-h: print this help and exit.
Note: The 'KAS-Analyzer trim' shell script mainly invoke the trim-galore, please refer to http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/ for more information."

# print help function.
printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "" 
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$adapterHelp"
    echo -e ""
    echo -e "$threadsHelp"
    echo -e ""
    echo -e "$fastqcHelp"
    echo -e ""
    echo -e "$qualityHelp"
    echo -e ""
    echo -e "$lengthHelp"
    echo -e ""
    echo -e "$read1Help"
    echo -e ""
    echo -e "$read2Help"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer trim' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
   printHelpAndExit
fi

# get the value of options.
while getopts 'hfa:t:q:l:1:2:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        a) adapter=$OPTARG ;;
	t) threads=$OPTARG ;;
        f) fastqc="true" ;;
	q) quality=$OPTARG ;;
	l) length=$OPTARG ;;
	1) read1=$OPTARG ;;
	2) read2=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

if ! type fastqc > /dev/null 2>&1 ;then
   echo "fastqc was not installed or not export to the \$PATH'"
   echo ""
   echo "Install fastqc with 'conda install -c bioconda fastqc'."
   echo ""
   exit 1
fi


if ! type trim_galore > /dev/null 2>&1 ;then
   echo "trim_galore was not installed or not export to the \$PATH'"
   echo ""
   echo "Install trim_galore with 'conda install -c bioconda trim-galore'."
   echo ""
   exit 1
fi

# Required options.
if test -z $read1 ;then
   echo ""
   echo "please input the single-end KAS-seq fastq file or read1 of paired-end KAS-seq fastq file. Compressed .fastq.gz is accepted."
   echo ""
   printHelpAndExit 0
fi

# test the read1 .fastq.gz integrity.
if [[ ${read1##*.} = gz ]] ;then
   if gzip -t $read1 ;then
      echo "$read1 is good"
   else 
      echo "$read1 is corrupt"
      exit 1
   fi
fi

# test the read2 .fastq.gz integrity, if specified.
if test -n "$read2" && [[ ${read2##*.} = gz ]] ;then
   if gzip -t $read2 ;then
      echo "$read2 is good"
   else
      echo "$read2 is corrupt"
      exit 1
   fi
fi

# setup the default parameters.
if test -z $threads ;then
   threads=1
fi

if test -z $length ;then
   length=30
fi

if test -z $adapter ;then
   echo "Adapter is not provided, 'KAS-Analyzer trim' will autodetect the adapter."

elif test -n "$adapter" && [[ $adapter != "illumina" ]] && [[ $adapter != "nextera" ]] && [[ $adapter != "small_rna" ]] ;then
   echo ""
   echo "Error: unsupported adapter types: $adapter. Please input the adapter types: illumina, nextera or small_rna."
   echo ""
fi



# Determine the type of KAS-seq sequencing data, single or paired.
if test -z $read2 ;then
   paired_or_single_end="single"
else 
   paired_or_single_end="paired"
fi

# Quality control analysis with fastqc.
if [[ $fastqc == "true" ]] && [[ $paired_or_single_end == "single" ]] ;then
   echo ""
   echo "Quality control analysis with fastqc for single-end KAS-seq data..."
   echo ""
   fastqc -t $threads $read1
   rm -f *_fastqc.zip

elif [[ $fastqc == "true" ]] && [[ $paired_or_single_end == "paired" ]] ;then
   echo ""
   echo "Quality control analysis with fastqc for paired-end KAS-seq data..."
   echo ""
   fastqc -t $threads $read1
   fastqc -t $threads $read2
   rm -f *_fastqc.zip

else
   echo "Quality control of KAS-seq raw fastq files was skipped..."

fi

# Perform quality and adapter trimming with Trim Galore.
echo ""
echo "Automate quality and adapter trimming with Trim Galore! ... "
echo ""

if [[ $paired_or_single_end == "single" ]] ;then
   if test -z $adapter ;then
      trim_galore -j $threads --fastqc --length $length $read1 
   else 
      trim_galore --${adapter} -j $threads --fastqc --length $length $read1
   fi

elif [[ $paired_or_single_end == "paired" ]]; then	
   if test -z $adapter ;then
      trim_galore -j $threads --fastqc --paired --length $length $read1 $read2    
   else
      trim_galore --${adapter} -j $threads --fastqc --paired --length $length $read1 $read2
   fi
fi
 
# Move the fastqc .html output to the fastqc directory; delete the report.txt and .zip file.
if test -e fastqc ;then
   rm -f *fastqc.zip
   cd fastqc
   mv ../*html ./
   cd ..

else
   mkdir -p fastqc
   rm -f *fastqc.zip
   cd fastqc
   mv ../*html ./
   cd ..
fi


mkdir -p Summary

if [[ $paired_or_single_end == "single" ]] ;then
   cd Summary
   mv ../${read1}_trimming_report.txt ./
   echo "Reads number before trimming" >> ${read1}_trimming_summary.txt
   grep "sequences processed in total" ${read1}_trimming_report.txt | awk '{print $1}' >> ${read1}_trimming_summary.txt
   echo "Trimmed reads number" >> ${read1}_trimming_summary.txt
   grep "Sequences removed because they became shorter than the length cutoff" ${read1}_trimming_report.txt | awk '{print $14}' >> ${read1}_trimming_summary.txt
   rm -f ${read1}_trimming_report.txt	

elif [[ $paired_or_single_end == "paired" ]] ;then
   cd Summary
   mv ../${read1}_trimming_report.txt ./	
   mv ../${read2}_trimming_report.txt ./
   echo "Reads number before trimming" >> ${read1}_trimming_summary.txt
   grep "sequences processed in total" ${read2}_trimming_report.txt | awk '{print $1}' >> ${read1}_trimming_summary.txt
   echo "Trimmed reads number" >> ${read1}_trimming_summary.txt 
   grep "Number of sequence pairs removed because at least one read was shorter than the length cutoff" ${read2}_trimming_report.txt | awk '{print $19}' >> ${read1}_trimming_summary.txt
   rm -f ${read1}_trimming_report.txt
   rm -f ${read2}_trimming_report.txt
fi  

echo "'KAS-Analyzer trim' run successfully!"
