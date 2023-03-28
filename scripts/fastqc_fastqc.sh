#!/bin/bash
# 'KAS-Analyzer fastqc' was developped by Ruitu Lyu on 12-10-2021.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-Analyzer fastqc [ -h/--help ] [ -t threads ] [ -c contaminants ] [ -o output dir ] [ -k KAS-seq ] "
exampleHelp="Example: nohup KAS-Analyzer fastqc -t 10 -k KAS-seq.rep1.fastq.gz,KAS-seq.rep2.fastq.gz,KAS-seq.rep3.fastq.gz &"
threadsHelp="-t [threads]: please input number of threads to be used for quality control check. Default: 1."
contaminantsHelp="-c [contaminants]: please specify a file which contains the list of contaminants (format: name[tab]sequence) to screen overrepresented sequences against. Default: no."
outputdirHelp="-o [output dir]: please specify the output directory with output files."
KASseqHelp="-k [KAS-seq]: please input the KAS-seq data that you want to know the quality control, like sequencing quality, duplicates, contaminants (adapter sequence)."
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-Analyzer fastqc' shell script is applied to check quality control and identify a potential type of problem in your KAS-seq data in non-interactive mode. It mainly invoke FASTQC, please refer to the FASTQC official website for more information."

printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$threadsHelp"
    echo -e ""
    echo -e "$contaminantsHelp"
    echo -e ""
    echo -e "$outputdirHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters or '--help' was provided, 'KAS-Analyzer fastqc' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'ht:c:o:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        t) threads=$OPTARG ;;
        c) contaminants=$OPTARG ;;
        o) outputdir=$OPTARG ;;
	k) KASseq=$OPTARG ;;
        ?) printHelpAndExit ;;
    esac
done

# Required options.
if test -z $KASseq ;then
   echo ""
   echo "please input the KAS-seq data that you want to know the quality control. -k [KAS-seq]"
   echo ""
   printHelpAndExit
fi

# setup the default parameters.
if test -z $threads ;then
   threads=1
fi

if test -z $contaminants ;then
   contaminants="off"
fi

# setup the $outputdir if specified.
if test -z $outputdir ;then
   outputdir="off"
else 
   mkdir -p $outputdir	
fi

# generate the KAS-seq sample list.
echo $KASseq > .KASseq.txt
KASseqlist=$(sed "s/,/ /g" .KASseq.txt)
rm -f .KASseq.txt

if [[ $contaminants == "off" ]] ;then
   if [[ $outputdir == "off" ]] ;then
      fastqc -t $threads $KASseqlist
   else 
      fastqc -t $threads -o $outputdir $KASseqlist
   fi

else 
   if [[ $outputdir == "off" ]] ;then
      fastqc -t $threads -c $contaminants $KASseqlist
   else      
      fastqc -t $threads -c $contaminants -o $outputdir $KASseqlist	
   fi
fi 


echo "'KAS-Analyzer fastqc' run successfully!"
