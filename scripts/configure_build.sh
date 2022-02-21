#!/bin/bash

# Stop on error
set -e

# Help information
usageHelp="Usage: KAS-pipe2 build [ -h ] [ -a aligner ] [ -g genome fasta ] [ -p index prefix ] [ -t threads ] [ -d index dir ]"
exampleHelp="Example: KAS-pipe2 build -a bowtie2 -g ./genome.fa -p hg19 -t 10 -d /Software/hg19_Bowtie2Index/ "
alignerHelp="-a [aligner]: please specify aligner name you want to use, e.g. bowtie, bowtie2 or bwa. REQUIRED."
fastaHelp="-g [genome fasta]: please input the path of reference genome fasta file. REQUIRED."
threadsHelp="-t [threads]: please specify the number of threads. Default: 1."
prefixHelp="-p [index prefix]: please input the prefix (basename) of the aligners' (bowtie, bowtie2 or bwa) reference genome index. Default: basename of fasta file ."
dirHelp="-d [index dir]: directory to save newly built reference genome index. REQUIRED."
helpHelp="-h\-help: print this help and exit."

# Function used to output the help information.
printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e "" 
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$alignerHelp"
    echo -e ""
    echo -e "$fastaHelp"
    echo -e ""
    echo -e "$threadsHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$dirHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 build' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'ha:g:p:t:d:' opt; do
    case $opt in
	p) prefix=$OPTARG ;;
	t) threads=$OPTARG ;;
        h) printHelpAndExit 0;;
        a) aligner=$OPTARG ;;
        g) fasta=$OPTARG ;;
        d) indexdir=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options. 
if test -z $aligner ;then
   echo ""
   echo "Please input aligner: bowtie2 or bwa."
   echo ""
   printHelpAndExit
fi

if test -z $fasta ;then
   echo ""
   echo "Please input the path of reference genome fasta file."
   echo ""
   printHelpAndExit
fi

if test -z $indexdir ;then
   echo ""
   echo "Please input the directory to save reference genome index."
   echo ""
   printHelpAndExit
fi

# setup the default parameters.
if test -z $threads ;then
   threads=1
fi

if test -z $prefix ;then
   prefix=$(basename ${fasta} .fa)
fi

# Build bowtie2 index with bowtie2-build.
if [[ "$aligner" == "bowtie2" ]]; then

   if ! type bowtie2 > /dev/null 2>&1 ;then
      echo "bowtie2 was not installed or not export to the \$PATH'"
      echo ""
      echo "Install bowtie2 with 'conda install -c bioconda bowtie2' or 'KAS-pipe2 install -t bowtie2'."
      echo ""
      exit 1
   fi            	
   
   mkdir -p $indexdir
   fasta_name=$(basename ${fasta} .fa).fa
   mv -n $fasta $indexdir
   cd $indexdir
   bowtie2-build --threads $threads $fasta_name $prefix

# Build bwa index with bwa index.
elif [[ "$aligner" == "bwa" ]]; then

   if ! type bwa > /dev/null 2>&1 ;then
   echo "bwa was not installed or not export to the \$PATH'"
   echo ""
   echo "Install bwa with 'conda install -c bioconda bwa' or 'KAS-pipe2 install -t bwa'."
   echo ""
   exit 1
   fi

   mkdir -p $indexdir
   fasta_name=$(basename ${fasta} .fa).fa
   mv -n $fasta $indexdir
   cd $indexdir
   bwa index -p $prefix $fasta_name
   mv -n $fasta_name ${prefix}.fa

else 
   echo ""	
   echo "Unsupported aligner $aligner. e.g. bowtie2 or bwa."
   echo ""
   printHelpAndExit

fi

echo "$aligner $fasta_name index built successfully."

echo "'KAS-pipe2 build' run successfully!"
