#!/bin/bash
# 'KAS-Analyzer motif' was developed by Ruitu Lyu on 12-22-2021.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-Analyzer motif [ -h/--help ] [ -t threads ] [ -o prefix ] [ -s assembly id ] [ -e enhancer position file ] [ -c control position file ] "
exampleHelp="Example: nohup KAS-Analyzer motif -o KAS-seq_enhancers_motifs -t 10 -s mm10 -e enhancers.bed -c control_background_peaks.bed &"
threadsHelp="-t [threads]: please specify the number of threads used for enriched TF motifs on transcribing enhancers. DEFAULT: 1."
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-Analyzer motif' output files. Default: basename of enhancers file."
assemblyidHelp="-s [assembly id]: please specify the genome assembly id of enhancers regulatory elements. -s [assembly id]. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
enhancersHelp="-e [enhaner]: please specify the enhancers position file in bed format. REQUIRED."
controlHelp="-c [control]: please specify the control position file in bed format, which was used as background for enriched TF motifs idenfication on enhancers. OPTIONAL"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-Analyzer motif' shell script is applied to identify enriched TF motifs on transcribing enhancers identified by KAS-seq data. For more detials, please refer to http://homer.ucsd.edu/homer/ngs/peakMotifs.html"

# print help function.
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
    echo -e "$enhancersHelp"
    echo -e ""
    echo -e "$controlHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer motif' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'ht:o:s:e:c:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        t) threads=$OPTARG ;;
	o) prefix=$OPTARG ;;
        s) assemblyid=$OPTARG ;;
        e) enhancers=$OPTARG ;;
        c) control=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $assemblyid ;then
   echo ""
   echo "Please specify the reference genome assembly id of enhancers regulatory elements. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -s [assembly id]"
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

if test -z $enhancers ;then
    echo ""
    echo "Please specify the enhancers position file in bed format used for TF motif identification. REQUIRED. -e [enhancer]"
    echo ""
    printHelpAndExit 0
fi

# setup the default of $prefix.
if test -z $prefix ;then
    prefix=$( basename ${enhancers} .bed )
fi

# setup the default of $threads.
if test -z $threads ;then
    $threads=1
fi

# test if homer was installed.
if ! type findMotifs.pl > /dev/null 2>&1 ;then
   echo "homer was not installed or not export to the \$PATH'"
   echo ""
   echo "Install homer findMotifs.pl with 'conda install -c bioconda homer' or 'KAS-Analyzer install -t homer'."
   echo ""
   exit 1
fi

if ! type bedtools > /dev/null 2>&1 ;then
   echo "bedtools was not installed or not export to the \$PATH'"
   echo ""
   echo "Install bedtools with 'conda install -c bioconda bedtools' or 'KAS-Analyzer install -t bedtools'."
   echo ""
   exit 1
fi

# the path of KAS-Analyzer KAS-seq shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

if test -z $control ;then
   echo "Extract sequences from $assemblyid genome for randomly selected 2000 enhancers defined in $enhancers."
   echo ""
   cat $enhancers | shuf | head -n 2000 | sortBed -i | awk '{printf("%s\t%d\t%d\t%s\n",$1,$2,$3,"enhancer"FNR)}' > ${enhancers}.2000.4bed
   bedtools getfasta -fi ${SH_SCRIPT_DIR}/../genome_index/${assemblyid}_BWAIndex/${assemblyid}.fa -bed ${enhancers}.2000.4bed -name > ${enhancers}.fa
   echo "done."
   echo ""

   echo "Starting to identify enriched motifs with HOMER from ${enhancers}.fa file."
   echo ""
   findMotifs.pl ${enhancers}.fa fasta ${prefix}_MotifOutput/ -len 8,10,12,14 -p $threads -mask > /dev/null 2>&1
   echo "done."
   echo ""

   echo "clean up intermediate files."
   echo ""
   rm -f ${enhancers}.2000.4bed
   rm -f ${enhancers}.fa
   echo "done."
   echo ""

else
   echo "Extract sequences from $assemblyid genome for randomly selected 2000 enhancers defined in $enhancers."
   echo ""
   cat $enhancers | shuf | head -n 2000 | sortBed -i | awk '{printf("%s\t%d\t%d\t%s\n",$1,$2,$3,"enhancer"FNR)}' $enhancers > ${enhancers}.2000.4bed
   bedtools getfasta -fi ${SH_SCRIPT_DIR}/../genome_index/${assemblyid}/${assemblyid}_BWAIndex/${assemblyid}.fa -bed ${enhancers}.2000.4bed -name > ${enhancers}.fa
   echo "done."
   echo ""

   echo "Extract sequences from $assemblyid genome for randomly selected 2000 background peaks defined in $control."
   echo ""
   cat $control | shuf | head -n 2000 | sortBed -i | awk '{printf("%s\t%d\t%d\t%s\n",$1,$2,$3,"control"FNR)}' > ${control}.2000.4bed
   bedtools getfasta -fi ${SH_SCRIPT_DIR}/../genome_index/${assemblyid}/${assemblyid}_BWAIndex/${assemblyid}.fa -bed ${control}.2000.4bed -name > ${control}.fa
   echo "done."
   echo ""

   echo "Starting to identify enriched motifs with HOMER from ${enhancers}.fa enhancers sequence file compared to ${control}.fa background sequence file."
   echo ""
   findMotifs.pl ${enhancers}.fa fasta ${prefix}_MotifOutput/ -fasta ${control}.fa -len 8,10,12,14 -p $threads -mask > /dev/null 2>&1
   echo "done."
   echo ""

   echo "clean up intermediate files."
   echo ""
   rm -f ${enhancers}.2000.4bed
   rm -f ${enhancers}.fa
   rm -f ${control}.2000.4bed
   rm -f ${control}.fa
   echo "done."
   echo ""

fi

echo "'KAS-Analyzer motif' run successfully!"
