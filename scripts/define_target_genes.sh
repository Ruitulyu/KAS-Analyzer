#!/bin/bash
# 'KAS-pipe2 targetgenes' was developed by Ruitu Lyu on 1-22-2022.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 targetgenes [ -h/--help ] [ -o prefix ] [ -s assembly id ] [ -f features ] [ -l length ] [ -p peaks ]"
exampleHelp="Example: 
Gene features:
nohup KAS-pipe2 targetgenes -o KAS-seq_peaks_target_genes -s mm10 -f promoter -p KAS-seq_peaks.bed &
Associated genes of enhancers:
nohup KAS-pipe2 targetgenes -o KAS-seq_ss_enhancers_asso_genes -s mm10 -f enhancer -l 50000 -p KAS-seq_ss_enhancers.bed &"
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 targetgenes' output files. DEFAULT: basename of peaks file."
assemblyidHelp="-s [assembly id]: please specify the genome assembly id of KAS-seq peaks. -s [assembly id]. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
featuresHelp="-f [features]: please specify the gene feagures used to define target genes. e.g. promoter, genebody, terminator, gene or enhancer. REQUIRED."
lengthHelp="-l [length]: please specify the length cutoff (length to TSS) to define enhancer's associated genes. DEFAULT: 50000."
peaksHelp="-p [peaks]: please specify the KAS-seq peaks, R-loops or enhancer to define their target genes. REQUIRED."
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 targetgenes' shell script is applied to define target or associated genes (promoter, genebody, terminator or gene) of KAS-seq peaks, R-loops or enhancers loci."

# print help function.
printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$featuresHelp"
    echo -e ""
    echo -e "$lengthHelp"
    echo -e ""
    echo -e "$peaksHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 targetgenes' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the values of options.
while getopts 'ho:s:f:l:p:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        o) prefix=$OPTARG ;;
	s) assemblyid=$OPTARG ;;
        f) features=$OPTARG ;;
	l) length=$OPTARG ;;
	p) peaks=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $assemblyid ;then
   echo ""
   echo "Please specify the genome assembly id of KAS-seq peaks, R-loop loci or enhancers. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED. -s [assembly id]"
   echo ""
   exit -1
fi

# supported assembly id.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
   echo ""
   echo "Error: unsupported assembly id: $assemblyid ; Please specify the genome assembly id of KAS-seq peaks or R-loop loci. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -s [assembly id]"
   echo ""
   exit -1
fi

if test -z $features ;then
   echo ""
   echo "Please specify the gene feagures used to define target genes. e.g. promoter, genebody, terminator, gene or enhancer. REQUIRED. -f [features]"
   echo ""
   exit -1
fi

if [[ $features != "promoter" ]] && [[ $features != "genebody" ]] && [[ $features != "terminator" ]] && [[ $features != "gene" ]] && [[ $features != "enhancer" ]] ;then
   echo ""
   echo "Error: unsupported gene features: $features; Please specify the gene feagures. e.g. promoter, genebody, terminator, gene or enhancer. REQUIRED. -f [features]"
   echo ""
   exit -1
fi 

if test -z $peaks ;then
   echo ""
   echo "Please specify the KAS-seq peaks, R-loops or enhancer to define their target or associated genes. REQUIRED. -p [peaks]"
   echo ""
   exit -1
fi

# setup the default options.
if test -z $prefix ;then
   prefix=$( basename ${peaks} .bed )
fi

if test -z $length ;then
   length=50000
fi 

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# define the protein coding target genes using intersectBed.
if [[ $features == "promoter" ]] || [[ $features == "genebody" ]] || [[ $features == "terminator" ]] || [[ $features == "gene" ]] ;then
   echo "Define the target $features of peaks: $peaks ..." 
   echo ""
   intersectBed -a $peaks -b ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.${features}.bed -wa -wb > ${prefix}_target_${features}.txt
   echo "done."
   echo ""

elif [[ $features == "enhancer" ]] ;then
   echo "Define the gene annotation based on the $length bp distance to TSS."
   echo ""
   awk -v x=${length} '{if(($2+$3)/2-x>=0) {printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,($2+$3)/2-x,($2+$3)/2+x,$4,$5,$6)} else{printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,"0",($2+$3)/2+x,$4,$5,$6)} }' ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.promoter.bed > ${assemblyid}_Refseq.promoter.${length}.bed
   echo "done."
   echo ""

   echo "Define the associated genes of enhancers: $peaks ..."
   echo ""
   colnum=$(awk '{print NF}' $peaks | head -n 1 )
   intersectBed -a $peaks -b ${assemblyid}_Refseq.promoter.${length}.bed -wa -wb -f 0.5 -loj | cut -f$((colnum+1)),$((colnum+2)),$((colnum+3)) --complement | sort -u | sortBed -i > ${prefix}_associated_genes.txt
   rm -f ${assemblyid}_Refseq.promoter.${length}.bed
   echo "done."
   echo ""
fi

echo "'KAS-pipe2 targetgenes' run successfully!"
