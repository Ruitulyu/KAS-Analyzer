#!/bin/bash
# 'KAS-Analyzer normalize' was developed by Ruitu Lyu on 12-10-2021.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-Analyzer normalize [ -h/--help ] [-m methods ] [ -k KAS-seq ] [ -r ratios ] [ -b ] [ -s assembly id ]"
exampleHelp="Example: nohup KAS-Analyzer normalize -m ratios -k KAS-seq_data.txt -r ratios.txt -b -s mm10 &"
methodsHelp="-m [methods]: please input the methods used for KAS-seq data normalization. e.g. ratios or RPKM. DEFAULT: ratios."
KASseqHelp="-k [KAS-seq_data.txt]: please input the text file containing the bedGraph files generated from 'KAS-Analyzer (sp)KAS-seq'. REQUIRED.
Example:
-m ratios:                            -m FPKM:
KAS-seq_WT.rep1.bg                    KAS-seq_WT.rep1.bam         
KAS-seq_WT.rep2.bg                    KAS-seq_WT.rep2.bam
KAS-seq_KO.rep1.bg                    KAS-seq_KO.rep1.bam
KAS-seq_KO.rep2.bg                    KAS-seq_KO.rep2.bam           ---KAS-seq_data.txt"
ratiosHelp="-r [ratios.txt]: please input the text file containing ratios that used to normalize KAS-seq data, which can be calculated based on mapped reads number or SpikeIn reads. The order and number of ratios should be the consistent with KAS-seq bedGraph files. REQUIRED when '-m ratios' was specified.
Example:
1.10
1.20
1.30
1.23                 ---ratios.txt"
bigWigHelp="-b: please specify if you want to convert the normalized bedGraph files into bigWig files. DEFAULT: off."
assemblyidHelp="-s [assembly id]: please input the reference genome assembly id of bedGraph files. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the normalized KAS-seq bedGraph files. REQUIRED only if -b was specified."
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-Analyzer normalize' shell script is applied to normalize spKAS-seq or KAS-seq data."

# print help function.
printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$methodsHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$ratiosHelp"
    echo -e ""
    echo -e "$bigWigHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer normalize' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'hm:k:r:bs:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
	m) methods=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        r) ratios=$OPTARG ;;
	b) bigWig="on" ;; 
	s) assemblyid=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing the bedGraph or bam files generated from 'KAS-Analyzer KAS-seq'. -k [KAS-seq.txt]"
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

# setup the genome size.
if [[ $assemblyid == "hg18" ]] ;then
    genomesize=2943582179
elif [[ $assemblyid == "hg19" ]] ;then
    genomesize=2960053816
elif [[ $assemblyid == "hg38" ]] ;then    
    genomesize=3113504866
elif [[ $assemblyid == "mm9" ]] ;then
    genomesize=2674874829
elif [[ $assemblyid == "mm10" ]] ;then    
    genomesize=2707406244
elif [[ $assemblyid == "mm39" ]] ;then
    genomesize=2709188173
elif [[ $assemblyid == "dm3" ]] ;then
    genomesize=165785384
elif [[ $assemblyid == "dm6" ]] ;then
    genomesize=145454841
elif [[ $assemblyid == "rn6" ]] ;then
    genomesize=2787632078
elif [[ $assemblyid == "rn7" ]] ;then    
    genomesize=2679540422
elif [[ $assemblyid == "ce10" ]] ;then    
    genomesize=102291793
elif [[ $assemblyid == "ce11" ]] ;then    
    genomesize=102292132
elif [[ $assemblyid == "danRer10" ]] ;then
    genomesize=1397189544
elif [[ $assemblyid == "danRer11" ]] ;then
    genomesize=1708203690
fi

# Default parameter.
if test -z $bigWig ;then
   bigWig="off"
fi   

if test -z $methods ;then
  methods="ratios"
fi

if [[ $methods == "ratios" ]] && test -z $ratios ;then
   echo ""
   echo "Please input the text file containing ratios that used to normalize KAS-seq data. -r [ratios.txt]"
   echo ""
   exit -1
fi

# get the path of 'KAS-Analyzer KAS-seq' shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# Test if the number of ratios is consistent with the number of samples.
number_of_samples=$(awk 'END {print NR}' $KASseq )

if test -n "$ratios" ;then
number_of_ratios=$(awk 'END {print NR}' $ratios )
fi

if test -n "$number_of_ratios" && [[ ${number_of_samples} != ${number_of_ratios} ]] ;then
   echo ""
   echo "Error: the number of ratios isn't consistent with the number of samples."
   echo ""
   exit
fi

# Normalize bedGraph files
if [[ $methods == "ratios" ]] ;then
echo "Normalize KAS-seq data using user provided ratios."
echo ""
   for ((i=1; i<=${number_of_samples}; i++))
   do
   sample_selected=$(sed -n ''$i'p' $KASseq)
   ratio_selected=$(sed -n ''$i'p' $ratios)
   KASseq_basename=$(basename ${sample_selected} .bg) 
   echo "Normalizing $sample_selected ..."
   echo ""
   awk -v ratio="$ratio_selected" '{printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$4*ratio)}' $sample_selected | intersectBed -a - -b ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -v > ${KASseq_basename}.nor.bg 
   echo "done."
   echo ""

   if [[ $bigWig == "on" ]] ;then 
      echo "Transfering ${KASseq_basename}.nor.bg into ${KASseq_basename}.nor.bigWig ..."
      echo ""
      sed '/^chrM/d' ${KASseq_basename}.nor.bg | sortBed -i - > ${KASseq_basename}.nor.sort.bg
      bedGraphToBigWig ${KASseq_basename}.nor.sort.bg ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes ${KASseq_basename}.nor.bigWig
      rm -f ${KASseq_basename}.nor.sort.bg
      echo "done."
      echo ""
   elif [[ $bigWig == "off" ]] ;then       
      echo "${KASseq_basename}.nor.bg will not be converted into ${KASseq_basename}.nor.bigWig." 	    
   fi 
   done

elif [[ $methods == "RPKM" ]] ;then
echo "Normalize KAS-seq data using RPKM method."
echo ""     
     if [[ $bigWig == "on" ]] ;then
     for ((i=1; i<=${number_of_samples}; i++))
     do
     sample_selected=$(sed -n ''$i'p' $KASseq)
     samtools index $sample_selected
     echo "Normalize $sample_selected using RPKM and output bigWig file ..."
     echo ""
     bamCoverage -b $sample_selected --outFileFormat bigwig -p 10 -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed --effectiveGenomeSize $genomesize --normalizeUsing RPKM -o ${sample_selected}.RPKM.bigWig > /dev/null 2>&1
     echo "done."
     echo ""
     done

     elif [[ $bigWig == "off" ]] ;then
     for ((i=1; i<=${number_of_samples}; i++))
     do
     sample_selected=$(sed -n ''$i'p' $KASseq)	     
     samtools index $sample_selected
     echo "Normalize $sample_selected using RPKM and output bedGraph file ..."    
     echo ""
     bamCoverage -b $sample_selected --outFileFormat bedgraph -p 10 -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed --effectiveGenomeSize $genomesize --normalizeUsing RPKM -o ${sample_selected}.RPKM.bg > /dev/null 2>&1
     echo "done."
     echo ""
     done
     fi
else      
echo "Error: unsupported normalization methods: $methods "
echo ""
exit

fi     

echo "'KAS-Analyzer normalize' run successfully!"
