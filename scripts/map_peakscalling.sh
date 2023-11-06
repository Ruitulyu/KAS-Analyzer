#!/bin/bash
# 'KAS-Analyzer peakscalling' was developed by Ruitu Lyu on 12-09-2021.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-Analyzer peakscalling [ -h ] [ -m peaks caller ] [ -k KAS-seq ] [ -c Control ] [ -b mode ] [ -o prefix ] [ -p pvalue or qvalue ] [ -g assembly id ]."
exampleHelp="Example: nohup KAS-Analyzer peakscalling -k KAS-seq.rep1.bed,KAS-seq.rep2.bed -c Control_Input.rep1.bed,Control_Input.rep2.bed -o KAS-seq -g hg19 &"
peakscallerHelp="-m [peaks caller]: please input the peaks caller (macs2, epic2 or macs2_and_epic2) that you want to use for KAS-seq peaks calling. DEFAULT: macs2_and_epic2."
KASseqHelp="-k [KAS-seq]: please input the KAS-seq bed or bam files. e.g. KAS-seq.rep1.bed,KAS-seq.rep2.bed or KAS-seq.rep1.bam,KAS-seq.rep2.bam. REQUIRED."
controlHelp="-c [Control]: please input the KAS-seq control bed or bam files. e.g. KAS-seq_Input.rep1.bed,KAS-seq_Input.rep2.bed or KAS-seq_Input.rep1.bam,KAS-seq_Input.rep2.bam. OPTIONAL."
modeHelp="-b [mode]: specify macs2 to perferm KAS-seq peaks calling with 'broad', 'sharp' or 'both' mode. epic2 for broad peaks; macs2 for sharp peaks; combined epic2&macs2 for broad and sharp peaks . DEFAULT: both."
prefixHelp="-o [prefix]: please input the prefix (basename), which will be used to generate the name of 'KAS-Analyzer peakscalling' output files. REQUIRED."
cutoffHelp="-p [FDR or qvalue]: please input the pvalue or qvalue for KAS-seq peaks calling with macs14 or macs2. DEFAULT: epic2 FDR: 0.05; macs2 qvalue: 0.01"
assemblyidHelp="-g [assembly id]: please specify the reference genome assembly id of KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
# genome size. e.g. human(hs): 2.7e9; mouse(mm): 1.87e9; C.elegans(ce): 9e7; fruitfly(dm): 1.2e8; rat(rn): 2.5e9; zebrafish(danRer): 1e9.
helpHelp="-h: print this help and exit.
Note: This shell script mainly invoke macs2 and epic2 for calling (sp)KAS-seq data peaks, please google their github pages for more information."

# print help function.
printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$peakscallerHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$controlHelp"
    echo -e ""
    echo -e "$modeHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$cutoffHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer peakscalling' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'hm:k:c:b:o:p:g:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        m) peakscaller=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        c) control=$OPTARG ;;
	b) mode=$OPTARG ;;
        o) prefix=$OPTARG ;;
        p) cutoff=$OPTARG ;;
        g) assemblyid=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# check if macs2\epic2\macs14 was installed.
if ! type macs2 > /dev/null 2>&1 ;then
   echo "macs2 was not installed or not export to the \$PATH'"
   echo ""
   echo "Install macs2 with 'conda install -c bioconda macs2' or refer the official website of 'macs2'."
   echo ""
   exit 1
fi

if ! type epic2 > /dev/null 2>&1 ;then
   echo "epic2 was not installed or not export to the \$PATH'"
   echo ""
   echo "Install epic2 with 'conda install -c bioconda epic2=0.0.52' or refer the official website of 'epic2'."
   echo ""
   exit 1
fi

# Required options.
if test -z $KASseq ;then
   echo ""
   echo "Please input the KAS-seq file. e.g. KAS-seq.bed or KAS-seq.bam. -t [KAS-seq]"
   echo ""
   exit -1
fi

if test -z $prefix ;then
   echo ""
   echo "Please input the prefix (basename) of peaks files. -o [prefix]"
   echo ""
   exit -1
fi

if test -z $assemblyid ;then
   echo ""	
   echo "Please specify the reference genome assembly id of KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -g [assembly id]"
   echo ""
   exit -1
fi

# genome size. e.g. human(hs): 2.7e9; mouse(mm): 1.87e9; C.elegans(ce): 9e7; fruitfly(dm): 1.2e8; rat(rn): 2.5e9; zebrafish(danRer): 1e9.

if [[ $assemblyid == "hg18" ]] || [[ $assemblyid == "hg19" ]] || [[ $assemblyid == "hg38" ]] ;then
	genomesize="2.7e9"
elif [[ $assemblyid == "mm9" ]] || [[ $assemblyid == "mm10" ]] || [[ $assemblyid == "mm39" ]] ;then	
        genomesize="1.87e9"
elif [[ $assemblyid == "ce10" ]] || [[ $assemblyid == "ce11" ]] ;then	
	genomesize="9e7"
elif [[ $assemblyid == "dm3" ]] || [[ $assemblyid == "dm6" ]] ;then	
	genomesize="1.2e8"
elif [[ $assemblyid == "rn6" ]] || [[ $assemblyid == "rn7" ]] ;then 	
	genomesize="2.5e9"
elif [[ $assemblyid == "danRer10" ]] || [[ $assemblyid == "danRer11" ]] ;then	
	genomesize="1e9"
else 
echo ""
echo "Error: unsupported assembly id: $assemblyid."
echo ""
exit -1
fi


# setup the default options.
if test -z $peakscaller ;then
   peakscaller="macs2_and_epic2"
fi

if test -z $mode ;then
   mode="both"
fi

if [[ $peakscaller == "epic2" ]] ;then 
    mode="broad"	
elif [[ $mode == "broad" ]] ;then
    peakscaller="epic2"	
fi 

if [[ $peakscaller == "macs2" ]] ;then
    mode="sharp"
elif [[ $mode == "sharp" ]] ;then
    peakscaller="macs2"
fi

if [[ $peakscaller == "macs2_and_epic2" ]] ;then
   mode="both"
elif [[ $mode == "both" ]] ;then
   peakscaller="macs2_and_epic2"
fi

# test unsuported peak calling mode.
if test -n "$mode" && [[ $mode != "broad" ]] && [[ $mode != "sharp" ]] && [[ $mode != "both" ]] ;then
    echo ""
    echo "Error: unsupported peak calling mode: $mode. e.g. sharp, broad or both. DEFAULT: both."
    echo ""
    exit -1
fi

# test unsuported peak caller.
if [[ $peakscaller != "macs2" ]] && [[ $peakscaller != "epic2" ]] && [[ $peakscaller != "macs2_and_epic2" ]] ;then
   echo ""
   echo "Error: unsupported peak caller: $peakscaller. e.g. macs2 or epic2 or macs2_and_epic2. DEFAULT: macs2_and_epic2."
   echo " "
   exit -1
fi

# get the path of 'KAS-Analyzer KAS-seq' shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# setup the default pvalue or qvalue.
if test -z $cutoff ;then
   if [[ $peakscaller == "macs2" ]] ;then
      cutoff="0.01"
   elif [[ $peakscaller == "epic2" ]] ;then
      cutoff="0.05"
   elif [[ $peakscaller == "macs2_and_epic2" ]] ;then
      macs2_cutoff="0.01"
      epic2_cutoff="0.05"      
   fi	
fi

# call KAS-seq peaks without control using macs2
if test -z $control && [[ $peakscaller == "macs2" ]] && [[ $mode == "sharp" ]] ;then	

   echo $KASseq > .KASseq.txt
   sed -i "s/\,/ /g" .KASseq.txt
   KASseq_list=$(cat .KASseq.txt)
   echo "Call sharp (sp)KAS-seq peaks for $KASseq_list using macs2 ..." 
   echo ""
   macs2 callpeak -t $KASseq_list -n $prefix -g $genomesize -q $cutoff 2> ${prefix}_output_${peakscaller}.log
   rm -f .KASseq.txt
   echo "done."
   echo ""

elif test -z $control && [[ $peakscaller == "epic2" ]] && [[ $mode == "broad" ]] ;then
   echo $KASseq > .KASseq.txt
   sed -i "s/\,/ /g" .KASseq.txt
   KASseq_list=$(cat .KASseq.txt)
   echo "Call broad (sp)KAS-seq peaks for $KASseq_list using epic2 ..." 
   echo ""
   epic2 -t $KASseq_list --output $prefix --chromsizes ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes --fragment-size 300 --false-discovery-rate-cutoff $cutoff 2> ${prefix}_output_${peakscaller}.log
   rm -f .KASseq.txt
   echo "done."
   echo ""

elif test -z $control && [[ $peakscaller == "macs2_and_epic2" ]] && [[ $mode == "both" ]] ;then

   echo "Control(Input) data were not provided, call sharp and broad KAS-seq peaks using both macs2 and epic2 can't be performed"
   echo ""
   exit -1

# call KAS-seq peaks with control using macs2 and epic2.
elif test -n "$control" && [[ $peakscaller == "macs2" ]] && [[ $mode == "sharp" ]] ;then
   echo $KASseq > .KASseq.txt
   echo $control > .control.txt
   sed -i "s/\,/ /g" .KASseq.txt
   sed -i "s/\,/ /g" .control.txt
   KASseq_list=$(cat .KASseq.txt)
   control_list=$(cat .control.txt)
   echo "Call sharp (sp)KAS-seq peaks for $KASseq_list compaired to $control_list using macs2 ..." 
   echo ""
   macs2 callpeak -t $KASseq_list -c $control_list -n $prefix -g $genomesize -q $cutoff 2> ${prefix}_output_${peakscaller}.log

   grep ^chr ${prefix}_peaks.xls > ${prefix}_macs2_peaks.bed

   rm -f ${prefix}_model.r
   rm -f ${prefix}_summits.bed
   rm -f ${prefix}_peaks.narrowPeak
   rm -f ${prefix}_peaks.xls
   rm -f .KASseq.txt
   rm -f .control.txt
   echo "done."
   echo ""

elif test -n "$control" && [[ $peakscaller == "epic2" ]] && [[ $mode == "broad" ]] ;then	
   echo $KASseq > .KASseq.txt
   echo $control > .control.txt
   sed -i "s/\,/ /g" .KASseq.txt
   sed -i "s/\,/ /g" .control.txt
   KASseq_list=$(cat .KASseq.txt)
   control_list=$(cat .control.txt)
   echo "Call broad (sp)KAS-seq peaks for $KASseq_list compared to $control_list using epic2 ..."
   echo ""
   epic2 -t $KASseq_list -c $control_list --output $prefix --chromsizes ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes --fragment-size 300 --false-discovery-rate-cutoff $cutoff 2> ${prefix}_output_${peakscaller}.log

   sed "1d" $prefix | awk 'BEGIN {print "#Chromosome\tStart\tEnd\tPValue\tScore\tChIPCount\tInputCount\tFDR\tlog2FoldChange"} $10>=0.58 {printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$7,$8,$9,$10)}' > ${prefix}_epic2_peaks.bed
   
   rm -f $prefix
   rm -f .KASseq.txt
   rm -f .control.txt
   echo "done."
   echo ""

elif test -n "$control" && [[ $peakscaller == "macs2_and_epic2" ]] && [[ $mode == "both" ]] ;then
   echo $KASseq > .KASseq.txt
   echo $control > .control.txt
   sed -i "s/\,/ /g" .KASseq.txt
   sed -i "s/\,/ /g" .control.txt
   KASseq_list=$(cat .KASseq.txt)
   control_list=$(cat .control.txt)	
 
   echo "Call broad and sharp (sp)KAS-seq peaks for $KASseq_list compared to $control_list using epic2 and macs2 ..." 

   echo "Call broad (sp)KAS-seq peaks for $KASseq_list compared to $control_list using epic2 ..."
   echo ""  
   epic2 -t $KASseq_list -c $control_list --output $prefix --chromsizes ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes --fragment-size 300 --false-discovery-rate-cutoff $epic2_cutoff 2> ${prefix}_output_epic2.log

   sed "1d" $prefix | awk '$10>=0.58 {printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$7,$8,$9,$10,"peaks"NR)}' > ${prefix}_epic2_peaks.bed
   echo "done."
   echo ""

   echo "Call sharp (sp)KAS-seq peaks for $KASseq_list compared to $control_list using macs2 ..."
   echo ""
   macs2 callpeak -t $KASseq_list -c $control_list -n $prefix -g $genomesize -q $macs2_cutoff 2> ${prefix}_output_macs2.log

   grep ^chr ${prefix}_peaks.xls | sed "1d" > ${prefix}_macs2_peaks.bed
   echo "done."
   echo ""

   echo "Filter macs2 KAS-seq peaks with at least 5 fold change compared to its relative input control ..."
   echo ""
   awk '$8>=5 {print $0}' ${prefix}_macs2_peaks.bed | awk '{if ($4<1000 && $2-1000>=0) printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\n",$1,$2,$3,$1,$2-1000,$2,$1,$3,$3+1000); else if ($4>=1000 && $2-3000>=0) printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\n",$1,$2,$3,$1,$2-3000,$2,$1,$3,$3+3000)}' > ${prefix}_macs2_peaks.middle.left.right.3bed

   cat $KASseq_list | sortBed -i > ${prefix}_KAS-seq_reads.bed

   intersectBed -a ${prefix}_macs2_peaks.middle.left.right.3bed -b ${prefix}_KAS-seq_reads.bed -wa -F 0.5 -c | awk '{printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n",$4,$5,$6,$1,$2,$3,$7,$8,$9,$10)}' | intersectBed -a - -b ${prefix}_KAS-seq_reads.bed -wa -F 0.5 -c | awk '{printf("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n",$7,$8,$9,$4,$5,$6,$1,$2,$3,$10,$11)}' | intersectBed -a - -b ${prefix}_KAS-seq_reads.bed -wa -F 0.5 -c | awk '{if ($6-$5>=1000) printf("%s\t%d\t%d\t%.3f\t%.3f\t%.3f\n",$4,$5,$6,$10/($6-$5),$11/3000,$12/3000); else if ($6-$5<1000) printf("%s\t%d\t%d\t%.3f\t%.3f\t%.3f\n",$4,$5,$6,$10/($6-$5),$11/1000,$12/1000)}' | awk '$4>=$5+$6 {printf("%s\t%d\t%d\n",$1,$2,$3)}' > ${prefix}_macs2_peaks.filter.3bed
  
   rm -f $prefix 
   
   rm -f ${prefix}_model.r
   rm -f ${prefix}_summits.bed
   rm -f ${prefix}_peaks.narrowPeak
   rm -f ${prefix}_peaks.xls
   
   rm -f ${prefix}_macs2_peaks.middle.left.right.3bed
   rm -f ${prefix}_KAS-seq_reads.bed

   rm -f .KASseq.txt
   rm -f .control.txt
   
   echo "done."
   echo ""
   

   echo "Output sharp KAS-seq peaks ..."
   echo ""

   intersectBed -a ${prefix}_macs2_peaks.bed -b ${prefix}_macs2_peaks.filter.3bed -wa -f 1 | intersectBed -a - -b ${prefix}_epic2_peaks.bed -wa -f 0.5 | sortBed -i | awk 'BEGIN {print "#chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)"} {printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9)}' > ${prefix}_KAS-seq_peaks.sharp.bed

   rm -f ${prefix}_macs2_peaks.bed
   rm -f ${prefix}_macs2_peaks.filter.3bed
   rm -f ${prefix}_output_macs2.log 
   echo "done."
   echo ""

   echo "Output broad KAS-seq peaks ..."
   echo ""
   bedtools subtract -a ${prefix}_epic2_peaks.bed -b ${prefix}_KAS-seq_peaks.sharp.bed | awk 'BEGIN {print "#Chromosome\tStart\tEnd\tPValue\tScore\tStrand\tChIPCount\tInputCount\tFDR\tlog2FoldChange"} {printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10)}' > ${prefix}_KAS-seq_peaks.broad.bed

   rm -f ${prefix}_epic2_peaks.bed
   rm -f ${prefix}_output_epic2.log

   echo "done."
   echo ""

else
printHelpAndExit 0

fi

echo "'KAS-Analyzer peakscalling' run successfully!"
