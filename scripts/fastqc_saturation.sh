#!/bin/bash
# 'KAS-pipe2 saturation' was developed by Ruitu Lyu on 12-16-2021.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 saturation [ -h/--help ] [ -o prefix ] [ -s assembly id ] [ -c control ] [ -k KAS-seq ]"
exampleHelp="Example: nohup KAS-pipe2 saturation -o KAS-seq_saturation -c KAS-seq_Input.bed -k KAS-seq.bed &"
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 saturation' output files. Default: basename of KAS-seq data."
assemblyidHelp="-s [assembly id]: please specify the genome assembly id of (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
controlHelp="-c [control]: please input the control data (input of (sp)KAS-seq data) containing uniquely mapped reads. e.g. -c KAS-seq_Input.bed. REQUIRED. Note: reads number of KAS-seq and Input should be similar."
KASseqHelp="-k [KAS-seq]: please input the KAS-seq data containing uniquely mapped reads. e.g. -k KAS-seq.bed. REQUIRED."
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 saturation' shell script is applied to evaluate the saturation of (sp)KAS-seq data."

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
    echo -e "$controlHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

if ! type macs2 > /dev/null 2>&1 ;then
   echo "macs2 was not installed or not export to the \$PATH'"
   echo ""
   echo "Install macs2 with 'conda install -c bioconda macs2' or refer the official website of 'macs2'."
   echo ""
   exit 1
fi

# if no parameters was provided, 'KAS-pipe2 saturation' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'ho:s:c:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        o) prefix=$OPTARG ;;
	s) assemblyid=$OPTARG ;;
        c) control=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $assemblyid ;then
   echo ""      
   echo "Please specify the reference genome assembly id of (sp)KAS-seq data. -s [assembly id]. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11."
   echo ""
   exit -1
fi

# if test -z $control ;then
#   echo ""
#   echo "Please input the (sp)KAS-seq control data containing uniquely mapped reads. e.g. -k KAS-seq_Input.bed. REQUIRED. -c [control]"
#   echo ""
#   printHelpAndExit 0
# fi

if test -z $KASseq ;then
   echo ""
   echo "Please input the (sp)KAS-seq data containing uniquely mapped reads. e.g. -k KAS-seq.bed. REQUIRED. -k [KAS-seq]"
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
   echo "Error: unsupported assembly id: $assemblyid. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11."
   echo ""
   exit -1
fi

# setup the default of $prefix.
if test -z $prefix ;then
   prefix=$( basename ${KASseq} .bed )
fi

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# permform saturation analysis with KAS-seq Input.
if test -n "$control" ;then

   echo "Call (sp)KAS-seq peaks for $KASseq vs $control using 100% unique mapped reads with macs2 ..."
   echo ""
   sed -i '/^chrM/d' $KASseq
   sed -i '/^chrM/d' $control
   macs2 callpeak -t $KASseq -c $control -n ${prefix} --broad --nomodel --verbose 0 -g $genomesize --broad-cutoff 0.01 -q 0.01
   echo "done."
   echo ""

   echo "Calculate the (sp)KAS-seq peaks coverage."
   echo ""
   peaks_coverage=$( awk 'BEGIN {sum=0} {sum+=$3-$2} END {print sum}' ${prefix}_peaks.broadPeak )
   echo "done."
   echo ""

   echo -e "Percentage\tPercent" > ${prefix}_saturation.txt

   rm -f ${prefix}_peaks.xls
   rm -f ${prefix}_model.r
   rm -f ${prefix}_peaks.gappedPeak
   rm -f ${prefix}_peaks.broadPeak
   echo "clean up the macs2 peaks calling intermediate files. done."
   echo ""

   echo "Generate random (sp)KAS-seq and Input uniquely mapped reads ..."
   cat $KASseq | shuf > ${prefix}_KAS-seq.shuf.bed
   cat $control | shuf > ${prefix}_control.shuf.bed
   echo "done."
   echo ""


   KASseq_num=$( wc -l ${prefix}_KAS-seq.shuf.bed | awk '{print $1}' )
   control_num=$( wc -l ${prefix}_control.shuf.bed | awk '{print $1}' )

   for ((i=1; i<=20; i++))
   do
   percentage_reads=$((i*5)) 	
   echo "Generate $percentage_reads% of $KASseq and $control uniquely mapped reads."
   echo ""
   head -n $[KASseq_num*i*5/100] ${prefix}_KAS-seq.shuf.bed | sortBed -i > ${prefix}_KAS-seq.subset.bed
   head -n $[control_num*i*5/100] ${prefix}_control.shuf.bed | sortBed -i > ${prefix}_control.subset.bed 
   echo "done."
   echo ""

   echo "Call (sp)KAS-seq peaks for $KASseq vs $control using $percentage_reads% of uniquely mapped reads with macs2 ..."
   echo ""
   macs2 callpeak -t ${prefix}_KAS-seq.subset.bed -c ${prefix}_control.subset.bed -n ${prefix}_subset --verbose 0 --broad --nomodel -g $genomesize --broad-cutoff 0.01 -q 0.01 
   echo "done."
   echo ""

   echo "Calculate the (sp)KAS-seq peaks coverage of $percentage_reads% of $KASseq and $control uniquely mapped reads and recall percentage."
   echo ""
   peaks_coverage_subset=$( awk 'BEGIN {sum=0} {sum+=$3-$2} END {print sum}' ${prefix}_subset_peaks.broadPeak )
   percentage=$( awk -v x=$peaks_coverage_subset -v y=$peaks_coverage 'BEGIN{printf "%d\n",x*100/y}' )
   echo "done."
   echo ""

   if [ $percentage -ge 100 ]; then
       let percentage=100
   fi

   echo -e "$((i*5))\t$percentage" >> ${prefix}_saturation.txt
   rm -f ${prefix}_subset_peaks.xls
   rm -f ${prefix}_subset_model.r
   rm -f ${prefix}_subset_peaks.gappedPeak
   rm -f ${prefix}_subset_peaks.broadPeak
   rm -f ${prefix}_KAS-seq.subset.bed
   rm -f ${prefix}_control.subset.bed
   echo "clean up the macs2 peaks calling intermediate files with $percentage_reads% of uniquely mapped reads. done."
   echo ""

   done

   rm -f ${prefix}_KAS-seq.shuf.bed
   rm -f ${prefix}_control.shuf.bed

elif test -z $control ;then
   # remove the reads on mitochondrial chromosome.
   sed -i '/^chrM/d' $KASseq

   echo "Call (sp)KAS-seq peaks for $KASseq using 100% unique mapped reads with macs2 ..."
   echo ""
   macs2 callpeak -t $KASseq -n ${prefix} --broad -g $genomesize --verbose 0 --nomodel --broad-cutoff 0.01 -q 0.01
   echo "done."
   echo ""

   echo "Calculate the (sp)KAS-seq peaks coverage."
   echo ""
   peaks_coverage=$( awk 'BEGIN {sum=0} {sum+=$3-$2} END {print sum}' ${prefix}_peaks.broadPeak )
   echo "done."
   echo ""
   
   # output the header in the ${prefix}_saturation.txt.
   echo -e "Percentage\tPercent" > ${prefix}_saturation.txt

   rm -f ${prefix}_peaks.xls
   rm -f ${prefix}_model.r
   rm -f ${prefix}_peaks.gappedPeak
   rm -f ${prefix}_peaks.broadPeak
   echo "Clean up the macs2 peaks calling intermediate files. done."
   echo ""

   echo "Generate random (sp)KAS-seq uniquely mapped reads ..."
   echo ""
   cat $KASseq | shuf > ${prefix}_KAS-seq.shuf.bed
   echo "done."
   echo ""

   KASseq_num=$( wc -l $KASseq | awk '{print $1}' )

   for ((i=1; i<=20; i++))
   do
   percentage_reads=$((i*5))
   echo "Generate $percentage_reads% of $KASseq uniquely mapped reads."	
   echo ""
   head -n $[KASseq_num*i*5/100] ${prefix}_KAS-seq.shuf.bed | sortBed -i > ${prefix}_KAS-seq.subset.bed
   echo "done."
   echo ""

   echo "Call (sp)KAS-seq peaks for $KASseq without input using $percentage_reads% of uniquely mapped reads with macs2 ..."
   echo ""
   macs2 callpeak -t ${prefix}_KAS-seq.subset.bed -n ${prefix}_subset --verbose 0 --broad --nomodel -g $genomesize --broad-cutoff 0.01 -q 0.01
   echo "done."
   echo ""

   echo "Calculate the (sp)KAS-seq peaks coverage of $percentage_reads% of $KASseq uniquely mapped reads and recall percentage."
   echo ""
   peaks_coverage_subset=$(awk 'BEGIN {sum=0} {sum+=$3-$2} END {print sum}' ${prefix}_subset_peaks.broadPeak)
   percentage=$( awk -v x=$peaks_coverage_subset -v y=$peaks_coverage 'BEGIN{printf "%d\n",x*100/y}' )
   echo "done."
   echo ""

   if [[ $percentage -ge 100 ]]; then
      let percentage=100
   fi
   
   # output the header in the ${prefix}_saturation.txt.
   echo -e "$((i*5))\t$percentage" >> ${prefix}_saturation.txt
   rm -f ${prefix}_subset_peaks.xls
   rm -f ${prefix}_subset_model.r
   rm -f ${prefix}_subset_peaks.gappedPeak
   rm -f ${prefix}_subset_peaks.broadPeak
   rm -f ${prefix}_KAS-seq.subset.bed

   echo "clean up the macs2 peaks calling intermediate files with $percentage_reads% of uniquely mapped reads without input data. done."
   done

   rm -f ${prefix}_KAS-seq.shuf.bed
fi 

echo "Plot the saturation analysis ..."
echo ""
Rscript --vanilla ${SH_SCRIPT_DIR}/../R/Saturation_lineplot.R ${prefix}_saturation.txt
echo "done."
echo ""

mv KAS-seq_saturation_plot.png ${prefix}_KAS-seq_saturation_plot.png
mv KAS-seq_saturation_plot.svg ${prefix}_KAS-seq_saturation_plot.svg
# rm -f ${prefix}_saturation.txt

echo "'KAS-pipe2 saturation' run successfully!"
