#!/bin/bash
# 'KAS-Analyzer genomicdist' was developed by Ruitu Lyu on 12-12-2021.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-Analyzer genomicdist [ -h/--help ] [ -o prefix ] [ -c ] [ -p peaks ] [ -s assembly id ]"
exampleHelp="Example: nohup KAS-Analyzer genomicdist -o KAS-seq_genomic_distribution -p KAS-seq_peaks.bed -s hg19 &"
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-Analyzer genomicdist' output files. DEFAULT: basename of KAS-seq peak file."
controlHelp="-c: please specify if the percentages of normal genomic feature distribution is generated, which is regard as a control. DEFAULT: off."
peaksHelp="-p [peaks]: please input the KAS-seq peak or differential KAS-seq peak file. REQUIRED."
assemblyid="-s [assembly id]: please specify the reference genome assembly id, e.g. Human: hg18, hg19, hg38, hs1; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. REQUIRED."
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-Analyzer genomicdist' shell script is applied to calculate and plot the percentages of (sp)KAS-seq peaks distribution on genomic features (Promoter(TSS +/-1kb), Exon, Intron, Terminal(TES+3kb) and Intergenic regions)."

printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$controlHelp"
    echo -e ""
    echo -e "$peaksHelp"
    echo -e ""
    echo -e "$assemblyid"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer genomicdist' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'ho:cp:s:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        o) prefix=$OPTARG ;;
        c) control="on" ;;
        p) peaks=$OPTARG ;;
        s) assemblyid=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $assemblyid ;then
   echo ""
   echo "Please specify the reference genome assembly id of your KAS-seq data. e.g. Human: hg18, hg19, hg38, hs1; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -s [assembly id]"
   echo ""
   exit -1
fi

if test -z $peaks ;then
   echo ""
   echo "Please input the KAS-seq peak or any peaks file, which will be used to calculate its distribution on genomic features. -p [peaks]"
   echo ""
   exit -1
fi

# setup the default parameters.
if test -z $control ;then
   control="off"
fi

if test -z $prefix ;then 
   prefix=$(basename ${peaks} .bed)
fi

# get the absolute path of 'KAS-Analyzer genomicdist' shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# test if $assemblyid is supported.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "hs1" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
   echo ""
   echo "Error: unsupported assembly id : ${assemblyid} in 'KAS-Analyzer genomicdist'. Supported assembly id: Human: hg18, hg19, hg38, hs1; Mouse: mm9,mm10,mm39; Fruitfly: dm3, dm6; Rat: rn6, rn7; C.elegans: ce10, ce11; Zebra fish: danRer10, danRer11. -s [assembly id]"
   echo ""
   exit -1
fi

# control distribution.
if [[ $control == "on" ]] ;then

   # get the number of 50bp bins on different genomic features. 
   echo "Calculate control 'Promoter' percentage ..."
   echo ""
   control_promoter=$( wc -l ${SH_SCRIPT_DIR}/../genomic_bins/${assemblyid}/${assemblyid}.50bp.promoter.bed | awk '{print $1}' )
   echo "done."
   echo ""
   
   echo "Calculate control 'Exon' percentage ..."
   echo ""
   control_exon=$( wc -l ${SH_SCRIPT_DIR}/../genomic_bins/${assemblyid}/${assemblyid}.50bp.exon.bed | awk '{print $1}' )
   echo "done."
   echo ""

   echo "Calculate control 'Intron' percentage ..."
   echo ""
   control_intron=$( wc -l ${SH_SCRIPT_DIR}/../genomic_bins/${assemblyid}/${assemblyid}.50bp.intron.bed | awk '{print $1}' )
   echo "done."
   echo ""

   echo "Calculate control 'Terminal3kb' percentage ..."
   echo ""
   control_terminal3kb=$( wc -l ${SH_SCRIPT_DIR}/../genomic_bins/${assemblyid}/${assemblyid}.50bp.terminal3kb.bed | awk '{print $1}' )
   echo "done."
   echo ""

   echo "Calculate control 'Intergenic' percentage ..."
   echo ""
   control_intergenic=$( wc -l ${SH_SCRIPT_DIR}/../genomic_bins/${assemblyid}/${assemblyid}.50bp.intergenic.bed | awk '{print $1}' )
   control_sum=$(($control_exon + $control_intergenic + $control_intron + $control_promoter + $control_terminal3kb))
   echo "done."
   echo ""
   
   # calculate the percentage.
   control_promoter_percentage=$( awk -v x=$control_promoter -v y=$control_sum 'BEGIN{printf "%.2f\n",x*100/y}' )
   control_exon_percentage=$( awk -v x=$control_exon -v y=$control_sum 'BEGIN{printf "%.2f\n",x*100/y}' )
   control_intron_percentage=$( awk -v x=$control_intron -v y=$control_sum 'BEGIN{printf "%.2f\n",x*100/y}' )
   control_terminal3kb_percentage=$( awk -v x=$control_terminal3kb -v y=$control_sum 'BEGIN{printf "%.2f\n",x*100/y}' )
   control_intergenic_percentage=$( awk -v x=$control_intergenic -v y=$control_sum 'BEGIN{printf "%.2f\n",x*100/y}' )

   # ouput to ${prefix}_control_distribution.txt.
   echo -e "Features\tPercentage" > ${assemblyid}_genomic_distribution.txt
   echo -e "Promoter\t${control_promoter_percentage}" >> ${assemblyid}_genomic_distribution.txt
   echo -e "Exon\t${control_exon_percentage}" >> ${assemblyid}_genomic_distribution.txt
   echo -e "Intron\t${control_intron_percentage}" >> ${assemblyid}_genomic_distribution.txt
   echo -e "Terminal3kb\t${control_terminal3kb_percentage}" >> ${assemblyid}_genomic_distribution.txt
   echo -e "Intergenic\t${control_intergenic_percentage}" >> ${assemblyid}_genomic_distribution.txt
   echo "$assemblyid genomic distribution percentage calculation done!"
   echo ""
fi

# distribution of KAS-seq peaks or any peaks.
# get the number of peaks 50bp bins on different genomic features.
echo "Calculate $peaks 'Promoter' percentage ..."
echo ""
peaks_promoter=$( intersectBed -a ${SH_SCRIPT_DIR}/../genomic_bins/${assemblyid}/${assemblyid}.50bp.promoter.bed -b $peaks -wa | sort -u |wc -l )
echo "done."
echo ""

echo "Calculate $peaks 'Exon' percentage ..."
echo ""
peaks_exon=$( intersectBed -a ${SH_SCRIPT_DIR}/../genomic_bins/${assemblyid}/${assemblyid}.50bp.exon.bed -b $peaks -wa | sort -u |wc -l )
echo "done."
echo ""

echo "Calculate $peaks 'Intron' percentage ..."
echo ""
peaks_intron=$( intersectBed -a ${SH_SCRIPT_DIR}/../genomic_bins/${assemblyid}/${assemblyid}.50bp.intron.bed -b $peaks -wa | sort -u |wc -l )
echo "done."
echo ""

echo "Calculate $peaks 'Intergenic' percentage ..."
echo ""
peaks_intergenic=$( intersectBed -a ${SH_SCRIPT_DIR}/../genomic_bins/${assemblyid}/${assemblyid}.50bp.intergenic.bed -b $peaks -wa | sort -u |wc -l )
echo "done."
echo ""

echo "Calculate $peaks 'Terminal3kb' percentage ..."
echo ""
peaks_terminal3kb=$( intersectBed -a ${SH_SCRIPT_DIR}/../genomic_bins/${assemblyid}/${assemblyid}.50bp.terminal3kb.bed -b $peaks -wa | sort -u |wc -l )
peaks_sum=$((peaks_exon + peaks_intergenic + peaks_intron + peaks_promoter + peaks_terminal3kb))
echo "done."
echo ""

# calculate the percentage of peaks 50bp bins on different genomic features.
peaks_promoter_percentage=$( awk -v x=$peaks_promoter -v y=$peaks_sum 'BEGIN{printf "%.2f\n",x*100/y}' )
peaks_exon_percentage=$( awk -v x=$peaks_exon -v y=$peaks_sum 'BEGIN{printf "%.2f\n",x*100/y}' )
peaks_intron_percentage=$( awk -v x=$peaks_intron -v y=$peaks_sum 'BEGIN{printf "%.2f\n",x*100/y}' )
peaks_terminal3kb_percentage=$( awk -v x=$peaks_terminal3kb -v y=$peaks_sum 'BEGIN{printf "%.2f\n",x*100/y}' )
peaks_intergenic_percentage=$( awk -v x=$peaks_intergenic -v y=$peaks_sum 'BEGIN{printf "%.2f\n",x*100/y}' )


echo -e "Features\tPercentage" > ${prefix}_peaks_${assemblyid}_genomic_distribution.txt
echo -e "Promoter\t${peaks_promoter_percentage}" >> ${prefix}_peaks_${assemblyid}_genomic_distribution.txt
echo -e "Exon\t${peaks_exon_percentage}" >> ${prefix}_peaks_${assemblyid}_genomic_distribution.txt
echo -e "Intron\t${peaks_intron_percentage}" >> ${prefix}_peaks_${assemblyid}_genomic_distribution.txt
echo -e "Terminal3kb\t${peaks_terminal3kb_percentage}" >> ${prefix}_peaks_${assemblyid}_genomic_distribution.txt
echo -e "Intergenic\t${peaks_intergenic_percentage}" >> ${prefix}_peaks_${assemblyid}_genomic_distribution.txt
echo "Genomic distribution percentage of $peaks calculation done."
echo ""

echo "Plot the pie chart of genomic distribution of KAS-seq peaks."
Rscript --vanilla ${SH_SCRIPT_DIR}/../R/Piechart_plot.R ${prefix}_peaks_${assemblyid}_genomic_distribution.txt
echo "done."
echo ""

mv KAS-seq_peaks_genomic_distribution_pie_chart.png ${prefix}_KAS-seq_genomic_distribution_pie_chart.png
mv KAS-seq_peaks_genomic_distribution_pie_chart.svg ${prefix}_KAS-seq_genomic_distribution_pie_chart.svg
# rm -f ${prefix}_peaks_${assemblyid}_genomic_distribution.txt

echo "'KAS-Analyzer genomicdist' run successfully!"
