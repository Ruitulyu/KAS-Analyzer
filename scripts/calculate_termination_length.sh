#!/bin/bash
# 'KAS-pipe2 termilength' was developed by Ruitu Lyu on 1-25-2022.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 termilength [ -h/--help ] [ -o prefix ] [ -t threads ] [ -b bin size ] [ -g assembly id ] [ -p peaks ] [ -l labels ] [ -k KAS-seq ]"
exampleHelp="Example: nohup KAS-pipe2 termilength -o KAS-seq_termination_length -t 10 -g mm10 -p peaks.txt -l labels.txt -k KAS-seq.txt &"
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 termilength' output files. REQUIET."
threadsHelp="-t [threads]: please specify the number of threads. DEFAULT: 1."
binsizeHelp="-b [bin size]: please specify the size of sliding bins used for transcription termination length calculation. DEFAULT:100."
peaksHelp="-p [peaks]: please input the text file containing the peaks file of (sp)KAS-seq data. REQUIET.
Example:
WT.peaks.rep1.bed
WT.peaks.rep2.bed              ---peaks.txt"
assemblyidHelp="-s [assembly id]: please specify the genome assembly id of (sp)KAS-seq data. -g [assembly id]. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
labelsHelp="-l [labels]: please input the text file containing the labels of (sp)KAS-seq data that show in 'KAS-pipe2 termilength' output files. DEFAULT: basename of KAS-seq file.
Example:
WT.rep1
WT.rep2                        ---labels.txt"
KASseqHelp="-k [KAS-seq]: please input the text file containing indexed (sp)KAS-seq bam files used for transcription termination length calculation. OPTIONAL.
Example:
KAS-seq.rep1.bam
KAS-seq.rep2.bam               ---KAS-seq.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 termilength' shell script is applied to calculate the transcription termination length of protein coding genes."

# print help.
printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$threadsHelp"
    echo -e ""
    echo -e "$binsizeHelp"
    echo -e ""
    echo -e "$peaksHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$labelsHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 KAS-seq' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

while getopts 'ho:t:p:b:g:l:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        o) prefix=$OPTARG ;;
        t) threads=$OPTARG ;;
        p) peaks=$OPTARG ;;
        b) binsize=$OPTARG ;;
        g) assemblyid=$OPTARG ;;
        l) labels=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $prefix ;then
   echo ""
   echo "Please provide the prefix (basename) of 'KAS-pipe2 termilength' output files. -o [prefix]"
   echo ""
   exit 1
fi

if test -z $assemblyid ;then
   echo ""
   echo "Please input the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra 
fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. -g [assembly id]"
   echo ""
   exit 1
fi

# supported assembly id.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
   echo ""
   echo "Error: unsupported assembly id: $assemblyid ; Please specify the reference genome assembly id of KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11."
   echo ""
   exit 1
fi

# setup the effective genome size.
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

if test -z $peaks ;then
   echo ""	
   echo "Please input the text file containing the peaks file of (sp)KAS-seq data. -p [peaks]."
   echo ""
   exit 1
fi 

if test -z $labels ;then
   echo ""
   echo "Please input the text file containing the labels of (sp)KAS-seq in 'KAS-pipe2 termination_length' output files. -l [labels]."
   echo ""
   exit 1
fi

# get the number of samples.
number_of_peaks=$( awk 'END {print NR}' $peaks )
number_of_labels=$( awk 'END {print NR}' $labels )

if [[ $number_of_labels != $number_of_peaks ]] ;then
   echo ""
   echo "Error: the number of labels isn't consistent with the number of peaks files."
   echo ""
   exit 1
fi

# setup the default options.
if test -z $threads ;then
   threads=1
fi

if test -z $binsize ;then
   binsize=100
fi

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)


if test -n "$KASseq" ;then
   
   number_of_samples=$( awk 'END {print NR}' $KASseq )
   if [[ $number_of_labels != $number_of_samples ]] ;then
      echo ""
      echo "Error: the number of labels isn't consistent with the number of indexed KAS-seq bam files."
      echo ""
      exit 1
   fi
	
   # normalize bam file into RPKM bigWig file.
   cat /dev/null > ${prefix}.KAS-seq.RPKM.bigWig.txt
   for ((i=1; i<=${number_of_samples}; i++))
   do
   sample_selected=$(sed -n ''$i'p' $KASseq)
   samtools index $sample_selected
   echo "Normalize $sample_selected using RPKM and output bigWig file ..."
   echo ""
   bamCoverage -b $sample_selected --outFileFormat bigwig -p $threads -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed --effectiveGenomeSize $genomesize --normalizeUsing RPKM -o ${sample_selected}.bigWig > /dev/null 2>&1
   echo ${sample_selected}.bigWig >> ${prefix}.KAS-seq.RPKM.bigWig.txt
   echo "done."
   echo ""
   done

   KASseq_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' ${prefix}.KAS-seq.RPKM.bigWig.txt)
   labels_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels)

   echo "Calculate the KAS-seq RPKM values on ${assemblyid} ${binsize}bp bins ..."
   echo ""
   distanceBetweenBins=$((-1*binsize/2))
   multiBigwigSummary bins -b $KASseq_list --labels $labels_list -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -bs $binsize -n $distanceBetweenBins -p $threads -out ${prefix}_on_${binsize}bp.bins.rmbl.npz --outRawCounts ${prefix}_on_${binsize}bp.bins.rmbl.tab

   sed "s/nan/0/g" ${prefix}_on_${binsize}bp.bins.rmbl.tab | sed "1d" > ${prefix}_on_${binsize}bp.bins.rmbl.bed
   echo "done."
   echo ""
   
   # filter the bins within KAS-seq peaks.
   for ((i=1; i<=${number_of_samples}; i++))
   do
   peak_selected=$( sed -n ''$i'p' $peaks )
   label_selected=$( sed -n ''$i'p' $labels )   
   echo "Calculate transcription termination length for $label_selected ..."
   echo ""
   awk -v x=$i '$(x+3)>=1 {printf("%s\t%d\t%d\t%d\n",$1,$2,$3,$(x+3))}' ${prefix}_on_${binsize}bp.bins.rmbl.bed | intersectBed -a - -b $peak_selected -wa -f 0.5 | sort -u | sortBed -i > ${prefix}_on_${label_selected}.bins.within_peaks.bed

   awk '{ if($6~"+"){printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2-1000,$3,$4,$5,$6)} else{printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$3+1000,$4,$5,$6)} }' ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.bed | awk '{ if($2<0) {printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,"0",$3,$4,$5,$6)} else{printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$3,$4,$5,$6)} }' > ${assemblyid}_Refseq.plus1kb.bed

   intersectBed -a ${prefix}_on_${label_selected}.bins.within_peaks.bed -b ${assemblyid}_Refseq.plus1kb.bed -v | mergeBed -i - -d 300 | awk '{printf("%s\t%d\t%d\t%d\n",$1,$2,$3,$3-$2)}' | intersectBed -a ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.terminator.bed -b - -wa -wb -f 0.3 | awk '{ if($6~"+"){printf("%s\t%d\t%d\t%s\t%d\n",$7,$2,$9,$4,$9-$2)} else{printf("%s\t%d\t%d\t%s\t%d\n",$7,$8,$3,$4,$3-$8)} }' | sort -k 10 -n -r | sort -u -k4,4 | sortBed -i > ${prefix}_${label_selected}.termination_length.without_header.txt 

   echo -e "chr\tstart\tend\tgene_name\ttermination_length" > ${prefix}_${label_selected}.header.txt
   cat ${prefix}_${label_selected}.header.txt ${prefix}_${label_selected}.termination_length.without_header.txt >  ${prefix}_${label_selected}.termination_length.txt
   
   rm -f ${prefix}_${label_selected}.termination_length.without_header.txt
   rm -f ${prefix}_${label_selected}.header.txt
   rm -f ${prefix}_on_${label_selected}.bins.within_peaks.bed
   rm -f ${assemblyid}_Refseq.plus1kb.bed
   echo "done."
   echo ""
   done
   
   rm -f ${prefix}_on_${binsize}bp.bins.rmbl.npz
   rm -f ${prefix}_on_${binsize}bp.bins.rmbl.tab
   rm -f ${prefix}_on_${binsize}bp.bins.rmbl.bed
   rm -f $KASseq_list
   rm -f ${prefix}.KAS-seq.RPKM.bigWig.txt

elif test -z $KASseq ;then
   
   distanceBetweenBins=$((binsize/2))
   echo "Generating ${binsize}bp bins with ${distanceBetweenBins}bp overlap for ${assemblyid} genome ..."
   echo ""
   bedtools makewindows -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes -w $binsize -s $distanceBetweenBins | intersectBed -a - -b ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -v > ${assemblyid}.bins.bed 
   echo "done."
   echo ""

   for ((i=1; i<=${number_of_peaks}; i++))
   do
   peak_selected=$( sed -n ''$i'p' $peaks )
   label_selected=$( sed -n ''$i'p' $labels )
   echo "Calculate transcription termination length for $label_selected ..."
   echo ""
   intersectBed -a ${assemblyid}.bins.bed -b $peak_selected -wa -f 0.5 | sort -u | sortBed -i > ${prefix}_on_${label_selected}.bins.within_peaks.bed

   awk '{ if($6~"+"){printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2-1000,$3,$4,$5,$6)} else{printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$3+1000,$4,$5,$6)} }' ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.bed | awk '{ if($2<0) {printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,"0",$3,$4,$5,$6)} else{printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$3,$4,$5,$6)} }' > ${assemblyid}_Refseq.plus1kb.bed

   intersectBed -a ${prefix}_on_${label_selected}.bins.within_peaks.bed -b ${assemblyid}_Refseq.plus1kb.bed -v | mergeBed -i - -d 300 | awk '{printf("%s\t%d\t%d\t%d\n",$1,$2,$3,$3-$2)}' | intersectBed -a ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.terminator.bed -b - -wa -wb -f 0.3 | awk '{ if($6~"+"){printf("%s\t%d\t%d\t%s\t%d\n",$7,$2,$9,$4,$9-$2)} else{printf("%s\t%d\t%d\t%s\t%d\n",$7,$8,$3,$4,$3-$8)} }' | sort -k 10 -n -r | sort -u -k4,4 | sortBed -i > ${prefix}_${label_selected}.termination_length.without_header.txt
   
   echo -e "chr\tstart\tend\tgene_name\ttermination_length" > ${prefix}_${label_selected}.header.txt
   cat ${prefix}_${label_selected}.header.txt ${prefix}_${label_selected}.termination_length.without_header.txt > ${prefix}_${label_selected}.termination_length.txt
   
   rm -f ${prefix}_${label_selected}.header.txt
   rm -f ${prefix}_${label_selected}.termination_length.without_header.txt
   rm -f ${prefix}_on_${label_selected}.bins.within_peaks.bed
   rm -f ${assemblyid}_Refseq.plus1kb.bed
   echo "done."
   echo ""
   done
   
   rm -f ${assemblyid}.bins.bed

fi 

echo "'KAS-pipe2 termilength' run successfully!"
