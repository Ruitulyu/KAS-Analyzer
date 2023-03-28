#!/bin/bash
# 'KAS-Analyzer PCA' was developed by Ruitu Lyu on 12-16-2021.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-Analyzer PCA [ -h/--help ] [ -o prefix ] [ -t threads ] [ -r regions ] [ -s assembly id ] [ -b bin size ] [ -p peaks] [ -l labels ] [ -k KAS-seq ]"
exampleHelp="Example: nohup KAS-Analyzer PCA -o KAS-seq_PCA -t 10 -r bin -s mm10 -l labels.txt -k KAS-seq.txt &"
prefixHelp="-o [KAS-seq_PCA]: please input the prefix (basename) of 'KAS-Analyzer PCA' output files. REQUIRED."
threadsHelp="-t [threads]: please specify the number of threads used for perform PCA analysis. DEFAULT: 1."
regionsHelp="-r [regions]: please specify the regions used to perform PCA analysis. e.g. promoter, genebody, peak or bin. DEFAULT: bin."
assemblyidHelp="-s [assemblyid]: please specify the reference genome assembly id of your (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
binsizeHelp="-b [bin size]: please specify the bin size for bins mode: '-r bins'. DEFAULT: 10000."
peaksHelp="-p [peaks file]: please input the custom regions file that used to perform PCA analysis. REQUIRED in 'peak' mode."
labelsHelp="-l [labels.txt]: please input the text file containing the labels of (sp)KAS-seq data that shown in the PCA plot. DEFAULT: basename of KAS-seq data.
Example:
0h
4h
8h
12h
24h                            ---labels.txt"
KASseqHelp="-k [KAS-seq.txt]: please input the text file containing indexed (sp)KAS-seq bam files. The order and number of (sp)KAS-seq data should be the consistent with the labels file if specified. REQUIRED.
Example:
KAS-seq.0h.bam
KAS-seq.4h.bam
KAS-seq.8h.bam
KAS-seq.12h.bam
KAS-seq.24h.bam                --KAS-seq.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-Analyzer PCA' shell script is applied to perform Principal Component Analysis (PCA) for (sp)KAS-seq data."

# print help function.
printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$threadsHelp"
    echo -e ""
    echo -e "$regionsHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$regionsHelp"
    echo -e ""
    echo -e "$binsizeHelp"
    echo -e ""
    echo -e "$peaksHelp"
    echo -e ""
    echo -e "$labelsHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer PCA' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'ho:t:r:s:b:p:l:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        o) prefix=$OPTARG ;;
	t) threads=$OPTARG ;;
	r) regions=$OPTARG ;;
	s) assemblyid=$OPTARG ;;
	b) binsize=$OPTARG ;;
        p) peaks=$OPTARG ;;
        l) labels=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# test if deeptools was installed.
if ! type deeptools > /dev/null 2>&1 ;then
   echo "deeptools was not installed or not export to the \$PATH'"
   echo ""
   echo "Install deeptools with 'conda install -c bioconda deeptools' or refer the official website of 'deeptools'."
   echo ""
   exit 1
fi

# Required options.
if test -z $prefix ;then
   echo ""
   echo "Please input the prefix (basename) of 'KAS-Analyzer PCA' output files. -o [prefix]"
   echo ""
   exit -1
fi

if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing indexed (sp)KAS-seq bam files. -k [KAS-seq]"
   echo ""
   exit -1
fi

if test -z $assemblyid ;then
   echo ""
   echo "Please specify the reference genome assembly id of KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -s [assembly id]"
   echo ""
   exit -1
fi

# supported assembly id.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
    echo ""
    echo "Error: unsupported assembly id: $assemblyid ; Please specify the reference genome assembly id of (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11."
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

# setup the default parameters.
if test -z $threads ;then
    threads=1
fi

if test -z $regions ;then
    regions="bin"
fi
# setup the default binsize. 
if test -n "$regions" ;then
   if [[ $regions == "peak" ]] && test -z $peaks ;then
      echo ""
      echo "Please input the custom regions file that used to perform PCA analysis. -p [peaks]"
      echo ""
      exit -1

   elif [[ $regions == "bin" ]] && test -z $binsize ;then
      binsize=10000
   fi
    
elif [[ $regions != "peak" ]] && [[ $regions != "bin" ]] && [[ $regions != "promoter" ]] && [[ $regions != "genebody" ]] ;then
   echo ""	    
   echo "Error: unsupported region types: $regions."
   echo ""
   exit 1
fi   

# get the number of KAS-seq data.
number_of_samples=$(awk 'END {print NR}' $KASseq)

# set up the label, if labels.txt is not provided.
if test -z $labels ;then
   for ((i=1; i<=${number_of_samples}; i++))
   do
   sample_selected=$(sed -n ''$i'p' $KASseq)
   label_basename=$(basename ${sample_selected} .bam)
   echo $label_basename >> ${prefix}.labels_basename.txt
   done
   labels = "${prefix}.labels_basename.txt"     

else
    number_of_labels=$( awk 'END {print NR}' $labels )
    if [[ $number_of_labels != $number_of_samples ]] ;then
       echo ""
       echo "Error:the number of labels isn't consistent with the number of (sp)KAS-seq data!" 
       echo ""
       exit 1
    fi	  
fi 
    	
# get the path of shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
      
# normalize bam file into RPKM bigWig file.
cat /dev/null > ${prefix}.KAS-seq.RPKM.bigWig.txt
for ((i=1; i<=${number_of_samples}; i++))
do
sample_selected=$(sed -n ''$i'p' $KASseq)
samtools index $sample_selected
echo "Normalize $sample_selected using RPKM and output bigWig file ..."
echo ""
bamCoverage -b $sample_selected --outFileFormat bigwig -p $threads -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed --effectiveGenomeSize $genomesize --normalizeUsing RPKM -o ${sample_selected}.bigWig  > /dev/null 2>&1
echo ${sample_selected}.bigWig >> ${prefix}.KAS-seq.RPKM.bigWig.txt
echo "done."
echo ""
done

# get the sample and labels list.
KASseq_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' ${prefix}.KAS-seq.RPKM.bigWig.txt)
labels_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $labels)

if [[ $regions == "bin" ]] ;then
   echo "Generating KAS-seq density matrix on ${assemblyid} ${binsize}bp bins"
   echo ""
   multiBigwigSummary bins -b $KASseq_list --labels $labels_list --binSize $binsize -p $threads --blackListFileName ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -out ${prefix}_on_${binsize}_bins.npz --outRawCounts ${prefix}_on_${binsize}_bins.tab > /dev/null 2>&1

   sed "s/nan/0/g" ${prefix}_on_${binsize}_bins.tab | sed "1d" > ${prefix}_on_${binsize}_bins.bed
   echo "done."
   echo ""

   # calculate the average value of every single row in the table.
   echo "Calculating KAS-seq density average for every bin in matrix ..."
   echo ""
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/(NF-3) }' ${prefix}_on_${binsize}_bins.bed > ${prefix}_on_${binsize}_bins.average
   echo "done."
   echo ""

   # filter the bins with averaged KAS-seq expressiong lower than 2.
   echo "Generating the final KAS-seq density matrix on ${assemblyid} $binsize bins ..."
   paste ${prefix}_on_${binsize}_bins.average ${prefix}_on_${binsize}_bins.bed | awk '$1>=2 {print $0}' > ${prefix}_on_${binsize}_bins.filter.bed
   cut -f1,2,3,4 --complement ${prefix}_on_${binsize}_bins.filter.bed | awk '{for(i=1;i<=NF;i++){printf "%.2f\t", $i}; printf "\n"}' > ${prefix}_on_${binsize}_bins.filter.matrix
   awk '{printf("%s\n","bin-"$2"-"$3"-"$4)}' ${prefix}_on_${binsize}_bins.filter.bed > ${prefix}_on_${binsize}_bins.filter.rowname

   paste ${prefix}_on_${binsize}_bins.filter.rowname ${prefix}_on_${binsize}_bins.filter.matrix > ${prefix}_on_${binsize}_bins.filter.without_header.txt

   echo -e "" > ${prefix}.header1.txt
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels > ${prefix}.header2.txt
   paste ${prefix}.header1.txt ${prefix}.header2.txt > ${prefix}.header.txt

   cat ${prefix}.header.txt ${prefix}_on_${binsize}_bins.filter.without_header.txt > ${prefix}_on_${assemblyid}_${regions}.txt
   echo "done."
   echo ""

   # clean up the intermediate files.
   echo "Clean up the intermediate files."
   rm -f ${prefix}_on_${binsize}_bins.tab
   rm -f ${prefix}_on_${binsize}_bins.npz
   rm -f ${prefix}_on_${binsize}_bins.bed
   rm -f ${prefix}_on_${binsize}_bins.average
   rm -f ${prefix}_on_${binsize}_bins.filter.bed
   rm -f ${prefix}_on_${binsize}_bins.filter.matrix
   rm -f ${prefix}_on_${binsize}_bins.filter.rowname
   rm -f ${prefix}_on_${binsize}_bins.filter.without_header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}.header.txt
   echo "done."
   echo ""

elif [[ $regions == "peak" ]] ;then
   # get the basename of peak file.
   peaks_basename=$(basename ${peaks} .bed)
   
   echo -e "Generating KAS-seq density matrix on ${peaks} ..."
   echo ""
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED $peaks --labels $labels_list -p $threads -out ${prefix}_on_${peaks_basename}.npz --outRawCounts ${prefix}_on_${peaks_basename}.tab > /dev/null 2>&1
   echo "done."
   echo ""

   echo "Generating the final KAS-seq density matrix on $peaks ..."
   sed "s/nan/0/g" ${prefix}_on_${peaks_basename}.tab | sed "1d" > ${prefix}_on_${peaks_basename}.bed
   cut -f1,2,3 --complement ${prefix}_on_${peaks_basename}.bed | awk '{for(i=1;i<=NF;i++){printf "%.2f\t", $i}; printf "\n"}' > ${prefix}_on_${peaks_basename}.matrix
   awk '{printf("%s\n","peak-"$1"-"$2"-"$3)}' ${prefix}_on_${peaks_basename}.bed > ${prefix}_on_${peaks_basename}.rowname

   paste ${prefix}_on_${peaks_basename}.rowname ${prefix}_on_${peaks_basename}.matrix > ${prefix}_on_${peaks_basename}.without_header.txt
   
   echo -e "" > ${prefix}.header1.txt
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels > ${prefix}.header2.txt
   paste ${prefix}.header1.txt ${prefix}.header2.txt > ${prefix}.header.txt
   cat ${prefix}.header.txt ${prefix}_on_${peaks_basename}.without_header.txt > ${prefix}_on_${assemblyid}_${regions}.txt
   
   echo "done."
   echo ""
   
   echo "Clean up the intermediate files."
   rm -f ${prefix}_on_${peaks_basename}.npz
   rm -f ${prefix}_on_${peaks_basename}.tab
   rm -f ${prefix}_on_${peaks_basename}.bed
   rm -f ${prefix}_on_${peaks_basename}.matrix
   rm -f ${prefix}_on_${peaks_basename}.rowname
   rm -f ${prefix}_on_${peaks_basename}.without_header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}.header.txt
   echo "done."
   echo ""

elif [[ $regions == "promoter" ]] || [[ $regions == "genebody" ]] ;then

   echo "Calculating KAS-seq density matrix on ${assemblyid} Refseq ${regions} ..."
   echo ""
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.${regions}.bed --labels $labels_list -p $threads -out ${prefix}_on_${assemblyid}_Refseq.${regions}.npz --outRawCounts ${prefix}_on_${assemblyid}_Refseq.${regions}.tab > /dev/null 2>&1

   sed "s/nan/0/g" ${prefix}_on_${assemblyid}_Refseq.${regions}.tab | sed "1d" | sortBed -i > ${prefix}_on_${assemblyid}_Refseq.${regions}.bed
   sortBed -i ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.${regions}.bed | awk '{printf("%s\n",$4"-"$1"-"$2"-"$3)}' > ${assemblyid}_Refseq.${regions}.genenames

   echo "done."
   echo ""

   # calculate the average KAS expression of every single row in the table.
   echo "Calculating the average KAS-seq density for every ${regions} in ${assemblyid}_Refseq.${regions}.bed ..."
   echo ""
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/(NF-3) }' ${prefix}_on_${assemblyid}_Refseq.${regions}.bed > ${prefix}_on_${assemblyid}_Refseq.${regions}.average
   echo "done."
   echo ""

   # filter the ${regions} with KAS expression over 2.
   echo "Filtering matrix of KAS-seq density on ${assemblyid}_Refseq.${regions}.bed ..."
   echo ""
   paste ${prefix}_on_${assemblyid}_Refseq.${regions}.average ${assemblyid}_Refseq.${regions}.genenames ${prefix}_on_${assemblyid}_Refseq.${regions}.bed | awk '$1>=2 {print $0}' > ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed
   echo "done."
   echo ""	

   echo "Generating the final KAS-seq density matrix on ${assemblyid}_Refseq.${regions}.bed ..."
   echo ""
   cut -f1,2,3,4,5 --complement ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed | awk '{for(i=1;i<=NF;i++){printf "%.2f\t", $i}; printf "\n"}' > ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.matrix
   awk '{printf("%s\n",$2)}' ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed > ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.genenames

   paste ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.genenames ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.matrix > ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt

   echo -e "" > ${prefix}.header1.txt
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels > ${prefix}.header2.txt
   paste ${prefix}.header1.txt ${prefix}.header2.txt > ${prefix}.header.txt

   cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt > ${prefix}_on_${assemblyid}_${regions}.txt
   echo "done."
   echo ""

   echo "clean up the intermediate files."
   echo ""
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.tab
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.bed
   rm -f ${assemblyid}_Refseq.${regions}.genenames
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.average
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.matrix
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.genenames
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}.header.txt
   echo "done."

fi

echo "Performing PCA analysis ...."
echo ""
Rscript --vanilla ${SH_SCRIPT_DIR}/../R/PCA_plot.R ${prefix}_on_${assemblyid}_${regions}.txt > /dev/null 2>&1
echo "done."
echo ""

mv KAS-seq_PCA_plot.png ${prefix}_KAS-seq_PCA_plot.png
mv KAS-seq_PCA_plot.svg ${prefix}_KAS-seq_PCA_plot.svg

mv KAS-seq_percentage_of_variances.png ${prefix}_KAS-seq_percentage_of_variances.png
mv KAS-seq_percentage_of_variances.svg ${prefix}_KAS-seq_percentage_of_variances.svg

rm -f $KASseq_list
rm -f ${prefix}.labels_basename.txt
rm -f ${prefix}.KAS-seq.RPKM.bigWig.txt
# rm -f ${prefix}_on_${assemblyid}_${regions}.txt

echo "'KAS-Analyzer PCA' run successfully!"
