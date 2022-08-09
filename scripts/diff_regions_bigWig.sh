#!/bin/bash
# 'KAS-pipe2 diff' was developed by Ruitu Lyu on 1-20-2022.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 diff [ -h/--help ] [ -t threads ] [ -o prefix ] [ -s assembly id ] [ -r regions ] [ -p peaks ] [ -f fold change ] [ -c comparison file ] [ -l labels ] [ -k KAS-seq ] "
exampleHelp="Example: nohup KAS-pipe2 diff -o KAS-seq_diff -t 10 -s mm10 -r gene -c comparision.txt -l labels.txt -k KAS-seq_data.txt &"
threadsHelp="-t [threads]: please specify the number of threads used for calculating KAS index. DEFAULT: 1."
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 diff' output files. DEFAULT: basename of KAS-seq txt file."
assemblyidHelp="-s [assembly id]: please specify the genome assembly id of KAS-seq data. -s [assembly id]. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce
11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
regionsHelp="-r [regions]: please specify the types of genomic regions used for differential KAS-seq analysis. e.g. promoter, genebody, gene, bin or peak. REQUIRED."
peaksHelp="-p [peaks]: please specify the custom regions file used for differential KAS-seq analysis, if regions type is set as 'peak'. OPTIONAL."
foldchangeHelp="-f [fold change]: please specify the fold change used for identify regions with differential KAS-seq density. DEFAULT: 2."
comparisonHelp="-c [comparison]: please input the text file containing the comparison information for differential KAS-seq analysis. REQUIRED.
	condition
WT.rep1	WT
WT.rep2	WT
KO.rep1	KO
KO.rep2	KO                  ---comparison.txt"
labelsHelp="-l [labels]: please input the text file containing the labels of (sp)KAS-seq data that show in output files, which need to be consistent with comparision file. REQUIRED.
Example:
WT.rep1
WT.rep2
KO.rep1
KO.rep2                     ---labels.txt"
KASseqHelp="-k [KAS-seq]: please input the text file containing the bigWig files of (sp)KAS-seq data that used for differential KAS-seq analysis. REQUIRED.
Example:
KAS-seq_WT.rep1.bigWig
KAS-seq_WT.rep2.bigWig
KAS-seq_KO.rep1.bigWig
KAS-seq_KO.rep2.bigWig      ---KAS-seq_data.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 diff' shell script is applied to perform differential KAS-seq analysis on promoter, genebody, gene, bin or custom regions."

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$threadsHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$regionsHelp"
    echo -e ""
    echo -e "$peaksHelp"
    echo -e ""
    echo -e "$foldchangeHelp"
    echo -e ""
    echo -e "$comparisonHelp"
    echo -e ""
    echo -e "$labelsHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 diff' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the options of 'KAS-pipe2 diff' shell script.
while getopts 'ht:o:s:r:p:f:c:l:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        t) threads=$OPTARG ;;
        o) prefix=$OPTARG ;;
        s) assemblyid=$OPTARG ;;
        r) regions=$OPTARG ;;
        p) peaks=$OPTARG ;;
	f) foldchange=$OPTARG ;;
	c) comparison=$OPTARG ;;
	l) labels=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# check bedtools, samtools and deeptools were installed in your system.
if ! type bedtools > /dev/null 2>&1 ;then
   echo "bedtools was not installed or not export to the \$PATH'"
   echo ""
   echo "Install bedtools with 'conda install -c bioconda bedtools' or refer the official website of 'bedtools'."
   echo ""
   exit 1
fi

if ! type samtools > /dev/null 2>&1 ;then
   echo "samtools was not installed or not export to the \$PATH'"
   echo ""
   echo "Install samtools with 'conda install -c bioconda samtools' or refer the official website of 'samtools'."
   echo ""
   exit 1
fi

if ! type deeptools > /dev/null 2>&1 ;then
   echo "deeptools was not installed or not export to the \$PATH'"
   echo ""
   echo "Install deeptools with 'conda install -c bioconda deeptools' or refer the official website of 'deeptools'."
   echo ""
   exit 1
fi

# Required options.
if test -z $assemblyid ;then
   echo ""
   echo "Please specify the reference genome assembly id of (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6;
 Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -s [assembly id]"
   echo ""
   exit -1
fi

# supported assembly id.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39
" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce
11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
    echo ""
    echo "Error: unsupported assembly id: $assemblyid; Please specify the reference genome assembly id of (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C
.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11."
    echo ""
    exit -1
fi

if test -z $regions ;then
    echo ""
    echo "Please specify the type of genomic regions used for differential KAS-seq analysis. e.g. promoter, genebody, gene, bin or peak. REQUIRED. -r [regions]"
    echo ""
    exit -1
fi

# supported regions types.
if [[ $regions != "promoter" ]] && [[ $regions != "genebody" ]] && [[ $regions != "bin" ]] && [[ $regions != "gene" ]] && [[ $regions != "peak" ]] ;then
   echo ""
   echo "Error: unsupported assembly id: $assemblyid ; Please specify the types of genomic regions. e.g. promoter, genebody, gene, bin or peak."
   echo ""
   exit -1
fi

if [[ $regions == "peak" ]] && test -z $peaks ;then
   echo ""
   echo "Please specify the custom regions file used for KAS density calculation. -p [peaks]"
   echo ""
   exit -1
fi

if test -z $comparison ;then
   echo ""
   echo "Please input the text file containing the comparison information for differential KAS-seq analysis."
   echo ""
   exit -1
fi   

if test -z $labels ;then
   echo ""
   echo "Please input the text file containing the labels of KAS-seq data, which need to be consistent with comparision file."
   echo ""
   exit -1
fi   

if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing bigWig files of (sp)KAS-seq data that used for differential KAS-seq analysis. -k [KAS-seq] "
   echo ""
   exit -1
fi

# test if the number of labels is same to the number of samples.
number_of_samples=$( awk 'END {print NR}' $KASseq )
number_of_labels=$( awk 'END {print NR}' $labels )

if [[ $number_of_labels != $number_of_samples ]] ;then
   echo ""
   echo "Error:the number of labels isn't consistent with the number of samples."
   echo ""
   exit -1

elif [[ $number_of_samples -lt 2 ]] ;then
   echo ""
   echo "Error: the number of samples must be greater than 2."
   echo "" 
   exit -1
fi

# setup the default options.
if test -z $threads ;then
   threads=1
fi

if test -z $foldchange ;then
   foldchange=2
fi

if test -z $prefix ;then
   prefix=$( basename ${KASSeq} .txt )
fi

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# get the list of KASseq bigWig files and labels.
KASseq_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $KASseq)
labels_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels)

# perform differential KAS-seq analysis for 1kb bins.
if [[ $regions == "bin" ]]; then
   echo "Generating the KAS-seq density matrix on ${assemblyid} 1kb bins..."
   echo ""
   multiBigwigSummary bins -b $KASseq_list --labels $labels_list --binSize 1000 -p $threads --blackListFileName ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -out ${prefix}_on_1kb_bins.npz --outRawCounts ${prefix}_on_1kb_bins.tab > /dev/null 2>&1
   sed "s/nan/0/g" ${prefix}_on_1kb_bins.tab | sed "1d" > ${prefix}_on_1kb_bins.bed
   echo "done."
   echo ""

   #calculate the average value of very single row in the table. 
   echo "Calculating the average KAS-seq signals on ${assemblyid} 1kb bins ..."
   echo ""
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print  a[j]/(NF-3) }' ${prefix}_on_1kb_bins.bed > ${prefix}_on_1kb_bins.average
   echo "done."
   echo ""

   #filter the bins with averaged KAS-seq density lower than 5.
   echo "Filter KAS density matrix on 1kb bins ..."
   echo ""
   paste ${prefix}_on_1kb_bins.average ${prefix}_on_1kb_bins.bed | awk '$1>=5 {print $0}' > ${prefix}_on_1kb_bins.filter.bed
   cut -f1,2,3,4 --complement ${prefix}_on_1kb_bins.filter.bed | awk '{for(i=1;i<=NF;i++){printf "%d\t", $i*10}; printf "\n"}' > ${prefix}_on_${assemblyid}_${regions}.matrix
   awk '{printf("%s\n","bin"FNR"linker"$2"linker"$3"linker"$4)}' ${prefix}_on_1kb_bins.filter.bed > ${prefix}_on_${assemblyid}_${regions}.rowname
   echo "done."
   echo ""

   echo "Generating the KAS-seq signal interger matrix with binname and header ..."
   echo ""
   paste ${prefix}_on_${assemblyid}_${regions}.rowname ${prefix}_on_${assemblyid}_${regions}.matrix > ${prefix}_on_${assemblyid}_${regions}.without_header.txt
   echo "done."
   echo ""

   if [[ $number_of_samples == 2 ]] ;then
      echo "KAS-seq replicates aren't provided, DEseq2 will not be used for differential KAS-seq analysis. Fold change will be used ..."
      echo ""
      awk -v x=$foldchange '($3+0.01)/($2+0.01)>=x {printf("%s\t%.1f\t%.1f\t%.2f\n",$1,$2/10,$3/10,log(($3+0.01)/($2+0.01))/log(2))}' ${prefix}_on_${assemblyid}_${regions}.without_header.txt > ${prefix}_gain_of_KAS-seq.without_header.txt
      awk -v x=$foldchange '($2+0.01)/($3+0.01)>=x {printf("%s\t%.1f\t%.1f\t%.2f\n",$1,$2/10,$3/10,log(($3+0.01)/($2+0.01))/log(2))}' ${prefix}_on_${assemblyid}_${regions}.without_header.txt > ${prefix}_loss_of_KAS-seq.without_header.txt
      sed -i "s/linker/\t/g" ${prefix}_gain_of_KAS-seq.without_header.txt
      sed -i "s/linker/\t/g" ${prefix}_loss_of_KAS-seq.without_header.txt

      echo -e "bins\tchr\tstart\tend" > ${prefix}.header1.txt
      awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' $labels > ${prefix}.header2.txt
      echo -e "log2foldchange" > ${prefix}.header3.txt
      paste ${prefix}.header1.txt ${prefix}.header2.txt ${prefix}.header3.txt > ${prefix}.header.txt

      cat ${prefix}.header.txt ${prefix}_gain_of_KAS-seq.without_header.txt > ${prefix}_gain_of_KAS-seq_1kb_bins.txt
      cat ${prefix}.header.txt ${prefix}_loss_of_KAS-seq.without_header.txt > ${prefix}_loss_of_KAS-seq_1kb_bins.txt
    
      echo "done."
      echo ""

   elif [[ $number_of_samples -ge 4 ]] ;then
      awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' $labels > ${prefix}.header.txt
      cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_${regions}.without_header.txt > ${prefix}_on_${assemblyid}_${regions}.txt

      echo "Performing differential KAS-seq analysis using DEseq2 ...."
      echo ""
      Rscript --vanilla ${SH_SCRIPT_DIR}/../R/DESeq2_diff_KAS-seq.R ${prefix}_on_${assemblyid}_${regions}.txt $comparison $foldchange
      echo "done."
      echo ""

      sed "s/\"//g" KAS-seq_DESeq2_output.csv | sed "s/\,/\t/g" > ${prefix}_KAS-seq_1kb_bins_DESeq2_output.txt
      sed "s/\"//g" DE.KAS-seq_DESeq2_Fold_padj0.05_output.csv | sed "s/\,/\t/g" > ${prefix}_DE.KAS-seq_1kb_bins_DESeq2_fold${foldchange}_padj0.05_output.txt

      rm -f KAS-seq_DESeq2_output.csv
      rm -f DE.KAS-seq_DESeq2_Fold_padj0.05_output.csv
   fi

   # clean up the intermediate files.    
   echo "clean up the intermediate files"
   echo ""
   rm -f ${prefix}_on_1kb_bins.tab
   rm -f ${prefix}_on_1kb_bins.npz
   rm -f ${prefix}_on_1kb_bins.bed
   rm -f ${prefix}_on_1kb_bins.average
   rm -f ${prefix}_on_1kb_bins.filter.bed
   rm -f ${prefix}_on_${assemblyid}_${regions}.matrix
   rm -f ${prefix}_on_${assemblyid}_${regions}.rowname
   rm -f ${prefix}_on_${assemblyid}_${regions}.without_header.txt
   rm -f ${prefix}.header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}.header3.txt
   rm -f ${prefix}_gain_of_KAS-seq.without_header.txt
   rm -f ${prefix}_loss_of_KAS-seq.without_header.txt
   echo "done."
   echo ""

elif [[ $regions == "peak" ]]; then
   echo "Generating the KAS-seq matrix on ${peaks}..."
   echo ""
   peaks_basename=$( basename ${peaks} .bed )
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED $peaks --labels $labels_list -p $threads -out ${prefix}_on_${peaks_basename}.npz --outRawCounts ${prefix}_on_${peaks_basename}.tab > /dev/null 2>&1
   sed "s/nan/0/g" ${prefix}_on_${peaks_basename}.tab | sed "1d" > ${prefix}_on_${peaks_basename}.bed
   echo "done."
   echo ""

   echo "Generating the KAS-seq signal interger matrix with rowname and header ..."
   echo ""
   # generate the KAS-seq signals integer matrix multiply by 10.
   cut -f1,2,3 --complement ${prefix}_on_${peaks_basename}.bed | awk '{for(i=1;i<=NF;i++){printf "%d\t", $i*10}; printf "\n"}' > ${prefix}_on_${peaks_basename}.matrix
   awk '{printf("%s\n","peak"NR"linker"$1"linker"$2"linker"$3)}' ${prefix}_on_${peaks_basename}.bed > ${prefix}_on_${peaks_basename}.rowname
   # generating the KAS-seq signals matrix without header.
   paste ${prefix}_on_${peaks_basename}.rowname ${prefix}_on_${peaks_basename}.matrix > ${prefix}_on_${assemblyid}_${regions}.without_header.txt
   echo "done."
   echo ""

   if [[ $number_of_samples == 2 ]] ;then
      echo "KAS-seq replicates aren't provided, DEseq2 will not be used for differential KAS-seq analysis. Fold change will be used ..."
      echo ""
      awk -v x=$foldchange '($3+0.01)/($2+0.01)>=x {printf("%s\t%.1f\t%.1f\t%.2f\n",$1,$2/10,$3/10,log(($3+0.01)/($2+0.01))/log(2))}' ${prefix}_on_${assemblyid}_${regions}.without_header.txt > ${prefix}_gain_of_KAS-seq.without_header.txt
      awk -v x=$foldchange '($2+0.01)/($3+0.01)>=x {printf("%s\t%.1f\t%.1f\t%.2f\n",$1,$2/10,$3/10,log(($3+0.01)/($2+0.01))/log(2))}' ${prefix}_on_${assemblyid}_${regions}.without_header.txt > ${prefix}_loss_of_KAS-seq.without_header.txt
      sed -i "s/"linker"/\t/g" ${prefix}_gain_of_KAS-seq.without_header.txt
      sed -i "s/"linker"/\t/g" ${prefix}_loss_of_KAS-seq.without_header.txt

      echo -e "peaks\tchr\tstart\tend" > ${prefix}.header1.txt
      awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' $labels > ${prefix}.header2.txt
      echo -e "log2foldchange" > ${prefix}.header3.txt
      paste ${prefix}.header1.txt ${prefix}.header2.txt ${prefix}.header3.txt > ${prefix}.header.txt

      cat ${prefix}.header.txt ${prefix}_gain_of_KAS-seq.without_header.txt > ${prefix}_gain_of_KAS-seq_${regions}.txt
      cat ${prefix}.header.txt ${prefix}_loss_of_KAS-seq.without_header.txt > ${prefix}_loss_of_KAS-seq_${regions}.txt

      echo "done."
      echo ""

   elif [[ $number_of_samples -ge 4 ]] ;then
      awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' $labels > ${prefix}.header.txt
      cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_${regions}.without_header.txt > ${prefix}_on_${assemblyid}_${regions}.txt

      echo "Performing differential KAS-seq analysis using DEseq2 ..."
      echo ""
      Rscript --vanilla ${SH_SCRIPT_DIR}/../R/DESeq2_diff_KAS-seq.R ${prefix}_on_${assemblyid}_${regions}.txt $comparison $foldchange
      echo "done."
      echo ""

      sed "s/\"//g" KAS-seq_DESeq2_output.csv | sed "s/\,/\t/g" > ${prefix}_KAS-seq_peaks_DESeq2_output.txt
      sed "s/\"//g" DE.KAS-seq_DESeq2_Fold_padj0.05_output.csv | sed "s/\,/\t/g" > ${prefix}_DE.KAS-seq_peaks_DESeq2_fold${foldchange}_padj0.05_output.txt

      rm -f KAS-seq_DESeq2_output.csv
      rm -f DE.KAS-seq_DESeq2_Fold_padj0.05_output.csv
   fi

   # clean up the intermediate files.
   echo "Clean up the intermediate files ..."
   echo ""
   rm -f ${prefix}_on_${peaks_basename}.npz
   rm -f ${prefix}_on_${peaks_basename}.tab
   rm -f ${prefix}_on_${peaks_basename}.bed
   rm -f ${prefix}_on_${peaks_basename}.matrix 
   rm -f ${prefix}_on_${peaks_basename}.rowname
   rm -f ${prefix}_on_${assemblyid}_${regions}.without_header.txt
   rm -f ${prefix}.header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}.header3.txt
   rm -f ${prefix}_gain_of_KAS-seq.without_header.txt
   rm -f ${prefix}_loss_of_KAS-seq.without_header.txt

   echo "done."
   echo ""

elif [[ $regions == "promoter" ]] || [[ $regions == "genebody" ]] ;then

   echo "Generating the KAS-seq matrix on ${assemblyid} Refseq ${regions} ..."
   echo ""
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.${regions}.bed --labels $labels_list -p $threads -out ${prefix}_on_${assemblyid}_Refseq.${regions}.npz --outRawCounts ${prefix}_on_${assemblyid}_Refseq.${regions}.tab

   sed "s/nan/0/g" ${prefix}_on_${assemblyid}_Refseq.${regions}.tab | sed "1d" | sortBed -i > ${prefix}_on_${assemblyid}_Refseq.${regions}.bed
   sortBed -i ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.${regions}.bed | awk '{printf("%s\n",$1"--"$2"--"$3"--"$4"--"$6)}' > ${assemblyid}_Refseq.${regions}.genenames
   echo "done."
   echo ""

   # calculate the average KAS index of every single row in the table.
   echo "Calculating the average KAS-seq signals on ${assemblyid}_Refseq.${regions}.bed ..."
   echo ""
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/(NF-3) }' ${prefix}_on_${assemblyid}_Refseq.${regions}.bed > ${prefix}_on_${assemblyid}_Refseq.${regions}.average
   echo "done."
   echo ""

# filter the ${regions} with normalized KAS-seq density over 5.
   echo "Filtering matrix of KAS signal on ${assemblyid}_Refseq.${regions}.bed ..."
   echo ""
   paste ${prefix}_on_${assemblyid}_Refseq.${regions}.average ${assemblyid}_Refseq.${regions}.genenames ${prefix}_on_${assemblyid}_Refseq.${regions}.bed | awk '$1>=5 {print $0}' > ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed
   echo "done."
   echo ""

   echo "Generating the KAS-seq signal interger matrix with genename ..."
   echo ""
   cut -f1,2,3,4,5 --complement ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed | awk '{for(i=1;i<=NF;i++){printf "%d\t", $i*10}; printf "\n"}' > ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.matrix
   awk '{printf("%s\n",$2)}' ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed > ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.genenames

   paste ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.genenames ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.matrix > ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt
   echo "done."
   echo ""

   if [[ $number_of_samples == 2 ]] ;then
       echo "KAS-seq replicates aren't provided, DEseq2 will not be used for differential KAS-seq analysis. Fold change will be used ..."
       echo ""
       awk -v x=$foldchange '($3+0.01)/($2+0.01)>=x {printf("%s\t%.1f\t%.1f\t%.2f\n",$1,$2/10,$3/10,log(($3+0.01)/($2+0.01))/log(2))}' ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt > ${prefix}_gain_of_KAS-seq.without_header.txt
       awk -v x=$foldchange '($2+0.01)/($3+0.01)>=x {printf("%s\t%.1f\t%.1f\t%.2f\n",$1,$2/10,$3/10,log(($3+0.01)/($2+0.01))/log(2))}' ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt > ${prefix}_loss_of_KAS-seq.without_header.txt
       sed -i "s/"--"/\t/g" ${prefix}_gain_of_KAS-seq.without_header.txt
       sed -i "s/"--"/\t/g" ${prefix}_loss_of_KAS-seq.without_header.txt

       echo -e "chr\tstart\tend\tgene\tstrand" > ${prefix}.header1.txt
       awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' $labels > ${prefix}.header2.txt
       echo -e "log2foldchange" > ${prefix}.header3.txt
       paste ${prefix}.header1.txt ${prefix}.header2.txt ${prefix}.header3.txt > ${prefix}.header.txt

       cat ${prefix}.header.txt ${prefix}_gain_of_KAS-seq.without_header.txt > ${prefix}_gain_of_KAS-seq_${regions}.txt
       cat ${prefix}.header.txt ${prefix}_loss_of_KAS-seq.without_header.txt > ${prefix}_loss_of_KAS-seq_${regions}.txt

       echo "done."
       echo ""

   elif [[ $number_of_samples -ge 4 ]] ;then	    
       awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' $labels > ${prefix}.header.txt
       cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt > ${prefix}_on_${assemblyid}_${regions}.txt
    
       echo "Performing differential KAS-seq analysis using DEseq2 ..."
       echo ""
       Rscript --vanilla ${SH_SCRIPT_DIR}/../R/DESeq2_diff_KAS-seq.R ${prefix}_on_${assemblyid}_${regions}.txt $comparison $foldchange
       echo "done."
       echo ""

       sed "s/\"//g" KAS-seq_DESeq2_output.csv | sed "s/\,/\t/g" > ${prefix}_KAS-seq_${regions}_DESeq2_output.txt
       sed "s/\"//g" DE.KAS-seq_DESeq2_Fold_padj0.05_output.csv | sed "s/\,/\t/g" > ${prefix}_DE.KAS-seq_${regions}_DESeq2_fold${foldchange}_padj0.05_output.txt

       rm -f KAS-seq_DESeq2_output.csv
       rm -f DE.KAS-seq_DESeq2_Fold_padj0.05_output.csv
   fi  

# clean up the intermediate files.
   echo "clean up the intermediate files"
   echo ""
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.tab
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.bed
   rm -f ${assemblyid}_Refseq.${regions}.genenames
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.average
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.matrix
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.genenames
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}.header3.txt
   rm -f ${prefix}.header.txt
   rm -f ${prefix}_gain_of_KAS-seq.without_header.txt
   rm -f ${prefix}_loss_of_KAS-seq.without_header.txt
   echo "done."
   echo ""

elif [[ $regions == "gene" ]] ;then

   echo "Generating the KAS-seq density matrix on ${assemblyid} Refseq ${regions}..." 
   echo ""
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.promoter.bed --labels $labels_list -p $threads -out ${prefix}_on_${assemblyid}_Refseq.promoter.npz --outRawCounts ${prefix}_on_${assemblyid}_Refseq.promoter.tab > /dev/null 2>&1

   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.genebody.bed --labels $labels_list -p $threads -out ${prefix}_on_${assemblyid}_Refseq.genebody.npz --outRawCounts ${prefix}_on_${assemblyid}_Refseq.genebody.tab > /dev/null 2>&1

   sed "s/nan/0/g" ${prefix}_on_${assemblyid}_Refseq.promoter.tab | sed "1d" | sortBed -i > ${prefix}_on_${assemblyid}_Refseq.promoter.bed
   sortBed -i ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.promoter.bed | awk '{printf("%s\t%s\n",$4,$6)}' > ${assemblyid}_Refseq.promoter.genenames

   sed "s/nan/0/g" ${prefix}_on_${assemblyid}_Refseq.genebody.tab | sed "1d" | sortBed -i > ${prefix}_on_${assemblyid}_Refseq.genebody.bed
   sortBed -i ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.genebody.bed | awk '{printf("%s\t%s\n",$4,$6)}' > ${assemblyid}_Refseq.genebody.genenames
   echo "done."
   echo ""

   # calculate the average KAS expression of every single row in the table.
   echo "Calculate the average KAS density on ${assemblyid} Refseq promoter and genebody ..."
   echo ""
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/(NF-3) }' ${prefix}_on_${assemblyid}_Refseq.promoter.bed > ${prefix}_on_${assemblyid}_Refseq.promoter.average
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/(NF-3) }' ${prefix}_on_${assemblyid}_Refseq.genebody.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.average
   echo "done."
   echo ""

   # filter the promoter or genebody with averaged KAS-seq density lower than 5.
   echo "Filter matrix of KAS density on ${assemblyid} Refseq promoter and genebody ..."
   echo ""
   paste ${prefix}_on_${assemblyid}_Refseq.promoter.average ${assemblyid}_Refseq.promoter.genenames ${prefix}_on_${assemblyid}_Refseq.promoter.bed | sort -k 1 -n -r | sort -u -k2,2 -k3,3 | awk '$1>=5 {print $0}' > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed
   paste ${prefix}_on_${assemblyid}_Refseq.genebody.average ${assemblyid}_Refseq.genebody.genenames ${prefix}_on_${assemblyid}_Refseq.genebody.bed | sort -k 1 -n -r | sort -u -k2,2 -k3,3 > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   echo "done."
   echo ""

   for ((i=1; i<=${number_of_samples}; i++))
   do
   echo "Combine KAS-seq density on promoter and genebody, and calculate the KAS-seq density on ${assemblyid} Refseq genes for the ${i}th sample." 
   echo ""
   awk -v x=$i '{printf("%s\t%d\n",$2"-"$3,$(x+6)*10)}' ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i}
   awk -v x=$i '{printf("%s\t%d\n",$2"-"$3,$(x+6)*10)}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}

   awk 'NR==FNR{a[$1]=$2;}NR!=FNR{print $0"\t"a[$1]}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i} ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i} | awk '{ if(NF==2){printf("%d\n",$2)} else{printf("%d\n",($2+$3)/2)} }' > ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i}

   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i}
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}

   echo "done."
   echo ""
   done

   echo "Generating the final KAS-seq density matrix on ${assemblyid}_Refseq.gene.bed for all samples ..."
   echo ""
   awk '{printf("%s\n",$2"--"$3)}' ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.genename
   paste ${prefix}_on_${assemblyid}_Refseq.promoter.filter.genename ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.* > ${prefix}_on_${assemblyid}_Refseq.genes.without_header.txt

   if [[ $number_of_samples == 2 ]] ;then
       echo "KAS-seq replicates aren't provided, DEseq2 will not be used for differential KAS-seq analysis. Fold change: $foldchange will be used ..."
       echo ""
       awk -v x=$foldchange '($3+0.01)/($2+0.01)>=x {printf("%s\t%.1f\t%.1f\t%.2f\n",$1,$2/10,$3/10,log(($3+0.01)/($2+0.01))/log(2))}' ${prefix}_on_${assemblyid}_Refseq.genes.without_header.txt > ${prefix}_gain_of_KAS-seq_genes.without_header.txt
       awk -v x=$foldchange '($2+0.01)/($3+0.01)>=x {printf("%s\t%.1f\t%.1f\t%.2f\n",$1,$2/10,$3/10,log(($3+0.01)/($2+0.01))/log(2))}' ${prefix}_on_${assemblyid}_Refseq.genes.without_header.txt > ${prefix}_loss_of_KAS-seq_genes.without_header.txt
       sed -i "s/"--"/\t/g" ${prefix}_gain_of_KAS-seq_genes.without_header.txt
       sed -i "s/"--"/\t/g" ${prefix}_loss_of_KAS-seq_genes.without_header.txt

       echo -e "genes\tstrand" > ${prefix}.header1.txt
       awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' $labels > ${prefix}.header2.txt
       echo -e "log2foldchange" > ${prefix}.header3.txt
       paste ${prefix}.header1.txt ${prefix}.header2.txt ${prefix}.header3.txt > ${prefix}.header.txt

       cat ${prefix}.header.txt ${prefix}_gain_of_KAS-seq_genes.without_header.txt > ${prefix}_gain_of_KAS-seq_genes.txt
       cat ${prefix}.header.txt ${prefix}_loss_of_KAS-seq_genes.without_header.txt > ${prefix}_loss_of_KAS-seq_genes.txt

       echo "done."
       echo ""

   elif [[ $number_of_samples -ge 4 ]] ;then
       awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' $labels > ${prefix}.header.txt
       cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_Refseq.genes.without_header.txt > ${prefix}_on_${assemblyid}_${regions}.txt

       echo "Performing differential KAS-seq analysis using DEseq2 ..."
       echo ""
       Rscript --vanilla ${SH_SCRIPT_DIR}/../R/DESeq2_diff_KAS-seq.R ${prefix}_on_${assemblyid}_${regions}.txt $comparison $foldchange
       echo "done."
       echo ""

       sed "s/\"//g" KAS-seq_DESeq2_output.csv | sed "s/\,/\t/g" > ${prefix}_KAS-seq_genes_DESeq2_output.txt
       sed "s/\"//g" DE.KAS-seq_DESeq2_Fold_padj0.05_output.csv | sed "s/\,/\t/g" > ${prefix}_DE.KAS-seq_genes_DESeq2_fold${foldchange}_padj0.05_output.txt

       rm -f KAS-seq_DESeq2_output.csv
       rm -f DE.KAS-seq_DESeq2_Fold_padj0.05_output.csv
   fi

# clean up the intermediate files.
   echo "clean up the intermediate files"
   echo ""
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.tab
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.tab
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.bed
   rm -f ${assemblyid}_Refseq.promoter.genenames
   rm -f ${assemblyid}_Refseq.genebody.genenames
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.average
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.average
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.*
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.filter.genename
   rm -f ${prefix}_on_${assemblyid}_Refseq.genes.without_header.txt
   rm -f ${prefix}_gain_of_KAS-seq_genes.without_header.txt
   rm -f ${prefix}_loss_of_KAS-seq_genes.without_header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}.header3.txt
   rm -f ${prefix}.header.txt
   echo "done."
   echo ""

fi 

echo "'KAS-pipe2 diff' run successfully!"
