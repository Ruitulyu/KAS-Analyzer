#!/bin/bash
# 'KAS-pipe2 KASexpre' was developed by Ruitu Lyu on 1-19-2022.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 KASexpre [ -h/--help ] [ -t threads ] [ -o prefix ] [ -s assembly id ] [ -r regions ] [ -p peaks ] [ -l labels ] [ -k KAS-seq ] "
exampleHelp="Example: nohup KAS-pipe2 KASexpre -o KAS-seq_expression -t 10 -s mm10 -r TSS -l labels.txt -k KAS-seq_data.txt &"
threadsHelp="-t [threads]: please specify the number of threads used for calculating KAS expression. Default: 1."
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 KASexpre' output files. Default: basename of KAS-seq data."
assemblyidHelp="-s [assembly id]: please specify the genome assembly id of KAS-seq data. -s [assembly id]. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
regionsHelp="-r [regions]: please specify the types of genomic regions. e.g. promoter, genebody, gene or peak. REQUIRED."
peaksHelp="-p [peaks]: please specify the custom regions file used for KAS index calculation. OPTIONAL."
labelsHelp="-l [labels]: please input the text file containing the labels of (sp)KAS-seq data that show in output files. REQUIRED.
Example:
WT.rep1
WT.rep2
KO.rep1
KO.rep2                       ---labels.txt"
KASseqHelp="-k [KAS-seq]: please input the text file containing the indexed bam files of (sp)KAS-seq data that used to calculate KAS expression. REQUIRED.
Example:
KAS-seq_WT.rep1.bam
KAS-seq_WT.rep2.bam
KAS-seq_KO.rep1.bam
KAS-seq_KO.rep2.bam           ---KAS-seq_data.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 KASexpre' shell script is applied to calculate normalized KAS-seq expression levels on promoter, genebody, genes or custom regions."

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
    echo -e "$labelsHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 KASindex' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the options of 'KAS-pipe2 KASindex' shell script.
while getopts 'ht:o:s:r:p:l:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        t) threads=$OPTARG ;;
        o) prefix=$OPTARG ;;
        s) assemblyid=$OPTARG ;;
        r) regions=$OPTARG ;;
        p) peaks=$OPTARG ;;
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
   echo "Please specify the reference genome assembly id of (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -s [assembly id]"
   echo ""
   exit 1
fi

# supported assembly id.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
    echo ""
    echo "Error: unsupported assembly id: $assemblyid ; Please specify the reference genome assembly id of (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11."
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

if test -z $regions ;then
    echo ""
    echo "Please specify the types of genomic regions used for KAS expression calculation. e.g. promoter, genebody, gene or peak. REQUIRED. -r [regions]"
    echo ""
    exit 1
fi

# supported regions types.
if [[ $regions != "promoter" ]] && [[ $regions != "genebody" ]] && [[ $regions != "gene" ]] && [[ $regions != "peak" ]] ;then
   echo ""
   echo "Error: unsupported assembly id: $assemblyid ; Please specify the types of genomic regions. e.g. promoter, genebody, gene or peak."
   echo ""
   exit 1
fi

if [[ $regions == "peak" ]] && test -z $peaks ;then
   echo ""
   echo "Please specify the custom regions file used for KAS expression calculation."
   echo ""
   exit 1
fi 

if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing the bigWig files of (sp)KAS-seq data that used to calculate KAS expression. -k [KAS-seq] "
   echo ""
   exit 1
fi   

# setup the default options.
if test -z $threads ;then
   threads=1
fi

if test -z $prefix ;then
   prefix=$( basename ${KASSeq} .txt )
fi

# get the number of samples.
number_of_samples=$( awk 'END {print NR}' $KASseq )

if test -n "$labels" ;then
number_of_labels=$( awk 'END {print NR}' $labels )
   if [[ $number_of_labels != $number_of_samples ]] ;then 
      echo "" 
      echo "Error: the number of labels isn't consistent with the number of samples." 
      echo ""
      exit 1
   fi

else
cat /dev/null > ${prefix}.labels_basename.txt	
for ((i=1; i<=${number_of_samples}; i++))
do
sample_selected=$( sed -n ''$i'p' $KASseq )
label_basename=$( basename ${sample_selected} .bam )
echo $label_basename >> .labels_basename.txt
done
labels="${prefix}.labels_basename.txt"

fi

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

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

if [[ $regions == "promoter" ]] || [[ $regions == "genebody" ]] ;then
   echo "Calculating KAS-seq expression on ${assemblyid} Refseq ${regions} ..."
   echo ""
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.${regions}.bed --labels $labels_list -p $threads -out ${prefix}_on_${assemblyid}_Refseq.${regions}.npz --outRawCounts ${prefix}_on_${assemblyid}_Refseq.${regions}.tab > /dev/null 2>&1

   sed "s/nan/0/g" ${prefix}_on_${assemblyid}_Refseq.${regions}.tab | sed "1d" | sortBed -i > ${prefix}_on_${assemblyid}_Refseq.${regions}.bed
   sortBed -i ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.${regions}.bed | awk '{printf("%s\t%s\n",$4,$6)}' > ${assemblyid}_Refseq.${regions}.genenames
   echo "done."
   echo ""

   # calculate the average KAS expression of every single row in the table.
   echo "Calculating the average KAS-seq expression on ${assemblyid}_Refseq.${regions}.bed ..."
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/NF }' ${prefix}_on_${assemblyid}_Refseq.${regions}.bed > ${prefix}_on_${assemblyid}_Refseq.${regions}.average
   echo ""

   # filter the ${regions} with KAS expression over 0.5.
   echo "Filtering matrix of normalized KAS expression on ${assemblyid}_Refseq.${regions}.bed ..."
   paste ${prefix}_on_${assemblyid}_Refseq.${regions}.average ${assemblyid}_Refseq.${regions}.genenames ${prefix}_on_${assemblyid}_Refseq.${regions}.bed | sort -k 1 -n -r | sort -u -k2,2 -k3,3 | awk '$1>=2 {print $0}' > ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed
   echo "done."
   echo ""

   echo "Generating the final normalized KAS expression on ${assemblyid}_Refseq.${regions}.bed ..."
   cut -f1,2,3,4,5,6 --complement ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed | awk '{for(i=1;i<=NF;i++){printf "%.2f\t", $i}; printf "\n"}' > ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.matrix
   awk '{printf("%s\t%d\t%d\t%s\t%s\n",$4,$5,$6,$2,$3)}' ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed > ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.genes

   paste ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.genes ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.matrix > ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt

   echo -e "chr\tstart\tend\tgenename\tstrand" > ${prefix}.header1.txt
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels > ${prefix}.header2.txt
   paste ${prefix}.header1.txt ${prefix}.header2.txt > ${prefix}.header.txt

   cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt > ${prefix}_on_${assemblyid}_Refseq_${regions}_KAS-seq_expression.txt
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
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.genes
   rm -f ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}.header.txt
   echo "done."
   echo ""

elif [[ $regions == "gene" ]] ;then

   echo "Calculate KAS expression on ${assemblyid} Refseq genes."
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
   echo "Calculate the average KAS expression on ${assemblyid} Refseq promoter and genebody ..."
   echo ""
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/NF }' ${prefix}_on_${assemblyid}_Refseq.promoter.bed > ${prefix}_on_${assemblyid}_Refseq.promoter.average
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/NF }' ${prefix}_on_${assemblyid}_Refseq.genebody.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.average
   echo "done."
   echo ""

   # filter the promoter or genebody with averaged KAS expression lower than 0.5.
   echo "Filter matrix of normalized KAS expression on ${assemblyid} Refseq promoter and genebody ..."
   paste ${prefix}_on_${assemblyid}_Refseq.promoter.average ${assemblyid}_Refseq.promoter.genenames ${prefix}_on_${assemblyid}_Refseq.promoter.bed | sort -k 1 -n -r | sort -u -k2,2 -k3,3 | awk '$1>=2 {print $0}' > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed

   paste ${prefix}_on_${assemblyid}_Refseq.genebody.average ${assemblyid}_Refseq.genebody.genenames ${prefix}_on_${assemblyid}_Refseq.genebody.bed | sort -k 1 -n -r | sort -u -k2,2 -k3,3 > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   echo "done."
   echo ""

   for ((i=1; i<=${number_of_samples}; i++))
   do
   echo "Calculate the normalized KAS-seq expression on ${assemblyid} Refseq genes for the ${i}th sample."	
   awk -v x=$i '{printf("%s\t%.2f\n",$2,$(x+6))}' ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i}
   awk -v x=$i '{printf("%s\t%.2f\n",$2,$(x+6))}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}

   # match the genename of promoter and genebody matrix
   awk 'NR==FNR{a[$1]=$2;}NR!=FNR{print $0"\t"a[$1]}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i} ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i} | awk '{ if(NF==2){printf("%.2f\n",$2)} else{printf("%.2f\n",($2+$3)/2)} }' > ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i}

   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i}
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}

   echo "done."
   echo ""
   done

   echo "Generating the final KAS expression on ${assemblyid}_Refseq.gene.bed for all samples ..."
   awk '{printf("%s\t%s\n",$2,$3)}' ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.genename
   paste ${prefix}_on_${assemblyid}_Refseq.promoter.filter.genename ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.* > ${prefix}_on_${assemblyid}_Refseq.gene.filter.without_header.txt

   echo -e "genename\tstrand" > ${prefix}.header1.txt
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels > ${prefix}.header2.txt
   paste ${prefix}.header1.txt ${prefix}.header2.txt > ${prefix}.header.txt

   cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_Refseq.gene.filter.without_header.txt > ${prefix}_on_${assemblyid}_Refseq_gene_KAS-seq_expression.txt
   echo "done."
   echo ""

   echo "Clean up the intermediate files."
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.tab
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.tab
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.bed
   rm -f ${assemblyid}_Refseq.promoter.genenames
   rm -f ${assemblyid}_Refseq.genebody.genenames
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.average
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.average
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.filter.genename
   rm -f ${prefix}_on_${assemblyid}_Refseq.gene.filter.without_header.txt
   rm -f ${prefix}.header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt

   echo "done."
   echo ""

elif [[ $regions == "peak" ]]; then

   echo "Calculate KAS expression on ${peaks}."
   echo ""
   peaks_basename=$( basename ${peaks} .bed )
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${peaks} --labels $labels_list -p $threads -out ${prefix}_on_${peaks_basename}.npz --outRawCounts ${prefix}_on_${peaks_basename}.tab > /dev/null 2>&1
   sed "s/nan/0/g" ${prefix}_on_${peaks_basename}.tab | sed "1d" | sortBed -i > ${prefix}_on_${peaks_basename}.bed
   echo "done."
   echo ""

   echo "Generate the final KAS expression on ${peaks} ..."
   cut -f1,2,3 --complement ${prefix}_on_${peaks_basename}.bed | awk '{for(i=1;i<=NF;i++){printf "%.2f\t", $i}; printf "\n"}' > ${prefix}_on_${peaks_basename}.matrix
   awk '{printf("%s\t%d\t%d\t%s\n",$1,$2,$3,"peak"FNR)}' ${prefix}_on_${peaks_basename}.bed > ${prefix}_on_${peaks_basename}.4bed
 
   paste ${prefix}_on_${peaks_basename}.4bed ${prefix}_on_${peaks_basename}.matrix > ${prefix}_on_${peaks_basename}.without_header.txt
   echo -e "chr\tstart\tend\tpeaks" > ${prefix}.header1.txt
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels > ${prefix}.header2.txt
   paste ${prefix}.header1.txt ${prefix}.header2.txt > ${prefix}.header.txt

   cat ${prefix}.header.txt ${prefix}_on_${peaks_basename}.without_header.txt > ${prefix}_on_${peaks_basename}_KAS-seq_expression.txt
   echo "done."
   echo ""

   echo "Clean up the intermediate files."
   rm -f ${prefix}.header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}_on_${peaks_basename}.npz
   rm -f ${prefix}_on_${peaks_basename}.tab
   rm -f ${prefix}_on_${peaks_basename}.bed
   rm -f ${prefix}_on_${peaks_basename}.matrix
   rm -f ${prefix}_on_${peaks_basename}.4bed
   rm -f ${prefix}_on_${peaks_basename}.without_header.txt
   echo "done."
   echo ""

fi 

rm -f $KASseq_list
rm -f ${prefix}.labels_basename.txt
rm -f ${prefix}.KAS-seq.RPKM.bigWig.txt

echo "'KAS-pipe2 KASexpre' run successfully!"
