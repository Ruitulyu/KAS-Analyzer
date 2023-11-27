#!/bin/bash
# 'KAS-Analyzer index' was developed by Ruitu Lyu on 1-22-2022.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-Analyzer index [ -h/--help ] [ -o prefix ] [ -t threads ] [ -s assembly id ] [ -i index types ] [ -l labels ] [ -k KAS-seq ]"
exampleHelp="Example: nohup KAS-Analyzer index -o KAS-seq_pausing_index -t 10 -s hg19 -i pausing -l labels.txt -k KAS-seq.txt &"
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-Analyzer index' output files. REQUIRED."
threadsHelp="-t [threads]: please specify the number of threads used for the pausing or termination index calculation. DEFAULT: 1."
assemblyidHelp="-s [assembly id]: please specify the genome assembly id of KAS-seq data. -s [assembly id]. e.g. Human: hg18, hg19, hg38, hs1; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
indexHelp="-i [index]: please specify the index type you want to calculate. e.g. pausing or termination. REQUIRED."
labelsHelp="-l [labels]: please input the text file containing the labels of (sp)KAS-seq that show in 'KAS-Analyzer index' output files. Default: basename of KAS-seq files.
Example:
WT_rep1
WT.rep2                        ---labels.txt"
KASseqHelp="-k [KAS-seq]: please input the text file containing indexed (sp)KAS-seq bam files used for pausing or termination index calculation. REQUIRED.
Example:
KAS-seq_WT_rep1.bam
KAS-seq_WT_rep2.bam            ---KAS-seq.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-Analyzer index' shell script is applied to calculate the pausing or termination index."

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
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$indexHelp"
    echo -e ""
    echo -e "$labelsHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer index' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

while getopts 'ho:t:s:i:l:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        o) prefix=$OPTARG ;;
	t) threads=$OPTARG ;;
	s) assemblyid=$OPTARG ;;
        i) index=$OPTARG ;;
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
if test -z $prefix ;then
   echo ""
   echo "Please input the prefix (basename) of 'KAS-Analyzer index' output files. -o [prefix]"
   echo ""
   exit 1
fi

if test -z $assemblyid ;then
   echo ""
   echo "Please specify the reference genome assembly id of (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -s [assembly id]"
   echo ""
   exit 1
fi

# supported assembly id.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "hs1" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
    echo ""
    echo "Error: unsupported assembly id: $assemblyid ; Please specify the reference genome assembly id of (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38, hs1; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11."
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
elif [[ $assemblyid == "hs1" ]] ;then
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

if test -z $index ;then
   echo ""
   echo "Please specify the index you want to calculate. e.g. pausing or termination. -i [index type] "
   echo ""
   exit 1
fi

# check unupported index types.
if [[ $index != "pausing" ]] && [[ $index != "termination" ]] ;then
   echo ""
   echo "Error: unsupported index types: $index ; Please specify the index types. e.g. pausing or termination."
   echo ""
   exit 1
fi    

if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing indexed (sp)KAS-seq bam files used for pausing or termination index calculation. -k [KAS-seq]"
   echo ""
   exit 1
fi

# setup the default parameters
number_of_samples=$( awk 'END {print NR}' $KASseq )

if test -n "$labels" ;then
number_of_labels=$( awk 'END {print NR}' $labels )	
   if [[ $number_of_labels != $number_of_samples ]] ;then
      echo ""
      echo "Error:the number of labels isn't consistent with the number of samples."
      echo ""
      printHelpAndExit 0	    
   fi
fi

if test -z $threads ;then
   threads=1
fi

# if labels is not provied, setup the labels as basename of bigWig files.
if test -z $labels ;then
cat /dev/null > ${prefix}.labels_basename.txt
for ((i=1; i<=${number_of_samples}; i++))
do
sample_selected=$( sed -n ''$i'p' $KASseq )
label_basename=$( basename ${sample_selected} .bam )
echo $label_basename >> ${prefix}.labels_basename.txt
done
labels="${prefix}.labels_basename.txt"
fi 

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# get the list of KAS-seq samples and labels.
cat /dev/null > ${prefix}.KAS-seq.RPKM.bigWig.txt
for ((i=1; i<=${number_of_samples}; i++))
do
sample_selected=$(sed -n ''$i'p' $KASseq)
echo "Generate new index of $sample_selected ..."
echo ""
samtools index $sample_selected
echo "done."
echo ""

echo "Normalize $sample_selected using RPKM and output bigWig file ..."
echo ""
bamCoverage -b $sample_selected --outFileFormat bigwig -p $threads -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed --effectiveGenomeSize $genomesize --normalizeUsing RPKM -o ${sample_selected}.bigWig > /dev/null 2>&1
echo ${sample_selected}.bigWig >> ${prefix}.KAS-seq.RPKM.bigWig.txt
echo "done."
echo ""
done

# get the list of KAS-seq bigWig files or labels.
KASseq_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' ${prefix}.KAS-seq.RPKM.bigWig.txt)
labels_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels)

if [[ $index == "pausing" ]] ;then
   echo "Calculating KAS-seq density on ${assemblyid} Refseq promoter or genebody ..."
   echo ""
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.promoter.bed --labels $labels_list -p $threads -out ${prefix}_on_${assemblyid}_Refseq.promoter.npz --outRawCounts ${prefix}_on_${assemblyid}_Refseq.promoter.tab > /dev/null 2>&1
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.genebody.bed --labels $labels_list -p $threads -out ${prefix}_on_${assemblyid}_Refseq.genebody.npz --outRawCounts ${prefix}_on_${assemblyid}_Refseq.genebody.tab > /dev/null 2>&1

   sed "s/nan/0/g" ${prefix}_on_${assemblyid}_Refseq.promoter.tab | sed "1d" | sortBed -i > ${prefix}_on_${assemblyid}_Refseq.promoter.bed
   sortBed -i ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.promoter.bed | awk '{printf("%s\t%s\n",$4,$6)}' > ${assemblyid}_Refseq.promoter.genenames

   sed "s/nan/0/g" ${prefix}_on_${assemblyid}_Refseq.genebody.tab | sed "1d" | sortBed -i > ${prefix}_on_${assemblyid}_Refseq.genebody.bed
   sortBed -i ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.genebody.bed | awk '{printf("%s\t%s\n",$4,$6)}' > ${assemblyid}_Refseq.genebody.genenames
   echo "done."
   echo ""

   echo "Filtering matrix of KAS density on ${assemblyid} Refseq promoter ..."
   echo ""
   paste ${assemblyid}_Refseq.promoter.genenames ${prefix}_on_${assemblyid}_Refseq.promoter.bed | sort -k 1 -n -r | sort -u -k1,1 -k2,2 > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed
   
   paste ${assemblyid}_Refseq.genebody.genenames ${prefix}_on_${assemblyid}_Refseq.genebody.bed | sort -k 1 -n -r | sort -u -k1,1 -k2,2 > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   echo "done."
   echo ""
   
   for ((i=1; i<=${number_of_samples}; i++))
   do
   echo "Calculate the pausing index on ${assemblyid} Refseq genes for the ${i}th sample."
   awk -v x=$i '$(x+5)>=2 {printf("%s\t%s\t%.2f\n",$1,$2,$(x+5))}' ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i}
   awk -v x=$i '{printf("%s\t%s\t%.2f\n",$1,$2,$(x+5))}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}

   awk 'NR==FNR{a[$1]=$3;}NR!=FNR{print $0"\t"a[$1]}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i} ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i} | awk '{ if(NF==3){printf("%s\t%s\t%.2f\n",$1,$2,"1")} else{printf("%s\t%s\t%.2f\n",$1,$2,$3/($4+0.1))} }' > ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i}

   echo -e "genename\tstrand" > ${prefix}.header1.txt
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels | awk -v x=$i '{printf("%s\n",$x)}' > ${prefix}.${i}.header2.txt
   paste ${prefix}.header1.txt ${prefix}.${i}.header2.txt > ${prefix}.${i}.header.txt

   labels_selected=$(sed -n ''$i'p' $labels)
   
   cat ${prefix}.${i}.header.txt ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i} > ${prefix}_${labels_selected}_on_${assemblyid}_Refseq_gene_pausing_index.txt

   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i}
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}
   rm -f ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i}
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.${i}.header2.txt
   rm -f ${prefix}.${i}.header.txt

   echo "done."
   echo ""
   done

   echo "Clean up the intermediate files."
   echo ""
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.tab
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.tab
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.bed
   rm -f ${assemblyid}_Refseq.promoter.genenames
   rm -f ${assemblyid}_Refseq.genebody.genenames
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   echo "done."
   echo ""

elif [[ $index == "termination" ]] ;then
   echo "Calculating KAS-seq density on ${assemblyid} Refseq terminator or genebody ..."
   echo ""
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.terminator.bed --labels $labels_list -p $threads -out ${prefix}_on_${assemblyid}_Refseq.terminator.npz --outRawCounts ${prefix}_on_${assemblyid}_Refseq.terminator.tab > /dev/null 2>&1
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.genebody.bed --labels $labels_list -p $threads -out ${prefix}_on_${assemblyid}_Refseq.genebody.npz --outRawCounts ${prefix}_on_${assemblyid}_Refseq.genebody.tab > /dev/null 2>&1

   sed "s/nan/0/g" ${prefix}_on_${assemblyid}_Refseq.terminator.tab | sed "1d" | sortBed -i > ${prefix}_on_${assemblyid}_Refseq.terminator.bed
   sortBed -i ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.terminator.bed | awk '{printf("%s\t%s\n",$4,$6)}' > ${assemblyid}_Refseq.terminator.genenames

   sed "s/nan/0/g" ${prefix}_on_${assemblyid}_Refseq.genebody.tab | sed "1d" | sortBed -i > ${prefix}_on_${assemblyid}_Refseq.genebody.bed
   sortBed -i ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.genebody.bed | awk '{printf("%s\t%s\n",$4,$6)}' > ${assemblyid}_Refseq.genebody.genenames
   echo "done."
   echo ""   

   echo "Filtering matrix of KAS density on ${assemblyid} Refseq genebody ..."
   paste ${assemblyid}_Refseq.terminator.genenames ${prefix}_on_${assemblyid}_Refseq.terminator.bed | sort -k 1 -n -r | sort -u -k1,1 -k2,2 > ${prefix}_on_${assemblyid}_Refseq.terminator.filter.bed

   paste ${assemblyid}_Refseq.genebody.genenames ${prefix}_on_${assemblyid}_Refseq.genebody.bed | sort -k 1 -n -r | sort -u -k1,1 -k2,2 > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   echo "done."
   echo ""

   for ((i=1; i<=${number_of_samples}; i++))
   do
   echo "Calculate the pausing index on ${assemblyid} Refseq genes for the ${i}th sample."
   awk -v x=$i '{printf("%s\t%s\t%.2f\n",$1,$2,$(x+5))}' ${prefix}_on_${assemblyid}_Refseq.terminator.filter.bed > ${prefix}_on_${assemblyid}_Refseq.terminator.filter.KAS-seq.${i}
   awk -v x=$i '$(x+5)>=2 {printf("%s\t%s\t%.2f\n",$1,$2,$(x+5))}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}

   # match the genename of terminator and genebody.
   awk 'NR==FNR{a[$1]=$3;}NR!=FNR{print $0"\t"a[$1]}' ${prefix}_on_${assemblyid}_Refseq.terminator.filter.KAS-seq.${i} ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i} | awk '{ if(NF==3){printf("%s\t%s\t%s\n",$1,$2,"nan")} else{printf("%s\t%s\t%.2f\n",$1,$2,$4/$3)} }' > ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i}
 
   echo -e "genename\tstrand" > ${prefix}.header1.txt
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels | awk -v x=$i '{printf("%s\n",$x)}' > ${prefix}.${i}.header2.txt
   paste ${prefix}.header1.txt ${prefix}.${i}.header2.txt > ${prefix}.${i}.header.txt

   labels_selected=$(sed -n ''$i'p' $labels)

   cat ${prefix}.${i}.header.txt ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i} > ${prefix}_${labels_selected}_on_${assemblyid}_Refseq_gene_termination_index.txt

   rm -f ${prefix}_on_${assemblyid}_Refseq.terminator.filter.KAS-seq.${i}
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}
   rm -f ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i}
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.${i}.header2.txt
   rm -f ${prefix}.${i}.header.txt

   echo "done."
   echo ""
   done

   echo "Clean up the intermediate files."
   echo ""
   rm -f ${prefix}_on_${assemblyid}_Refseq.terminator.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.terminator.tab
   rm -f ${prefix}_on_${assemblyid}_Refseq.terminator.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.tab
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.bed
   rm -f ${assemblyid}_Refseq.terminator.genenames
   rm -f ${assemblyid}_Refseq.genebody.genenames
   rm -f ${prefix}_on_${assemblyid}_Refseq.terminator.filter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   echo "done."
   echo ""

fi 

rm -f $KASseq_list
rm -f ${prefix}.labels_basename.txt
rm -f ${prefix}.KAS-seq.RPKM.bigWig.txt

echo "'KAS-Analyzer index' run successfully."
