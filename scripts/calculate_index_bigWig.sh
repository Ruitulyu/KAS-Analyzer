#!/bin/bash
# 'KAS-pipe2 index' was developed by Ruitu Lyu on 1-22-2022.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 index [ -h/--help ] [ -o prefix ] [ -t threads ] [ -s assembly id ] [ -i index types ] [ -l labels ] [ -k KAS-seq ]"
exampleHelp="Example: nohup KAS-pipe2 index -o KAS-seq_pausing_index -t 10 -s hg19 -i pausing -l labels.txt -k KAS-seq.txt &"
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 index' output files. REQUIRED."
threadsHelp="-t [threads]: please specify the number of threads used for the pausing or termination index calculation. Default: 1."
assemblyidHelp="-s [assembly id]: please specify the genome assembly id of KAS-seq data. -s [assembly id]. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
indexHelp="-i [index]: please specify the index type you want to calculate. e.g. pausing or termination. REQUIRED."
labelsHelp="-l [labels]: please input the text file containing the labels of (sp)KAS-seq that show in 'KAS-pipe2 index' output files. Default: basename of KAS-seq files.
Example:
WT_rep1
WT.rep2                        ---labels.txt"
KASseqHelp="-k [KAS-seq]: please input the text file containing KAS-seq bigWig files used for pausing or termination index calculation. REQUIRED.
Example:
KAS-seq_WT_rep1.bigWig
KAS-seq_WT_rep2.bigWig         ---KAS-seq.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 index' shell script is applied to calculate the pausing or termination index."

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

# if no parameters was provided, 'KAS-pipe2 index' will print the help.
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

# Required options.
if test -z $prefix ;then
   echo ""
   echo "Please input the prefix (basename) of 'KAS-pipe2 index' output files. -o [prefix]"
   echo ""
   exit -1
fi

if test -z $assemblyid ;then
   echo ""
   echo "Please specify the reference genome assembly id of (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -s [assembly id]"
   echo ""
   exit -1
fi

# supported assembly id.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
    echo ""
    echo "Error: unsupported assembly id: $assemblyid ; Please specify the reference genome assembly id of (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11."
    echo ""
    exit -1
fi

if test -z $index ;then
   echo ""
   echo "Please specify the index you want to calculate. e.g. pausing or termination. -i [index type] "
   echo ""
   exit -1
fi

# check unupported index types.
if [[ $index != "pausing" ]] && [[ $index != "termination" ]] ;then
   echo ""
   echo "Error: unsupported index types: $index ; Please specify the index types. e.g. pausing or termination."
   echo ""
   exit -1
fi    

if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing indexed (sp)KAS-seq bam files used for pausing or termination index calculation. -k [KAS-seq]"
   echo ""
   exit -1
fi

# setup the default parameters
number_of_samples=$( awk 'END {print NR}' $KASseq )

if test -n "$labels" ;then
number_of_labels=$( awk 'END {print NR}' $labels )	
   if [[ $number_of_labels != $number_of_samples ]] ;then
      echo ""
      echo "Error:the number of labels isn't consistent with the number of samples."
      echo ""
      exit -1
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
label_basename=$( basename ${sample_selected} .bigWig )
echo $label_basename >> ${prefix}.labels_basename.txt
done
labels="${prefix}.labels_basename.txt"
fi 

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# get the list of KAS-seq bigWig files or labels.
KASseq_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $KASseq)
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

   # calculate the average KAS density of all samples.
   echo "Calculate the average KAS density of all samples on ${assemblyid} Refseq promoter or genebody ..."
   echo ""
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/NF }' ${prefix}_on_${assemblyid}_Refseq.promoter.bed > ${prefix}_on_${assemblyid}_Refseq.promoter.average
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/NF }' ${prefix}_on_${assemblyid}_Refseq.genebody.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.average
   echo "done."
   echo ""

   echo "Filtering matrix of KAS density on ${assemblyid} Refseq promoter ..."
   paste ${prefix}_on_${assemblyid}_Refseq.promoter.average ${assemblyid}_Refseq.promoter.genenames ${prefix}_on_${assemblyid}_Refseq.promoter.bed | sort -k 1 -n -r | sort -u -k2,2 -k3,3 | awk '$1>=5 {print $0}' > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed
   
   paste ${prefix}_on_${assemblyid}_Refseq.genebody.average ${assemblyid}_Refseq.genebody.genenames ${prefix}_on_${assemblyid}_Refseq.genebody.bed | sort -k 1 -n -r | sort -u -k2,2 -k3,3 > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   echo "done."
   echo ""
   
   for ((i=1; i<=${number_of_samples}; i++))
   do
   echo "Calculate the pausing index on ${assemblyid} Refseq genes for the ${i}th sample."
   awk -v x=$i '{printf("%s\t%s\t%.2f\n",$2,$3,$(x+6))}' ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i}
   awk -v x=$i '{printf("%s\t%s\t%.2f\n",$2,$3,$(x+6))}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}

   awk 'NR==FNR{a[$1]=$3;}NR!=FNR{print $0"\t"a[$1]}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i} ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i} | awk '{ if(NF==3){printf("%.2f\n",1)} else{printf("%.2f\n",$3/($4+0.01))} }' > ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i}

   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i}
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}

   echo "done."
   echo ""
   done
   
   echo "Generating the final pausing index on ${assemblyid}_Refseq.promoter.bed for all samples ..."
   awk '{printf("%s\t%s\n",$2,$3)}' ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.genename
   paste ${prefix}_on_${assemblyid}_Refseq.promoter.filter.genename ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.* > ${prefix}_on_${assemblyid}_Refseq.gene.filter.without_header.txt

   echo -e "genename\tstrand" > ${prefix}.header1.txt
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels > ${prefix}.header2.txt
   paste ${prefix}.header1.txt ${prefix}.header2.txt > ${prefix}.header.txt

   cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_Refseq.gene.filter.without_header.txt > ${prefix}_on_${assemblyid}_Refseq_gene_pausing_index.txt
   echo "done."
   echo ""

   echo "Clean up the intermediate files."
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.tab
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.tab
   rm -f ${assemblyid}_Refseq.promoter.genenames
   rm -f ${assemblyid}_Refseq.genebody.genenames
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.average
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.average
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.*
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.filter.genename
   rm -f ${prefix}_on_${assemblyid}_Refseq.gene.filter.without_header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}.header.txt
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

   # calculate the average KAS density of all samples.
   echo "Calculate the average KAS density of all samples on ${assemblyid} Refseq terminator or genebody ..."
   echo ""
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/NF }' ${prefix}_on_${assemblyid}_Refseq.terminator.bed > ${prefix}_on_${assemblyid}_Refseq.terminator.average
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/NF }' ${prefix}_on_${assemblyid}_Refseq.genebody.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.average
   echo "done."
   echo ""

   echo "Filtering matrix of KAS density on ${assemblyid} Refseq genebody ..."
   paste ${prefix}_on_${assemblyid}_Refseq.terminator.average ${assemblyid}_Refseq.terminator.genenames ${prefix}_on_${assemblyid}_Refseq.terminator.bed | sort -k 1 -n -r | sort -u -k2,2 -k3,3 > ${prefix}_on_${assemblyid}_Refseq.terminator.filter.bed

   paste ${prefix}_on_${assemblyid}_Refseq.genebody.average ${assemblyid}_Refseq.genebody.genenames ${prefix}_on_${assemblyid}_Refseq.genebody.bed | sort -k 1 -n -r | sort -u -k2,2 -k3,3 | awk '$1>=2 {print $0}' > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   echo "done."
   echo ""

   for ((i=1; i<=${number_of_samples}; i++))
   do
   echo "Calculate the pausing index on ${assemblyid} Refseq genes for the ${i}th sample."
   awk -v x=$i '{printf("%s\t%s\t%.2f\n",$2,$3,$(x+6))}' ${prefix}_on_${assemblyid}_Refseq.terminator.filter.bed > ${prefix}_on_${assemblyid}_Refseq.terminator.filter.KAS-seq.${i}
   awk -v x=$i '{printf("%s\t%s\t%.2f\n",$2,$3,$(x+6))}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}

   awk 'NR==FNR{a[$1]=$3;}NR!=FNR{print $0"\t"a[$1]}' ${prefix}_on_${assemblyid}_Refseq.terminator.filter.KAS-seq.${i} ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i} | awk '{ if(NF==3){printf("%.2f\n",nan)} else{printf("%.2f\n",$4/($3+0.1))} }' > ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i}
 
   rm -f ${prefix}_on_${assemblyid}_Refseq.terminator.filter.KAS-seq.${i}
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}

   echo "done."
   echo ""
   done

   echo "Generating the final termination index on ${assemblyid}_Refseq.terminator.bed for all samples ..."
   awk '{printf("%s\t%s\n",$2,$3)}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.genename
   paste ${prefix}_on_${assemblyid}_Refseq.genebody.filter.genename ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.* > ${prefix}_on_${assemblyid}_Refseq.gene.filter.without_header.txt

   echo -e "genename\tstrand" > ${prefix}.header1.txt
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels > ${prefix}.header2.txt
   paste ${prefix}.header1.txt ${prefix}.header2.txt > ${prefix}.header.txt

   cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_Refseq.gene.filter.without_header.txt > ${prefix}_on_${assemblyid}_Refseq_gene_termination_index.txt
   echo "done."
   echo ""

   echo "Clean up the intermediate files."
   rm -f ${prefix}_on_${assemblyid}_Refseq.terminator.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.terminator.tab
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.npz
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.tab
   rm -f ${assemblyid}_Refseq.terminator.genenames
   rm -f ${assemblyid}_Refseq.genebody.genenames
   rm -f ${prefix}_on_${assemblyid}_Refseq.terminator.average
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.average
   rm -f ${prefix}_on_${assemblyid}_Refseq.terminator.filter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   rm -f ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.*
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.genename
   rm -f ${prefix}_on_${assemblyid}_Refseq.gene.filter.without_header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}.header.txt
   echo "done."
   echo ""

fi 

rm -f ${prefix}.labels_basename.txt

echo "'KAS-pipe2 index' run successfully."
