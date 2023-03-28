#!/bin/bash
# 'KAS-Analyzer SST_enhancer' was developed by Ruitu Lyu on 1-22-2022.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-Analyzer SST_enhancer [ -h/--help ] [ -o prefix ] [ -t threads ] [ -s assembly id ] [ -e enhancer ] [ -p peaks ] [ -k KAS-seq ] "
exampleHelp="Example: nohup KAS-Analyzer SST_enhancer -o KAS-seq_SST_enhancers -t 10 -s mm10 -e H3K27ac_enhancers.bed -p KAS-seq_peaks.bed -k KAS-seq.rep1.bam,KAS-seq.rep2.bam &"
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-Analyzer SST_enhancer' output files. DEFAULT: basename of enhancer file."
threadsHelp="-t [threads]: please specify the number of threads used for single-stranded transcribing enhancers (SST_enhancers) identification. DEFAULT: 1."
assemblyidHelp="-s [assembly id]: please specify the genome assembly id. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
enhancerHelp="-e [enhancer]: please specify the enhancer file used for single-stranded transcribing enhancers (SST_enhancers) identification. Enhancer file can be H3K27ac, P300 or Med1 ChIP-seq peaks file. REQUIRED."
peaksHelp="-p [peaks]: please specify the (sp)KAS-seq peaks file. REQUIRED."
KASseqHelp="-k [KAS-seq]: please specify the indexed bam file of KAS-seq data used for single-stranded transcribing enhancers (SST_enhancers) identification. e.g. KAS-seq.rep1.bam,KAS-seq.rep2.bam. REQUIRED."
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-Analyzer SST_enhancer' shell script is applied to identify single-stranded transcribing enhancers (SST_enhancers)."

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
    echo -e "$enhancerHelp"
    echo -e ""
    echo -e "$peaksHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer SST_enhancer' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
   printHelpAndExit
fi

while getopts 'ho:t:s:e:p:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        o) prefix=$OPTARG ;;
        t) threads=$OPTARG ;;
        s) assemblyid=$OPTARG ;;
        e) enhancer=$OPTARG ;;
	p) peaks=$OPTARG ;;
	k) KASseq=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $assemblyid ;then
   echo ""
   echo "Please specify the reference genome assembly id of enhancer file. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -s [assembly id]"
   echo ""
   exit -1
fi

# supported assembly id.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
    echo ""
    echo "Error: unsupported assembly id: $assemblyid ; Please specify the assembly id of enhancer file. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -s [assembly id]"
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

if test -z $enhancer ;then
    echo ""
    echo "Please specify the enhancer file used for single stranded transcribing enhancers (SST_enhancers) identification. REQUIRED. -e [enhancer]"
    echo ""
    exit -1
fi

if test -z $peaks ;then
    echo ""
    echo "Please specify the (sp)KAS-seq peaks file. REQUIRED. -p [peaks]"
    echo ""
    exit -1
fi

if test -z $KASseq ;then
    echo ""
    echo "Please specify the indexed bam file of KAS-seq data used for single-stranded transcribing enhancers (SST_enhancers) identification. REQUIRED. -k [KAS-seq]"
    echo ""
    exit -1
fi

# setup the default options.
if test -z $threads ;then
   threads=1
fi

if test -z $prefix ;then
   prefix=$( basename ${enhancer} .bed )
fi

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# get the number of samples and the list of indexed bam files.
echo $KASseq > ${prefix}.KASseq.txt
sed -i "s/\,/ /g" ${prefix}.KASseq.txt
KASseq_list=$(cat ${prefix}.KASseq.txt)
number_of_samples=$( awk 'END {print NF}' ${prefix}.KASseq.txt )

# normalized the indexed bam files into bigWig files with RPKM.
cat /dev/null > ${prefix}.KAS-seq.RPKM.bigWig.txt
for ((i=1; i<=${number_of_samples}; i++))
do
sample_selected=$( awk -v x=$i '{print $x}' ${prefix}.KASseq.txt )
echo "Generate the index of $sample_selected."
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

KASseq_bigWig_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' ${prefix}.KAS-seq.RPKM.bigWig.txt)

# get the labels list.
echo "Get the list of labels of KAS-seq data."
echo ""
cat /dev/null > ${prefix}.labels.txt
for ((i=1; i<=${number_of_samples}; i++))
do
sample_selected=$( awk -v x=$i '{print $x}' ${prefix}.KASseq.txt )
labels_basename=$( basename ${sample_selected} .bam )
echo ${labels_basename} >> ${prefix}.labels.txt
done

labels_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' ${prefix}.labels.txt)
echo "done."
echo ""

# filter the distal KAS-seq overlapped enhancer list.
echo "Filter the distal (sp)KAS-seq overlapped enhancers from ${enhancer}."
echo ""
intersectBed -a $enhancer -b ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}.promoter.bed -v | intersectBed -a - -b $peaks -wa -f 0.2 | sort -u | sortBed -i > ${enhancer}.KAS.distal.bed
echo "done."
echo ""

# generate the left and right shores of enhancer. 
echo "Generate the left and right shores of enhancer regions: ${enhancer}.KAS.distal.bed ."
echo ""
awk '{printf("%s\t%d\t%d\t%d\n",$1,$2,$3,$3-$2)}' ${enhancer}.KAS.distal.bed | awk '{ if($2-$4<0){printf("%s\t%d\t%d\t%s\n",$1,0,$2,"enhancer"FNR)} else{printf("%s\t%d\t%d\t%s\n",$1,$2-$4,$2,"enhancer"FNR)} }' > ${enhancer}.KAS.distal.left.bed
awk '{printf("%s\t%d\t%d\t%d\n",$1,$2,$3,$3-$2)}' ${enhancer}.KAS.distal.bed | awk '{printf("%s\t%d\t%d\t%s\n",$1,$3,$3+$4,"enhancer"FNR)}' | intersectBed -a - -b ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes.bed -wa -wb | awk '{ if($7-$3>=0){printf("%s\t%d\t%d\t%s\n",$1,$2,$3,$4)} else{printf("%s\t%d\t%d\t%s\n",$1,$2,$7,$4)} }' > ${enhancer}.KAS.distal.right.bed
echo "done."
echo ""

echo "Calculating KAS-seq density on distal enhancers and enhancer left or right shores ..."
echo ""
multiBigwigSummary BED-file --bwfiles $KASseq_bigWig_list --BED ${enhancer}.KAS.distal.bed --labels $labels_list -p $threads -out ${prefix}_on_${enhancer}.KAS.distal.npz --outRawCounts ${prefix}_on_${enhancer}.KAS.distal.tab

multiBigwigSummary BED-file --bwfiles $KASseq_bigWig_list --BED ${enhancer}.KAS.distal.left.bed --labels $labels_list -p $threads -out ${prefix}_on_${enhancer}.KAS.distal.left.npz --outRawCounts ${prefix}_on_${enhancer}.KAS.distal.left.tab

multiBigwigSummary BED-file --bwfiles $KASseq_bigWig_list --BED ${enhancer}.KAS.distal.right.bed --labels $labels_list -p $threads -out ${prefix}_on_${enhancer}.KAS.distal.right.npz --outRawCounts ${prefix}_on_${enhancer}.KAS.distal.right.tab

sed "s/nan/0/g" ${prefix}_on_${enhancer}.KAS.distal.tab | sed "1d" | sortBed -i > ${prefix}_on_${enhancer}.KAS.distal.bed
sed "s/nan/0/g" ${prefix}_on_${enhancer}.KAS.distal.left.tab | sed "1d" | sortBed -i > ${prefix}_on_${enhancer}.KAS.distal.left.bed
sed "s/nan/0/g" ${prefix}_on_${enhancer}.KAS.distal.right.tab | sed "1d" | sortBed -i > ${prefix}_on_${enhancer}.KAS.distal.right.bed

awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/(NF-3) }' ${prefix}_on_${enhancer}.KAS.distal.bed > ${prefix}_on_${enhancer}.KAS.distal.average
awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/(NF-3) }' ${prefix}_on_${enhancer}.KAS.distal.left.bed > ${prefix}_on_${enhancer}.KAS.distal.left.average
awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/(NF-3) }' ${prefix}_on_${enhancer}.KAS.distal.right.bed > ${prefix}_on_${enhancer}.KAS.distal.right.average

awk '{printf("%s\t%d\t%d\n",$1,$2,$3)}' ${prefix}_on_${enhancer}.KAS.distal.bed > ${prefix}_on_${enhancer}.KAS.distal.3bed

echo "done."
echo ""

echo "Identify the single-stranded transcribing enhancers ..."
echo ""
paste ${prefix}_on_${enhancer}.KAS.distal.3bed ${prefix}_on_${enhancer}.KAS.distal.average ${prefix}_on_${enhancer}.KAS.distal.left.average ${prefix}_on_${enhancer}.KAS.distal.right.average | awk '$4*2/($5+$6+0.1)>=1.5 {print $0}' | awk '{printf("%s\t%d\t%d\t%s\t%.2f\t%.2f\n",$1,$2,$3,"ss_enhancer"FNR,$4,$4*2/($5+$6+0.01))}' > ${prefix}_SST_enhancers.txt
echo "done."
echo ""

echo "Clean up the intermediate files."
rm -f ${prefix}.KASseq.txt
rm -f ${prefix}.labels.txt
rm -f ${prefix}.KAS-seq.RPKM.bigWig.txt
rm -f $KASseq_bigWig_list
rm -f ${enhancer}.KAS.distal.bed
rm -f ${enhancer}.KAS.distal.left.bed
rm -f ${enhancer}.KAS.distal.right.bed
rm -f ${prefix}_on_${enhancer}.KAS.distal.npz
rm -f ${prefix}_on_${enhancer}.KAS.distal.tab
rm -f ${prefix}_on_${enhancer}.KAS.distal.bed
rm -f ${prefix}_on_${enhancer}.KAS.distal.left.npz
rm -f ${prefix}_on_${enhancer}.KAS.distal.left.tab
rm -f ${prefix}_on_${enhancer}.KAS.distal.left.bed
rm -f ${prefix}_on_${enhancer}.KAS.distal.right.npz
rm -f ${prefix}_on_${enhancer}.KAS.distal.right.tab
rm -f ${prefix}_on_${enhancer}.KAS.distal.right.bed
rm -f ${prefix}_on_${enhancer}.KAS.distal.3bed
rm -f ${prefix}_on_${enhancer}.KAS.distal.average
rm -f ${prefix}_on_${enhancer}.KAS.distal.left.average
rm -f ${prefix}_on_${enhancer}.KAS.distal.right.average
echo "done."
echo ""

echo "'KAS-Analyzer SST_enhancer' run successfully."
