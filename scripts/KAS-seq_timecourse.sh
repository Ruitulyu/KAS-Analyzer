#!/bin/bash
# 'KAS-pipe2 TC' was developed by Ruitu Lyu on 1-19-2022.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 TC [ -h/--help ] [ -t threads ] [ -o prefix ] [ -g assembly id ] [ -r regions ] [ -s bin size ] [ -b ] [ -p peaks ] [ -d differential analysis ] [ -a annotation ] [ -l labels ] [ -k KAS-seq ] "
exampleHelp="Example: 
Case-only:
nohup KAS-pipe2 TC -o KAS-seq_timecourse -t 10 -g mm10 -r bin -d case_only -a Case_only.annotation.txt -l labels.txt -k KAS-seq_data.txt &

Case-control:
nohup KAS-pipe2 TC -o KAS-seq_timecourse -t 10 -g mm10 -r bin -d case_control -a Case_control.annotation.txt -l labels.txt -k KAS-seq_data.txt &"
threadsHelp="-t [threads]: please specify the number of threads used for calculating KAS index. DEFAULT: 1."
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 KASindex' output files. REQUIRED."
assemblyidHelp="-g [assembly id]: please specify the genome assembly id of KAS-seq data. -g [assembly id]. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
regionsHelp="-r [regions]: please specify the types of genomic regions. e.g. promoter, genebody, bin, gene or peak. REQUIRED."
binsizeHelp="-s [bin size]: please specify the size of bins if -r is set to 'bin'. DEFAULT:1000."
batchHelp="-b: please specify if KAS-seq data have batch effects. DEFAULT: off."
peaksHelp="-p [peaks]: please specify the custom regions file for perform differential time course KAS-seq analysis. OPTIONAL."
diffanalysisHelp="-d [diff analysis type]: please specify the types of differential time course KAS-seq analysis. e.g. case_only or case_control. REQUIRED."
annotationHelp="-a [annotation]: please specify the annotation file. REQUIRED. Example:
1) Case-only:                                                or 2) Case-Control:
	Sample	Condition	Time	Batch				Sample  Condition       Time    Batch          or	Batch
rep1.0min	rep1.0min	case	1	B_NULL		WT.rep1.0min	WT.rep1.0min	control	1	B_NULL		B1
rep2.0min	rep2.0min	case	1	B_NULL		WT.rep2.0min	WT.rep2.0min	control	1	B_NULL		B2	
rep3.0min	rep3.0min	case	1	B_NULL		WT.rep3.0min	WT.rep3.0min	control	1	B_NULL		B3
rep1.15min	rep1.15min	case	2	B_NULL		WT.rep1.15min	WT.rep1.15min	control	2	B_NULL		B1
rep2.15min	rep2.15min	case	2	B_NULL		WT.rep2.15min	WT.rep2.15min	control	2	B_NULL		B2
rep3.15min	rep3.15min	case	2	B_NULL		WT.rep3.15min	WT.rep3.15min	control	2	B_NULL		B3
rep1.30min	rep1.30min	case	3	B_NULL		WT.rep1.30min	WT.rep1.30min	control	3	B_NULL		B1
rep2.30min	rep2.30min	case	3	B_NULL		WT.rep2.30min	WT.rep2.30min	control	3	B_NULL		B2
rep3.30min	rep3.30min	case	3	B_NULL		WT.rep3.30min	WT.rep3.30min	control	3	B_NULL		B3
rep1.60min	rep1.60min	case	4	B_NULL		KO.rep1.0min	KO.rep1.0min	case	1	B_NULL		B1
rep2.60min	rep2.60min	case	4	B_NULL		KO.rep2.0min	KO.rep2.0min	case	1	B_NULL		B2
rep3.60min	rep3.60min	case	4	B_NULL		KO.rep3.0min	KO.rep3.0min	case	1	B_NULL		B3
								KO.rep1.15min	KO.rep1.15min	case	2	B_NULL		B1		
								KO.rep2.15min	KO.rep2.15min	case	2	B_NULL		B2
                                                                KO.rep3.15min	KO.rep3.15min	case	2	B_NULL		B3
								KO.rep1.30min	KO.rep1.30min	case	3	B_NULL		B1
								KO.rep2.30min	KO.rep2.30min	case	3	B_NULL		B2
								KO.rep3.30min	KO.rep3.30min	case	3	B_NULL		B3          ---annotation.txt"
labelsHelp="-l [labels]: please input the text file containing the labels of KAS-seq data for differential time course (TC) analysis. REQUIRED. Example:
1) Case-only:                   or 2) Case-Control:
rep1.0min                          WT.rep1.0min
rep2.0min                          WT.rep2.0min
rep3.0min                          WT.rep3.0min
rep1.15min                         WT.rep1.15min
rep2.15min                         WT.rep2.15min
rep3.15min                         WT.rep3.15min
rep1.30min                         WT.rep1.30min
rep2.30min                         WT.rep2.30min
rep3.30min                         WT.rep3.30min
rep1.60min                         KO.rep1.0min
rep2.60min                         KO.rep2.0min
rep3.60min                         KO.rep3.0min
                                   KO.rep1.15min
                                   KO.rep2.15min
                                   KO.rep3.15min
                                   KO.rep1.30min
				   KO.rep2.30min
				   KO.rep3.30min                                ---labels.txt"
KASseqHelp="-k [KAS-seq]: please input the text file containing indexed bam files of KAS-seq data for differential time course (TC) KAS-seq analysis. REQUIRED. Example:
1) Case-only:                   or 2) Case-Control:
treat.rep1.0min.bam                WT.treat.rep1.0min.bam
treat.rep2.0min.bam                WT.treat.rep2.0min.bam
treat.rep3.0min.bam                WT.treat.rep3.0min.bam
treat.rep1.15min.bam               WT.treat.rep1.15min.bam
treat.rep2.15min.bam               WT.treat.rep2.15min.bam
treat.rep3.15min.bam               WT.treat.rep3.15min.bam
treat.rep1.30min.bam               WT.treat.rep1.30min.bam
treat.rep2.30min.bam               WT.treat.rep2.30min.bam
treat.rep3.30min.bam               WT.treat.rep3.30min.bam
treat.rep1.60min.bam               KO.treat.rep1.0min.bam
treat.rep2.60min.bam               KO.treat.rep2.0min.bam
treat.rep3.60min.bam               KO.treat.rep3.0min.bam           
                                   KO.treat.rep1.15min.bam
                                   KO.treat.rep2.15min.bam
                                   KO.treat.rep3.15min.bam
                                   KO.treat.rep1.30min.bam
                                   KO.treat.rep2.30min.bam                                   
                                   KO.treat.rep3.30min.bam                      ---KAS-seq_data.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 TC' shell script is applied to performed differential analysis for 'case only' or 'case-control' time course(TC) KAS-seq on promoter, genebody, bin, genes or custom regions."

# print help.
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
    echo -e "$binsizeHelp"
    echo -e ""
    echo -e "$batchHelp"
    echo -e ""
    echo -e "$peaksHelp"
    echo -e ""
    echo -e "$diffanalysisHelp"
    echo -e ""
    echo -e "$annotationHelp"
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

while getopts 'ht:o:g:r:s:bp:d:a:l:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        t) threads=$OPTARG ;;
        o) prefix=$OPTARG ;;
        g) assemblyid=$OPTARG ;;
        r) regions=$OPTARG ;;
        s) binsize=$OPTARG ;;
        b) batch="on" ;;
        p) peaks=$OPTARG ;;
        d) difftypes=$OPTARG ;;
        a) annotation=$OPTARG ;;
	l) labels=$OPTARG ;;
	k) KASseq=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $assemblyid ;then
   echo ""
   echo "Please input the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. -g [assembly id]"
   echo ""
   exit 1
fi

# supported assembly id.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39
" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce
11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
    echo ""
    echo "Error: unsupported assembly id: $assemblyid ; Please specify the reference genome assembly id of (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C
.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11."
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

if test -z $prefix ;then
   echo ""
   echo "Please provide the prefix (basename) of 'KAS-pipe2 TC' output files. -o [prefix]"
   echo ""
   exit 1
fi

if test -z $regions ;then
   echo ""
   echo "Please specify the the types of genomic regions for differential TC KAS-seq analysis. e.g. promoter, genebody, bin, gene or peak. -r [regions] "
   echo ""
   exit 1
fi

# supported regions types.
if [[ $regions != "bin" ]] && [[ $regions != "promoter" ]] && [[ $regions != "genebody" ]] && [[ $regions != "gene" ]] && [[ $regions != "peak" ]] ;then
   echo ""
   echo "Error: unsupported type of genomic regions: ${assemblyid}; please specify the types of genomic regions for differential time course (TC) KAS-seq analysis. e.g. promoter, genebody, gene or peak."
   echo ""
   exit 1
fi

if [[ $regions == "peak" ]] && test -z $peaks ;then
   echo ""      
   echo "Please specify the bed file of custom regions for differential time course (TC) KAS-seq analysis. -p [peaks]."
   echo ""
   exit 1
fi

if test -z $difftypes ;then
   echo ""
   echo "Please specify the types of differential time course (TC) KAS-seq analysis. e.g. case_only or case_control. -d [analysis type]."
   echo ""
   exit 1
fi

# supported types of differential TC KAS-seq analysis.
if [[ $difftypes != "case_only" ]] && [[ $difftypes != "case_control" ]] ;then
   echo ""
   echo "Error: unsupported differential TC KAS-seq analysis types: $difftypes; Please specify the types of differential TC KAS-seq analysis. e.g. case_only or case_control."
   echo ""
   exit 1
fi   


if test -z $annotation ;then
   echo ""
   echo "Please specify the annotation file. -a [annotation]"
   echo ""
   exit 1
fi 

if test -z $labels ;then
   echo ""
   echo "Please input the text file containing the labels of KAS-seq data for differential time course (TC) KAS-seq analysis. -l [labels]."
   echo ""
   exit 1
fi

if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing indexed bam files of KAS-seq data for differential time course (TC) KAS-seq analysis. -k [KAS-seq]"
   echo ""
   exit 1
fi

# setup the default options.
if test -z $threads ;then
   threads=1
fi

if [[ $regions == "bin" ]] && test -z $binsize ;then
   binsize=1000
fi

if test -z $batch ;then
   batch="off"
fi

# get the number of samples.
number_of_samples=$( awk 'END {print NR}' $KASseq )
number_of_labels=$( awk 'END {print NR}' $labels )

if [[ $number_of_labels != $number_of_samples ]] ;then
   echo ""
   echo "Error: the number of labels isn't consistent with the number of samples."
   echo ""
   exit 1
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
labels_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $labels)

if [[ $regions == "bin" ]] ;then

   echo "Generating the KAS-seq RPKM values matrix on ${assemblyid} ${binsize}bp bins..."
   echo ""
   multiBigwigSummary bins -b $KASseq_list --labels $labels_list --binSize $binsize -p $threads --blackListFileName ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -out ${prefix}_on_${binsize}_bins.npz --outRawCounts ${prefix}_on_${binsize}_bins.tab > /dev/null 2>&1
   sed "s/nan/0/g" ${prefix}_on_${binsize}_bins.tab | sed "1d" > ${prefix}_on_${binsize}_bins.bed
   echo "done."
   echo ""   
   	
   #calculate the average value of very single row in the table. 
   echo "Calculating the average KAS-seq RPKM values on ${assemblyid} ${binsize}bp bins ..."
   echo ""
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print  a[j]/NF }' ${prefix}_on_${binsize}_bins.bed > ${prefix}_on_${binsize}_bins.average
   echo "done."
   echo ""

   #filter the bins with averaged KAS-seq expression over 2.
   echo "Filter KAS-seq KAS-seq RPKM values matrix on ${binsize}bp bins ..."
   echo ""
   paste ${prefix}_on_${binsize}_bins.average ${prefix}_on_${binsize}_bins.bed | awk '$1>=2 {print $0}' > ${prefix}_on_${binsize}_bins.filter.bed
   cut -f1,2,3,4 --complement ${prefix}_on_${binsize}_bins.filter.bed | awk '{for(i=1;i<=NF;i++){printf "%d\t", $i*10}; printf "\n"}' > ${prefix}_on_${assemblyid}_${regions}.matrix
   awk '{printf("%s\n","bin"FNR"-"$2"-"$3"-"$4)}' ${prefix}_on_${binsize}_bins.filter.bed > ${prefix}_on_${assemblyid}_${regions}.rowname
   echo "done."
   echo ""

   echo "Generating the KAS-seq RPKM values interger matrix with binname and header ..."
   echo ""
   paste ${prefix}_on_${assemblyid}_${regions}.rowname ${prefix}_on_${assemblyid}_${regions}.matrix > ${prefix}_on_${assemblyid}_${regions}.without_header.txt

   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' $labels > ${prefix}.header.txt
   cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_${regions}.without_header.txt > ${prefix}_on_${assemblyid}_${regions}.txt
   echo "done."
   echo ""

   # clean up the intermediate files.
   echo "clean up the intermediate files"
   echo ""
   rm -f ${prefix}_on_${binsize}_bins.tab
   rm -f ${prefix}_on_${binsize}_bins.npz
   rm -f ${prefix}_on_${binsize}_bins.bed
   rm -f ${prefix}_on_${binsize}_bins.average
   rm -f ${prefix}_on_${binsize}_bins.filter.bed
   rm -f ${prefix}_on_${assemblyid}_${regions}.matrix
   rm -f ${prefix}_on_${assemblyid}_${regions}.rowname
   rm -f ${prefix}_on_${assemblyid}_${regions}.without_header.txt
   rm -f ${prefix}.header.txt
   echo "done."
   echo ""

elif [[ $regions == "peak" ]] ;then
  
   echo "Generating the KAS-seq matrix on ${peaks}..."
   echo ""	
   peaks_basename=$( basename ${peaks} .bed )
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED $peaks --labels $labels_list -p $threads -out ${prefix}_on_${peaks_basename}.npz --outRawCounts ${prefix}_on_${peaks_basename}.tab > /dev/null 2>&1
   sed "s/nan/0/g" ${prefix}_on_${peaks_basename}.tab | sed "1d" > ${prefix}_on_${peaks_basename}.bed
   echo "done."
   echo ""

   echo "Generating the KAS-seq RPKM values interger matrix with rowname and header ..."
   echo ""
   # generate the KAS-seq signals integer matrix multiply by 10.
   cut -f1,2,3 --complement ${prefix}_on_${peaks_basename}.bed | awk '{for(i=1;i<=NF;i++){printf "%d\t", $i*10}; printf "\n"}' > ${prefix}_on_${peaks_basename}.matrix
   awk '{printf("%s\n","peak"NR"-"$1"-"$2"-"$3)}' ${prefix}_on_${peaks_basename}.bed > ${prefix}_on_${peaks_basename}.rowname
   # generating the KAS-seq signals matrix without header.
   paste ${prefix}_on_${peaks_basename}.rowname ${prefix}_on_${peaks_basename}.matrix > ${prefix}_on_${assemblyid}_${regions}.without_header.txt
   
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' $labels > ${prefix}.header.txt
   cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_${regions}.without_header.txt > ${prefix}_on_${assemblyid}_${regions}.txt
   echo "done."
   echo ""

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
   echo "done."
   echo ""

elif [[ $regions == "promoter" ]] || [[ $regions == "genebody" ]] ;then

   echo "Calculating KAS-seq RPKM values of ${assemblyid} Refseq ${regions} ..."
   echo ""
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.${regions}.bed --labels $labels_list -p $threads -out ${prefix}_on_${assemblyid}_Refseq.${regions}.npz --outRawCounts ${prefix}_on_${assemblyid}_Refseq.${regions}.tab > /dev/null 2>&1
   sed "s/nan/0/g" ${prefix}_on_${assemblyid}_Refseq.${regions}.tab | sed "1d" | sortBed -i > ${prefix}_on_${assemblyid}_Refseq.${regions}.bed
   sortBed -i ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.${regions}.bed | awk '{printf("%s\n",$1"--"$2"--"$3"--"$4"--"$6)}' > ${assemblyid}_Refseq.${regions}.genenames
   echo "done."
   echo ""

   # calculate the average KAS RPKM values of every single row in the table.
   echo "Calculating the averaged KAS-seq RPKM values of ${assemblyid}_Refseq.${regions}.bed ..."
   echo ""
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/NF }' ${prefix}_on_${assemblyid}_Refseq.${regions}.bed > ${prefix}_on_${assemblyid}_Refseq.${regions}.average
   echo "done."
   echo ""

   # filter the ${regions} with KAS expression over 2.
   echo "Filtering matrix of KAS-seq RPKM values on ${assemblyid}_Refseq.${regions}.bed ..."
   echo ""
   paste ${prefix}_on_${assemblyid}_Refseq.${regions}.average ${assemblyid}_Refseq.${regions}.genenames ${prefix}_on_${assemblyid}_Refseq.${regions}.bed | awk '$1>=2 {print $0}' > ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed
   echo "done."
   echo ""

   echo "Generating the KAS-seq RPKM values interger matrix with genenames ..."
   echo ""
   cut -f1,2,3,4,5 --complement ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed | awk '{for(i=1;i<=NF;i++){printf "%d\t", $i*10}; printf "\n"}' > ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.matrix
   awk '{printf("%s\n",$2)}' ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.bed > ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.genenames

   paste ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.genenames ${prefix}_on_${assemblyid}_Refseq.${regions}.filter.matrix > ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt
   echo "done."
   echo ""

   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' $labels > ${prefix}.header.txt
   cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_Refseq.${regions}.without_header.txt > ${prefix}_on_${assemblyid}_${regions}.txt

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
   rm -f ${prefix}.header.txt
   echo "done."
   echo ""

elif [[ $regions == "gene" ]] ;then
   
   echo "Generating the KAS-seq RPKM values matrix on ${assemblyid} Refseq ${regions}..." 
   echo ""
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.promoter.bed --labels $labels_list -p $threads -out ${prefix}_on_${assemblyid}_Refseq.promoter.npz --outRawCounts ${prefix}_on_${assemblyid}_Refseq.promoter.tab > /dev/null 2>&1

   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.genebody.bed --labels $labels_list -p $threads -out ${prefix}_on_${assemblyid}_Refseq.genebody.npz --outRawCounts ${prefix}_on_${assemblyid}_Refseq.genebody.tab > /dev/null 2>&1

   sed "s/nan/0/g" ${prefix}_on_${assemblyid}_Refseq.promoter.tab | sed "1d" | sortBed -i > ${prefix}_on_${assemblyid}_Refseq.promoter.bed
   sortBed -i ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.promoter.bed | awk '{printf("%s\t%s\n",$4,$6)}' > ${assemblyid}_Refseq.promoter.genenames

   sed "s/nan/0/g" ${prefix}_on_${assemblyid}_Refseq.genebody.tab | sed "1d" | sortBed -i > ${prefix}_on_${assemblyid}_Refseq.genebody.bed
   sortBed -i ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.genebody.bed | awk '{printf("%s\t%s\n",$4,$6)}' > ${assemblyid}_Refseq.genebody.genenames
   echo "done."
   echo ""	

   # calculate the averaged KAS-seq RPKM values of every single row in the table.
   echo "Calculate the average KAS-seq RPKM values on ${assemblyid} Refseq promoter and genebody ..."
   echo ""
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/NF }' ${prefix}_on_${assemblyid}_Refseq.promoter.bed > ${prefix}_on_${assemblyid}_Refseq.promoter.average
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/NF }' ${prefix}_on_${assemblyid}_Refseq.genebody.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.average
   echo "done."
   echo ""

   echo "Filter matrix of KAS-seq RPKM values on ${assemblyid} Refseq promoter and genebody ..."
   echo ""
   paste ${prefix}_on_${assemblyid}_Refseq.promoter.average ${assemblyid}_Refseq.promoter.genenames ${prefix}_on_${assemblyid}_Refseq.promoter.bed | sort -k 1 -n -r | sort -u -k2,2 -k3,3 | awk '$1>=2 {print $0}' > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed
   paste ${prefix}_on_${assemblyid}_Refseq.genebody.average ${assemblyid}_Refseq.genebody.genenames ${prefix}_on_${assemblyid}_Refseq.genebody.bed | sort -k 1 -n -r | sort -u -k2,2 -k3,3 > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed
   echo "done."
   echo ""

   echo "Generating the final KAS-seq RPKM values matrix on ${assemblyid}_Refseq.gene.bed for all samples ..."
   echo ""
   awk '{printf("%s\n",$2"--"$3)}' ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed > ${prefix}_on_${assemblyid}_Refseq.genes.matrix.0

   for ((i=1; i<=${number_of_samples}; i++))
   do
   echo "Combine KAS-seq RPKM values on promoter and genebody, and calculate the KAS-seq RPKM values on ${assemblyid} Refseq genes for the ${i}th sample."
   echo ""
   awk -v x=$i '{printf("%s\t%d\n",$2"-"$3,$(x+6)*10)}' ${prefix}_on_${assemblyid}_Refseq.promoter.filter.bed > ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i}
   awk -v x=$i '{printf("%s\t%d\n",$2"-"$3,$(x+6)*10)}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.bed > ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}

   awk 'NR==FNR{a[$1]=$2;}NR!=FNR{print $0"\t"a[$1]}' ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i} ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i} | awk '{ if(NF==2){printf("%d\n",$2)} else{printf("%d\n",($2+$3)/2)} }' > ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i}
  
   paste ${prefix}_on_${assemblyid}_Refseq.genes.matrix.$((i-1)) ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i} > ${prefix}_on_${assemblyid}_Refseq.genes.matrix.${i}
   
   rm -f ${prefix}_on_${assemblyid}_Refseq.promoter.filter.KAS-seq.${i}
   rm -f ${prefix}_on_${assemblyid}_Refseq.genebody.filter.KAS-seq.${i}
   rm -f ${prefix}_on_${assemblyid}_Refseq.gene.filter.KAS-seq.${i}

   echo "done."
   echo ""
   done

   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' $labels > ${prefix}.header.txt
   cat ${prefix}.header.txt ${prefix}_on_${assemblyid}_Refseq.genes.matrix.${number_of_samples} > ${prefix}_on_${assemblyid}_${regions}.txt
   echo "Final KAS-seq RPKM values matrix generation done."
   echo ""

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
   rm -f ${prefix}_on_${assemblyid}_Refseq.genes.matrix.*
   rm -f ${prefix}.header.txt
   echo "done."
   echo ""

fi

if [[ $difftypes == "case_only" ]] ;then

   if [[ $batch == "off" ]] ;then       
      echo "Perform 'case_only' differential time course KAS-seq analysis without batch effect normalization ..."
      echo ""
      echo ${prefix}_on_${assemblyid}_${regions}.txt
      echo $annotation
      exit
      Rscript --vanilla ${SH_SCRIPT_DIR}/../R/ImpulseDE2_TC_diff_KAS-seq.case_only.R ${prefix}_on_${assemblyid}_${regions}.txt $annotation
      echo "done."
      echo ""
   elif [[ $batch == "on" ]] ;then
      echo "Perform 'case_only' differential time course KAS-seq analysis with batch effect normalization ..."	   
      echo ""
      Rscript --vanilla ${SH_SCRIPT_DIR}/../R/ImpulseDE2_TC_diff_KAS-seq.case_only.batch-effect.R ${prefix}_on_${assemblyid}_${regions}.txt $annotation
      echo "done."
      echo ""
   fi   
elif [[ $difftypes == "case_control" ]]	;then

   if [[ $batch == "off" ]] ;then      
      echo "Perform 'case_control' differential time course KAS-seq analysis without batch effect normalization ... "
      echo ""
      Rscript --vanilla ${SH_SCRIPT_DIR}/../R/ImpulseDE2_TC_diff_KAS-seq.case_control.R ${prefix}_on_${assemblyid}_${regions}.txt $annotation
      echo "done."
      echo ""
   elif [[ $batch == "on" ]] ;then
      echo "Perform 'case_control' differential time course KAS-seq analysis with batch effect normalization ..."
      echo ""
      Rscript --vanilla ${SH_SCRIPT_DIR}/../R/ImpulseDE2_TC_diff_KAS-seq.case_control.batch-effect.R ${prefix}_on_${assemblyid}_${regions}.txt $annotation
      echo "done."
      echo ""
   fi	
fi	


rm -f $KASseq_list
rm -f ${prefix}.KAS-seq.RPKM.bigWig.txt

echo "'KAS-pipe2 TC' run successfully!"
