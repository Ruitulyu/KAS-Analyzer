#!/bin/bash
# 'KAS-pipe2 KAS-seq' was developped by Ruitu Lyu on 12-08-2021.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-pipe2 KAS-seq [ -h ] [ -a aligner ] [ -t threads ] [ -i index path ] [ -e extend length ] [ -o prefix ] [ -s assembly id ] [ -1 read1 ] [ -2 read2 ]"
exampleHelp="Example:
       Single-end:
       nohup KAS-pipe2 KAS-seq -a bowtie2 -t 10 -i /absolute path/hg19_Bowtie2Index/hg19 -e 150 -o KAS-seq -s hg19 -1 KAS-seq.trim.fastq.gz &
       Paired-end:
       nohup KAS-pipe2 KAS-seq -a bowtie2 -t 10 -i /absolute path/hg19_Bowtie2Index/hg19 -o KAS-seq -s hg19 -1 KAS-seq.trim.R1.fastq.gz -2 KAS-seq.trim.R2.fastq.gz &
       Note: Bowtie2 Index example: /absolute path/Bowtie2Index/hg19; bwa Index example: /absolute path/BWAIndex/hg19.fa"
alignerHelp="-a [aligner]: please specify the aligner (bowtie2 or bwa) you want to use to map KAS-seq data. Default: bowtie2."
threadsHelp="-t [threads]: please input the number of threads used for KAS-seq data mapping. Default: 1."
indexpathHelp="-i [index path]: please input the absolute path of reference genome index for aligner. REQUIRED."
extendlengthHelp="-e [extendlengthHelp]: please input the extend length for single-end KAS-seq data. Default: 150."
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 KAS-seq' output files. REQUIRED."
assemblyidHelp="-s [assembly id]: please specify the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. REQUIRED."
read1Help="-1 [read1]: please input trimmed single-end KAS-seq fastq file or read 1 of paired-end KAS-seq raw fastq files; Compressed .fastq.gz is accepted. REQUIRED."
read2Help="-2 [read2]: please input trimmed read2 of paired-end KAS-seq raw fastq files. Compressed .fastq.gz is accepted."
helpHelp="-h: print this help and exit.
Note: The 'KAS-pipe2 KAS-seq' shell script mainly invoke the specific aligner (bowtie, bowtie2 or bwa) for KAS-seq data mapping, please refer to their official websites for more information."

# print help.
printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "" 
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$alignerHelp"
    echo -e ""
    echo -e "$threadsHelp"
    echo -e ""
    echo -e "$indexpathHelp"
    echo -e ""
    echo -e "$extendlengthHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$read1Help"
    echo -e ""
    echo -e "$read2Help"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 KAS-seq' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

while getopts 'ha:t:i:e:o:s:1:2:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        a) aligner=$OPTARG ;;
        t) threads=$OPTARG ;;
        i) indexpath=$OPTARG ;;
	e) extendlength=$OPTARG ;;
        o) prefix=$OPTARG ;;
        s) assemblyid=$OPTARG ;;
        1) read1=$OPTARG ;;
        2) read2=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done


# Required options. 
if test -z $indexpath ;then
   echo ""
   echo "Please provide the absolute path of reference genome index for aligner. -i [index path]"
   echo ""
   exit 1
fi

if test -z $prefix ;then
   echo ""
   echo "Please provide the prefix (basename) of 'KAS-pipe2 KAS-seq' output files. -o [prefix]"
   echo ""
   exit 1
fi

if test -z $assemblyid ;then
   echo ""	
   echo "Please input the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. -s [assembly id]"
   echo ""
   exit 1
fi

if test -z $read1 ;then
   echo ""
   echo "Please input the single-end KAS-seq fastq file or read1 of paired-end KAS-seq data; .fastq.gz is accepted. -1 [read1]"
   echo ""
   exit 1
fi

# setup the default options.

if test -z $threads ;then
   threads=1
fi

if test -z $aligner ;then
   aligner="bowtie2"
elif [[ $aligner != "bowtie2" ]] && [[ $aligner != "bwa" ]] ;then
   echo ""	
   echo "Error: unsupported aligner: ${aligner}. Please input the aligner: bowtie2 or bwa."
   echo ""
   exit 1
fi


if test -z $extendlength ;then
   extendlength=150
fi

# test if aligners was installed.
if ! type $aligner > /dev/null 2>&1 ;then
   echo "$aligner was not installed or not export to the \$PATH'"
   echo ""
   echo "Install ${aligner} with 'conda install -c bioconda ${aligner}' or refer the official website of '${aligner}'."
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

if ! type bedtools > /dev/null 2>&1 ;then
   echo "bedtools was not installed or not export to the \$PATH'"
   echo ""
   echo "Install bedtools with 'conda install -c bioconda bedtools' or refer the official website of 'bedtools'."
   echo ""
   exit 1
fi


# Determine single or paired KAS-seq data.
if test -z $read2 ;then
   paired_or_single_end="single"
else
   paired_or_single_end="paired"
fi

# get the path of 'KAS-pipe2 KAS-seq' shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)


echo ""
echo "Welcome to analyze KAS-seq data with 'KAS-pipe2'..."
echo ""

# Map single-end KAS-seq data.
if [[ $paired_or_single_end == "single" ]]; then
   echo "Map single-end KAS-seq data."
   echo ""
   
   # map single-end KAS-seq data with bowtie2.
   if [[ $aligner == "bowtie2" ]] ;then      
      # create $prefix folder and move read1 KAS-seq data into $prefix folder
      mkdir -p $prefix
      cd $prefix
      mv ../$read1 ./

      # ouput the numebr of raw reads into mapping_summary txt file.
      # test if $read1 is .gz compressed file.
      if [[ ${read1##*.} == gz ]] ;then
         reads_num=$( zcat $read1 | wc -l | awk '{print $1/4}' )
      else
         reads_num=$( wc -l $read1 | awk '{print $1/4}' )
      fi

      echo "Number of KAS-seq reads: $reads_num" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
 
      echo "Map $read1 to $assemblyid reference genome with $aligner."
      echo ""
      # map $read1 to reference genome with bowtie2.
      bowtie2 -p $threads -x $indexpath $read1 -S ${prefix}.sam >> .${prefix}_SE_KAS-seq_mapping_statistics.txt 2>&1
      mapping_ratios=$( grep "overall alignment rate" .${prefix}_SE_KAS-seq_mapping_statistics.txt | awk '{print $1}' )
      mono_mapped_reads_num=$( grep "aligned exactly 1 time" .${prefix}_SE_KAS-seq_mapping_statistics.txt | awk '{print $1}' )
      multi_mapped_reads_num=$( grep "aligned >1 times" .${prefix}_SE_KAS-seq_mapping_statistics.txt | awk '{print $1}' )
      mapped_reads_num=$(( mono_mapped_reads_num + multi_mapped_reads_num ))
      echo "Number of mapped reads. ${aligner}: $mapped_reads_num" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      echo "Mapping ratios: $mapping_ratios" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      rm -f .${prefix}_SE_KAS-seq_mapping_statistics.txt
      echo "$aligner mapping done."
      echo ""
      
      # sort sam file with 'samtools sort' and transfer to bam file.
      echo "Sort ${prefix}.sam file with 'samtools sort' and transfer into ${prefix}_sorted.bam file."
      echo ""
      samtools sort -@ 3 ${prefix}.sam -o ${prefix}_sorted.bam
      samtools index ${prefix}_sorted.bam
      echo "'samtools sort' done."
      echo ""

      # remove duplicates with samtools rmdup.
      echo "" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      echo "Remove duplicates in ${prefix}_sorted.bam with 'samtools rmdup'"
      echo ""
      samtools rmdup -s ${prefix}_sorted.bam ${prefix}_rmdup.bam >> .${prefix}_SE_KAS-seq_deduplication_ratios.txt 2>&1
      samtools index ${prefix}_rmdup.bam
      unique_reads_num=$( samtools view -h ${prefix}_rmdup.bam | grep chr | grep ^@ -v | awk '{print $1}' | sort -u | wc -l )
      echo "Number of unique mapped reads: $unique_reads_num" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      duplication_ratios=$( awk '{print $6}' .${prefix}_SE_KAS-seq_deduplication_ratios.txt | awk '{printf("%s\n",$1*100"%")}' )
      echo "Duplication ratios. samtools: $duplication_ratios" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      rm -f .${prefix}_SE_KAS-seq_deduplication_ratios.txt
      echo "'samtools rmdup' done."
      echo ""
      
      # transfer bam into bed.
      echo "Transfer ${prefix}_rmdup.bam into ${prefix}.bed with bamToBed."
      echo ""
      bamToBed -i ${prefix}_rmdup.bam > ${prefix}.bed
      intersectBed -a ${prefix}.bed -b ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -v > ${prefix}.filter.bed 
      mv ${prefix}.filter.bed ${prefix}.bed
      echo "'bamToBed' done."
      echo ""
      
      # extend the deduplicated mapped reads to extendlength.  
      echo "Extend the deduplicated reads in ${prefix}.bed to ${extendlength}."
      echo ""
      awk '$3-'${extendlength}'>0 {if ($6~"+") printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$2+'${extendlength}',$4,$5,$6); else if ($6~"-") printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$3-'${extendlength}',$3,$4,$5,$6)}' ${prefix}.bed > ${prefix}.ext${extendlength}.bed
      echo "'extend' done."
      echo ""

      echo "Transfer ${prefix}.ext${extendlength}.bed into ${prefix}.ext${extendlength}.bg with genomeCoverageBed."
      echo ""
      genomeCoverageBed -bg -i ${prefix}.ext${extendlength}.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.ext${extendlength}.bg 
      echo "'genomeCoverageBed' done."
      echo ""
     
      # clean up the .sam, sorted .bam and unextended .bed files.
      rm -f ${prefix}.sam
      rm -f ${prefix}_rmdup.bam.bai
      rm -f ${prefix}_rmdup.bam
      rm -f ${prefix}.bed
      echo "Clean up intermediate files. done."
      echo ""

      # move deduplicated bam files and index into Bam_files.
      cd ..
      mkdir -p Bam_files
      cd Bam_files
      mv ../${prefix}/${prefix}_sorted.bam ./
      mv ../${prefix}/${prefix}_sorted.bam.bai ./
      cd ..

      # move bedGraph files into BedGraph files.
      mkdir -p BedGraph_files
      cd BedGraph_files
      mv ../${prefix}/${prefix}.ext${extendlength}.bg ./
      cd ..

      # move extended bed files into Bed_files.
      mkdir -p Bed_files
      cd Bed_files
      mv ../${prefix}/${prefix}.ext${extendlength}.bed ./
      cd ..

      # move mapping summary files into Mapping_summary.
      mkdir -p Summary
      cd Summary
      mv ../${prefix}/${prefix}_SE_KAS-seq_mapping_summary.txt ./
      cd ..
      echo "Move output files into folders. done."
      echo ""

   # map single-end KAS-seq data with bwa.
   elif [[ $aligner == "bwa" ]]; then

      # create $prefix folder and move read1 KAS-seq data into $prefix folder
      mkdir -p $prefix
      cd $prefix
      mv ../$read1 ./
	      
      # ouput the number of raw reads into mapping_summary txt file.
      # test if $read1 is .gz compressed file.
      if [[ ${read1##*.} == gz ]] ;then
         reads_num=$( zcat $read1 | wc -l | awk '{print $1/4}' )
      else
         reads_num=$( wc -l $read1 | awk '{print $1/4}' )
      fi

      echo "Number of KAS-seq reads: $reads_num"  >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_SE_KAS-seq_mapping_summary.txt

      echo "Map $read1 to $assemblyid reference genome with $aligner."
      echo ""
      # map $read1 to reference genome with bwa.
      bwa mem -t $threads $indexpath $read1 > ${prefix}.sam
      mapped_reads_num=$( grep chr HeLa-S3_spKAS-seq.rep1.sam | grep ^@ -v | awk '{print $1}' | sort -u | wc -l )
      echo "Number of mapped reads. ${aligner}: $mapped_reads_num" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      mapping_ratio=$( awk -v a=$reads_num -v b=$mapped_reads_num 'BEGIN {printf("%.2f\n",b*100/a)}' | awk '{print $1"%"}' )
      echo "Mapping ratios: $mapping_ratio" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      echo "$aligner mapping done."
      echo ""

      echo "Sort ${prefix}.sam file with 'samtools sort' and transfer into ${prefix}_sorted.bam file."
      echo ""
      samtools sort -@ 3 ${prefix}.sam -o ${prefix}_sorted.bam
      samtools index ${prefix}_sorted.bam 
      echo "'samtools sort' done."
      echo ""

      # remove duplicates with samtools rmdup.
      echo "Remove duplicates in ${prefix}_sorted.bam with 'samtools rmdup'"
      echo ""
      samtools rmdup -s ${prefix}_sorted.bam ${prefix}_rmdup.bam >> .${prefix}_SE_KAS-seq_deduplication_ratios.txt 2>&1
      samtools index ${prefix}_rmdup.bam
      unique_reads_num=$( samtools view -h ${prefix}_rmdup.bam | grep chr | grep ^@ -v | awk '{print $1}' | sort -u | wc -l )
      echo "Number of unique mapped reads: $unique_reads_num" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      duplication_ratios=$( awk '{print $6}' .${prefix}_SE_KAS-seq_deduplication_ratios.txt | awk '{printf("%s\n",$1*100"%")}' )
      echo "Duplication ratios. samtools: $duplication_ratios" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      rm -f .${prefix}_SE_KAS-seq_deduplication_ratios.txt
      echo "'samtools rmdup' done."
      echo ""
 
      # transfer bam into bed.
      echo "Transfer ${prefix}_rmdup.bam into ${prefix}.bed with bamToBed."
      echo ""
      bamToBed -i ${prefix}_rmdup.bam > ${prefix}.bed
      intersectBed -a ${prefix}.bed -b ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -v > ${prefix}.filter.bed
      mv ${prefix}.filter.bed ${prefix}.bed
      echo "'bamToBed' done."
      echo ""

      # extend the deduplicated mapped reads to extendlength.
      echo "Extend the deduplicated reads in ${prefix}.bed to ${extendlength}."
      echo ""
      awk '$3-'${extendlength}'>0 {if ($6~"+") printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$2+'${extendlength}',$4,$5,$6); else if ($6~"-") printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$3-'${extendlength}',$3,$4,$5,$6)}' ${prefix}.bed > ${prefix}.ext${extendlength}.bed
      echo "'extend' done."
      echo ""

      echo "Transfer ${prefix}.ext${extendlength}.bed into ${prefix}.ext${extendlength}.bg with genomeCoverageBed."
      echo ""
      genomeCoverageBed -bg -i ${prefix}.ext${extendlength}.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.ext${extendlength}.bg
      echo "'genomeCoverageBed' done."
      echo ""

      # clean up the .sam, sorted .bam and unextended .bed files.
      rm -f ${prefix}.sam
      rm -f ${prefix}_rmdup.bam.bai
      rm -f ${prefix}_rmdup.bam
      rm -f ${prefix}.bed
      echo "Clean up intermediate files. done."
      echo ""

      # move deduplicated bam files and index into Bam_files.
      cd ..
      mkdir -p Bam_files
      cd Bam_files
      mv ../${prefix}/${prefix}_sorted.bam ./
      mv ../${prefix}/${prefix}_sorted.bam.bai ./
      cd ..

      # move bedGraph files into BedGraph files.
      mkdir -p BedGraph_files
      cd BedGraph_files
      mv ../${prefix}/${prefix}.ext${extendlength}.bg ./
      cd ..

      # move extended bed files into Bed_files.
      mkdir -p Bed_files
      cd Bed_files
      mv ../${prefix}/${prefix}.ext${extendlength}.bed ./
      cd ..
      
      # move mapping summary files into Mapping_summary.
      mkdir -p Summary
      cd Summary
      mv ../${prefix}/${prefix}_SE_KAS-seq_mapping_summary.txt ./
      cd ..      
      echo "Move output files into folders. done."
      echo ""
   fi

# Map paired-end KAS-seq data.
elif [[ $paired_or_single_end == "paired" ]] ;then
   echo "Map paired-end KAS-seq data."
   echo ""
   # map paired-end KAS-seq data with bowtie2.
   if [[ $aligner == "bowtie2" ]]; then
      
      # create $prefix folder and move read1 and read2 KAS-seq data into $prefix folder
      mkdir -p $prefix
      cd $prefix
      mv ../${read1} ./
      mv ../${read2} ./

      # ouput the numebr of raw reads into mapping_summary txt file.
      # test if $read1 is .gz compressed file.
      if [[ ${read1##*.} == gz ]] ;then
         reads_num=$( zcat $read1 | wc -l | awk '{print $1/4}' )
      else
         reads_num=$( wc -l $read1 | awk '{print $1/4}' )
      fi

      echo "Number of KAS-seq reads pairs: $reads_num" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_PE_KAS-seq_mapping_summary.txt

      echo "Map ${read1} ${read2} to ${assemblyid} reference genome with ${aligner}."
      echo ""
      # map $read1 $read2 to reference genome with bowtie2.
      bowtie2 -X 1500 -p $threads -x $indexpath -1 $read1 -2 $read2 -S ${prefix}.sam >> .${prefix}_PE_KAS-seq_mapping_statistics.txt 2>&1 
      mapping_ratios=$( grep "overall alignment rate" .${prefix}_PE_KAS-seq_mapping_statistics.txt | awk '{print $1}' )
      mapped_reads_num=$( awk -v a=$reads_num -v b=$mapping_ratios 'BEGIN {printf("%d\n",a*b/100)}' )
      echo "Number of mapped reads. ${aligner}: $mapped_reads_num" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "Mapping ratios: $mapping_ratios" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      rm -f .${prefix}_PE_KAS-seq_mapping_statistics.txt
      echo "$aligner mapping done."
      echo ""

      # sort sam file with 'samtools sort' and transfer to bam file.
      echo "sort ${prefix}.sam file with 'samtools sort' and transfer into ${prefix}_sorted.bam file."
      echo ""
      sed '/^@PG/d' ${prefix}.sam | samtools sort -@ 3 - -o ${prefix}_sorted.bam
      samtools index ${prefix}_sorted.bam
      echo "'samtools sort' done."
      echo ""

      # remove duplicates with picard MarkDuplicates.
      echo "Remove duplicates in ${prefix}_sorted.bam with 'picard MarkDuplicates'"
      echo ""
      java -Xms512m -Xmx5g -jar ${SH_SCRIPT_DIR}/../src/picard.jar MarkDuplicates INPUT=${prefix}_sorted.bam OUTPUT=${prefix}_rmdup.bam METRICS_FILE=${prefix}.PCR_duplicates REMOVE_DUPLICATES=true
      samtools index ${prefix}_rmdup.bam
      echo "'picard MarkDuplicates' done."
      echo ""

      # estimate the exact size of KAS-seq fragments with mapped read1 and read2.
      echo "Combine 'properly paired' alignments into a single BED interval. ${prefix}_rmdup.bam into ${prefix}.bed."
      echo ""
      samtools view -h ${prefix}_rmdup.bam | ${SH_SCRIPT_DIR}/../src/SAMtoBED  -i - -o  ${prefix}.bed -x -v >> .${prefix}_PE_KAS-seq_fragment_length.txt 2>&1 
      paired_alignments_num=$( grep "Paired alignments (fragments):" .${prefix}_PE_KAS-seq_fragment_length.txt | awk '{printf("%d\n",$4/2)}' )
      unpaired_alignments_num=$( grep "Unpaired alignments:" .${prefix}_PE_KAS-seq_fragment_length.txt | awk '{printf("%d\n",$3/2)}' )
      alignments_num=$(( paired_alignments_num + unpaired_alignments_num ))
      echo "Number of unique mapped reads: $alignments_num" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
 
      duplication_ratios=$( grep "Unknown Library" ${prefix}.PCR_duplicates | awk '{printf("%.2f\n",$9*100)}' | awk '{print $1"%"}' )
      echo "Duplication ratios. picard: $duplication_ratios" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      rm -f ${prefix}.PCR_duplicates

      sed -i "/^Warning/d" .${prefix}_PE_KAS-seq_fragment_length.txt
      fragment_length=$( grep "Average fragment length:" .${prefix}_PE_KAS-seq_fragment_length.txt | awk '{printf("%s\n",$4)}' )
      echo "Length of DNA fragments: $fragment_length" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      rm -f .${prefix}_PE_KAS-seq_fragment_length.txt

      echo "'SAMtoBED' combines "properly paired" alignments done."
      echo "" 

      echo "Remove blacklist reads ..."
      echo ""
      bedSort ${prefix}.bed ${prefix}.sort.bed
      sed -i '/^chrM/d' ${prefix}.sort.bed
      intersectBed -a ${prefix}.sort.bed -b ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -v > ${prefix}.bed
      echo "done."
      echo ""

      echo "Attach strand information into combined alignments."
      echo ""
      samtools view -hbf 64 ${prefix}_rmdup.bam > ${prefix}_rmdup.R1.bam
      bamToBed -i ${prefix}_rmdup.R1.bam > ${prefix}_rmdup.R1.bed
      sed -i '/^chrM/d' ${prefix}_rmdup.R1.bed
      sed -i "s/\/1//g" ${prefix}_rmdup.R1.bed

      intersectBed -a ${prefix}.bed -b ${prefix}_rmdup.R1.bed -wa -wb -F 1 | awk '$4==$8 {printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$3,$8,$9,$10)}' > ${prefix}.6bed
      mv ${prefix}.6bed ${prefix}.bed

      echo "done."
      echo ""

      echo "Transfer ${prefix}.bed into ${prefix}.bg with genomeCoverageBed."
      echo ""
      genomeCoverageBed -bg -i ${prefix}.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.bg
      echo "'genomeCoverageBed' done."
      echo ""

      # clean up the .sam, sorted .bam and unextended .bed files.
      rm -f ${prefix}.sam
      rm -f ${prefix}_rmdup.bam.bai
      rm -f ${prefix}_rmdup.bam
      rm -f ${prefix}_rmdup.R1.bam
      rm -f ${prefix}_rmdup.R1.bed
      rm -f ${prefix}.sort.bed
      echo "Clean up intermediate files. done."
      echo ""

      # move deduplicated bam files and index into Bam_files.
      cd ..
      mkdir -p Bam_files
      cd Bam_files
      mv ../${prefix}/${prefix}_sorted.bam ./
      mv ../${prefix}/${prefix}_sorted.bam.bai ./
      cd ..

      # move bedGraph files into BedGraph files.
      mkdir -p BedGraph_files
      cd BedGraph_files
      mv ../${prefix}/${prefix}.bg ./
      cd ..

      # move extended bed files into Bed_files.
      mkdir -p Bed_files
      cd Bed_files
      mv ../${prefix}/${prefix}.bed ./
      cd ..

      # move mapping summary files into Mapping_summary.
      mkdir -p Summary
      cd Summary
      mv ../${prefix}/${prefix}_PE_KAS-seq_mapping_summary.txt ./
      cd ..
      echo "Move output files into folders. done."
      echo ""

   # map paired-end KAS-seq data with bwa.
   elif [[ $aligner == "bwa" ]] ;then
      
      # create $prefix folder and move read1 and read2 KAS-seq data into $prefix folder.
      mkdir -p $prefix
      cd $prefix
      mv ../$read1 ./
      mv ../$read2 ./

      # ouput the numebr of raw reads into mapping_summary txt file.
      # test if $read1 is .gz compressed file.
      if [[ ${read1##*.} == gz ]] ;then
         reads_num=$( zcat $read1 | wc -l | awk '{print $1/4}' )
      else
         reads_num=$( wc -l $read1 | awk '{print $1/4}' )
      fi
      
      echo "Number of KAS-seq reads pairs: $reads_num" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_PE_KAS-seq_mapping_summary.txt

      echo "Map ${read1} ${read2} to ${assemblyid} reference genome with ${aligner}."
      echo ""
      # map $read1 $read2 to reference genome with bwa.
      bwa mem -M -t $threads $indexpath $read1 $read2 > ${prefix}.sam
      mapped_reads_num=$( grep chr ${prefix}.sam | grep ^@ -v | awk '{print $1}' | sort -u | wc -l )
      echo "Number of mapped reads. ${aligner}: $mapped_reads_num" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      mapping_ratio=$( awk -v a=$reads_num -v b=$mapped_reads_num 'BEGIN {printf("%.2f\n",b*100/a)}' | awk '{print $1"%"}' )
      echo "Mapping ratios: $mapping_ratio" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "$aligner mapping done."
      echo ""

      # sort sam file with 'samtools sort' and transfer to bam file.
      # remove @PG lines from .sam file.
      echo "sort ${prefix}.sam file with 'samtools sort' and transfer into ${prefix}_sorted.bam file."
      echo ""
      sed '/^@PG/d' ${prefix}.sam | samtools sort -@ 3 - -o ${prefix}_sorted.bam
      samtools index ${prefix}_sorted.bam  
      echo "'samtools sort' done."
      echo ""

      # remove duplicates with picard MarkDuplicates.
      echo "Remove duplicates in ${prefix}_sorted.bam with 'picard MarkDuplicates'"
      echo ""
      java -Xms512m -Xmx5g -jar ${SH_SCRIPT_DIR}/../src/picard.jar MarkDuplicates INPUT=${prefix}_sorted.bam OUTPUT=${prefix}_rmdup.bam METRICS_FILE=${prefix}.PCR_duplicates REMOVE_DUPLICATES=true
      samtools index ${prefix}_rmdup.bam
      echo "'picard MarkDuplicates' done."
      echo ""

      # estimate the exact size of fragments from map read1 and read2.
      echo "Combine 'properly paired' alignments into a single BED interval. ${prefix}_rmdup.bam into ${prefix}.bed."
      echo ""
      samtools view -h  ${prefix}_rmdup.bam | ${SH_SCRIPT_DIR}/../src/SAMtoBED  -i - -o ${prefix}.bed -x -v >> .${prefix}_PE_KAS-seq_fragment_length.txt 2>&1
      paired_alignments_num=$( grep "Paired alignments (fragments):" .${prefix}_PE_KAS-seq_fragment_length.txt | awk '{printf("%d\n",$4/2)}' )
      unpaired_alignments_num=$( grep "Unpaired alignments:" .${prefix}_PE_KAS-seq_fragment_length.txt | awk '{printf("%d\n",$3/2)}' )
      alignments_num=$(( paired_alignments_num + unpaired_alignments_num ))
      echo "Number of unique mapped reads: $alignments_num" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      
      duplication_ratios=$( grep "Unknown Library" ${prefix}.PCR_duplicates | awk '{printf("%.2f\n",$9*100)}' | awk '{print $1"%"}' )
      echo "Duplication ratios. picard: $duplication_ratios" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      rm -f ${prefix}.PCR_duplicates

      sed -i "/^Warning/d" .${prefix}_PE_KAS-seq_fragment_length.txt
      fragment_length=$( grep "Average fragment length:" .${prefix}_PE_KAS-seq_fragment_length.txt | awk '{printf("%s\n",$4)}' )
      echo "Length of DNA fragments: $fragment_length" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      rm -f .${prefix}_PE_KAS-seq_fragment_length.txt
            
      echo "'SAMtoBED' combine 'properly paired' alignments done."
      echo "" 

      echo "Remove blacklist reads ..."
      echo ""
      bedSort ${prefix}.bed ${prefix}.sort.bed
      sed -i '/^chrM/d' ${prefix}.sort.bed
      intersectBed -a ${prefix}.sort.bed -b ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -v > ${prefix}.bed
      echo "done."
      echo ""

      echo "Attach strand information into combined alignments."
      echo ""
      samtools view -hbf 64 ${prefix}_rmdup.bam > ${prefix}_rmdup.R1.bam
      bamToBed -i ${prefix}_rmdup.R1.bam > ${prefix}_rmdup.R1.bed
      sed -i '/^chrM/d' ${prefix}_rmdup.R1.bed
      sed -i "s/\/1//g" ${prefix}_rmdup.R1.bed
      intersectBed -a ${prefix}.bed -b ${prefix}_rmdup.R1.bed -wa -wb -F 1 | awk '$4==$8 {printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$3,$8,$9,$10)}' > ${prefix}.6bed
      mv ${prefix}.6bed ${prefix}.bed
      echo "done."
      echo ""

      echo "Transfer ${prefix}.bed into ${prefix}.bg with genomeCoverageBed."
      echo ""
      genomeCoverageBed -bg -i ${prefix}.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.bg
      echo "'genomeCoverageBed' done."
      echo ""

      # clean up the .sam, sorted .bam and unextended .bed files.
      rm -f ${prefix}.sam
      rm -f ${prefix}_rmdup.bam.bai
      rm -f ${prefix}_rmdup.bam
      rm -f ${prefix}_rmdup.R1.bam
      rm -f ${prefix}_rmdup.R1.bed
      rm -f ${prefix}.sort.bed
      echo "Clean up intermediate files. done."
      echo "" 

      # move deduplicated bam files and index into Bam_files.
      cd ..
      mkdir -p Bam_files
      cd Bam_files
      mv ../${prefix}/${prefix}_sorted.bam ./
      mv ../${prefix}/${prefix}_sorted.bam.bai ./
      cd ..

      # move bedGraph files into BedGraph files.
      mkdir -p BedGraph_files
      cd BedGraph_files
      mv ../${prefix}/${prefix}.bg ./
      cd ..

      # move extended bed files into Bed_files.
      mkdir -p Bed_files
      cd Bed_files
      mv ../${prefix}/${prefix}.bed ./
      cd ..

      # move mapping summary files into Mapping_summary.
      mkdir -p Summary
      cd Summary
      mv ../${prefix}/${prefix}_PE_KAS-seq_mapping_summary.txt ./
      cd ..
      echo "Move output files into folders. done."
      echo ""
   fi
fi

echo "'KAS-pipe2 KAS-seq' run successfully!"
