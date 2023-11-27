#!/bin/bash
# 'KAS-Analyzer spKAS-seq' was developped by Ruitu Lyu on 12-08-2021.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-Analyzer spKAS-seq [ -h ] [ -t threads ] [ -i index path ] [ -u ] [ -r ] [ -f fold change ] [ -b bin size ] [ -e extend length ] [ -o prefix ] [ -s assembly id ] [ -1 read1 ] [ -2 read2 ]

Note: we strongly recommend paired-end sequencing for strand specific KAS-seq (spKAS-seq) data to accurately measure the fragments size."
exampleHelp="Example:
       Single-end:
       nohup KAS-Analyzer spKAS-seq -t 10 -i /absolute path/hg19_Bowtie2Index/hg19 -o spKAS-seq -r -s hg19 -1 spKAS-seq.trim.R1.fastq.gz &
       Paired-end:
       nohup KAS-Analyzer spKAS-seq -t 10 -i /absolute path/hg19_Bowtie2Index/hg19 -o spKAS-seq -r -s hg19 -1 spKAS-seq.trim.R1.fastq.gz -2 spKAS-seq.trim.R2.fastq.gz &
       Note: Bowtie2 Index example: /absolute path/Bowtie2Index/hg19."
threadsHelp="-t [threads]: please specify the number of threads used for spKAS-seq data mapping. DEFAULT: 1."
indexpathHelp="-i [index path]: please input the absolute path of reference genome index for aligner. Note: the path to a folder followed by a prefix of genome index. REQUIRED."
uniqueHelp="-u: please specify to filter the unique mapped reads. DEFAULT: off."
rloopsHelp="-r: please specify if identify R-loops regions with spKAS-seq data. DEFAULT: off."
foldchangeHelp="-f [fold change cutoff]: please specify the fold change cutoff of spKAS-seq reads difference between plus and minus strands used for R-loops identification. DEFAULT: 2."
binsizeHelp="-b [bin size]: please specify the size of bins used to identify R-loops. DEFAULT: 500."
extendlengthHelp="-e [extendlengthHelp]: please input the extend length for single-end spKAS-seq data. DEFAULT: 150."
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-Analyzer spKAS-seq' output files. REQUIRED."
assemblyidHelp="-s [assembly id]: please input the reference genome assembly id, e.g. Human: hg18, hg19, hg38, hs1; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. REQUIRED."
read1Help="-1 [read1]: please input trimmed single-end spKAS-seq fastq file or read1 of paired-end spKAS-seq fastq files; compressed .fastq.gz is accepted. REQUIRED."
read2Help="-2 [read2]: please input trimmed read2 of paired-end spKAS-seq raw fastq files. compressed fastq.gz is accepted."
helpHelp="-h: print this help and exit.
Note: The 'KAS-Analyzer spKAS-seq' shell script mainly invoke bowtie2 for spKAS-seq data mapping and R-loops identification, please refer to their official websites for more information."

# print help function.
printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$threadsHelp"
    echo -e ""
    echo -e "$indexpathHelp"
    echo -e ""
    echo -e "$uniqueHelp"
    echo -e ""
    echo -e "$rloopsHelp"
    echo -e ""
    echo -e "$foldchangeHelp"
    echo -e ""
    echo -e "$binsizeHelp"
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

# if no parameters was provided, 'KAS-Analyzer spKAS-seq' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'hrt:i:uf:b:e:o:s:1:2:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        r) rloops="true" ;;
        t) threads=$OPTARG ;;
        i) indexpath=$OPTARG ;;
	u) unique="on" ;;
	f) foldchange=$OPTARG ;;
	b) binsize=$OPTARG ;;
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
   echo "Please provide the path of reference genome index for bowtie2. -i [index path]"
   echo "Note: Bowtie2 Index example: /absolute path/Bowtie2Index/hg19; bwa Index example: /absolute path/BWAIndex/hg19.fa"
   echo ""
   exit -1
fi

if test -z $prefix ;then
   echo ""
   echo "Please provide the prefix (basename) of 'KAS-Analyzer spKAS-seq' output files. -o [prefix]"
   echo ""
   exit -1
fi

if test -z $assemblyid ;then
   echo ""	
   echo "Please input the reference genome assembly id, e.g. Human: hg18, hg19, hg38, hs1; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the bowtie2 reference genome index. -s [assembly id]."
   echo ""
   exit -1
fi

if test -z $read1 ;then
   echo ""
   echo "Please input the trimmed single-end spKAS-seq fastq file or read1 of paired-end spKAS-seq fastq file; compressed .fastq.gz is accepted. -1 [read1]"
   echo ""
   exit -1
fi


# setup the default options.
if test -z $threads ;then
   threads=1
fi

if test -z $unique ;then
   unique="off"
fi

if test -z $rloops ;then
   rloops="false"
fi

if test -z $foldchange ;then
   foldchange=2
fi

if test -z $binsize ;then
   binsize=500
fi

if test -z $extendlength ;then
   extendlength=150
fi

# the path of KAS-Analyzer spKAS-seq script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# supported assembly id.
if [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "hs1" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]] ;then
   echo ""
   echo "Error: unsupported assembly id: $assemblyid. Supported assembly id: Human: hg18, hg19, hg38, hs1; Mouse: mm9,mm10,mm39; Fruitfly: dm3, dm6; Rat: rn6, rn7; C.elegans: ce10, ce11; Zebra fish: danRer10, danRer11."
   echo ""
   exit 1
fi

# test if aligners was installed.
if ! type bowtie2 >/dev/null 2>&1; then
   echo "bowtie2 was not installed or not export to the \$PATH'"
   echo ""
   echo "Install bowtie2 with 'conda install -c bioconda bowtie2' or refer the official website of 'bowtie2'."
   echo ""
   exit 1
fi

if ! type samtools >/dev/null 2>&1; then
   echo "samtools was not installed or not export to the \$PATH'"
   echo ""
   echo "Install samtools with 'conda install -c bioconda samtools' or refer the official website of 'samtools'."
   echo ""
   exit 1
fi

if ! type bedtools >/dev/null 2>&1; then
   echo "bedtools was not installed or not export to the \$PATH'"
   echo ""
   echo "Install bedtools with 'conda install -c bioconda bedtools' or refer the official website of 'bedtools'."
   echo ""
   exit 1
fi

if ! type macs2 >/dev/null 2>&1; then
   echo "macs2 was not installed or not export to the \$PATH'"
   echo ""
   echo "Install macs2 with 'conda install -c bioconda macs2' or refer the official website of 'macs2'."
   echo ""
   exit 1
fi

# Determine single or paired KAS-seq data.
if test -z $read2 ;then
   paired_or_single_end="single"
else
   paired_or_single_end="paired"
fi


echo ""
echo "Welcome to analyze spKAS-seq data with 'KAS-Analyzer spKAS-seq'... "
echo ""

# Map single-end spKAS-seq data.
if [[ $paired_or_single_end == "single" ]] ;then
   echo "Map single-end $read1 spKAS-seq data to $assemblyid genome with bowtie2 ..."
   echo ""

   # create $prefix folder and move read1 spKAS-seq data into $prefix folder
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

   echo "Number of spKAS-seq reads: $reads_num" >> ${prefix}_SE_spKAS-seq_mapping_summary.txt
   echo "" >> ${prefix}_SE_spKAS-seq_mapping_summary.txt
   bowtie2 -p $threads -x $indexpath $read1 -S ${prefix}.sam >> .${prefix}_SE_spKAS-seq_mapping_statistics.txt 2>&1
   mapping_ratios=$( grep "overall alignment rate" .${prefix}_SE_spKAS-seq_mapping_statistics.txt | awk '{print $1}' )
   mono_mapped_reads_num=$( grep "aligned exactly 1 time" .${prefix}_SE_spKAS-seq_mapping_statistics.txt | awk '{print $1}' )
   multi_mapped_reads_num=$( grep "aligned >1 times" .${prefix}_SE_spKAS-seq_mapping_statistics.txt | awk '{print $1}' )
   mapped_reads_num=$(( mono_mapped_reads_num + multi_mapped_reads_num ))
   echo "Number of mapped reads. bowtie2: $mapped_reads_num" >> ${prefix}_SE_spKAS-seq_mapping_summary.txt
   echo "" >> ${prefix}_SE_spKAS-seq_mapping_summary.txt
   echo "Mapping ratios: $mapping_ratios" >> ${prefix}_SE_spKAS-seq_mapping_summary.txt
   echo "" >> ${prefix}_SE_spKAS-seq_mapping_summary.txt
   rm -f .${prefix}_SE_spKAS-seq_mapping_statistics.txt
   echo "bowtie2 mapping done."
   echo ""

   # sort sam file with 'samtools sort' and transfer to bam file.
   echo "Sort ${prefix}.sam file with 'samtools sort' and transfer into ${prefix}_sorted.bam file."
   echo ""
   samtools sort -@ 3 ${prefix}.sam -o ${prefix}_sorted.bam
   samtools index ${prefix}_sorted.bam
   echo "'samtools sort' done."
   echo ""

   # remove duplicates with samtools rmdup.
   echo "Remove duplicates in ${prefix}_sorted.bam with 'samtools rmdup'"
   echo ""
   samtools rmdup -s ${prefix}_sorted.bam ${prefix}_rmdup.bam >> .${prefix}_SE_spKAS-seq_deduplication_ratios.txt 2>&1
   samtools index ${prefix}_rmdup.bam
   deduplicated_reads_num=$( samtools view -h ${prefix}_rmdup.bam | grep chr | grep ^@ -v | awk '{print $1}' | sort -u | wc -l )
   echo "Number of deduplicated mapped reads: $deduplicated_reads_num" >> ${prefix}_SE_spKAS-seq_mapping_summary.txt
   echo "" >> ${prefix}_SE_spKAS-seq_mapping_summary.txt
   duplication_ratios=$( awk '{print $6}' .${prefix}_SE_spKAS-seq_deduplication_ratios.txt | awk '{printf("%s\n",$1*100"%")}' )
   echo "Duplication ratios. samtools: $duplication_ratios" >> ${prefix}_SE_spKAS-seq_mapping_summary.txt
   rm -f .${prefix}_SE_spKAS-seq_deduplication_ratios.txt
   echo "'samtools rmdup' done."
   echo ""
   
   if [[ $unique == "on" ]] ;then
      echo "Filter the unique mapped reads ..."
      echo ""
      samtools view -q 10 -b ${prefix}_rmdup.bam | bamToBed -i - | awk '$3-'${extendlength}'>0 {if ($6~"+") printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$2+'${extendlength}',$4,$5,$6); else if ($6~"-") printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$3-'${extendlength}',$3,$4,$5,$6)}' | intersectBed -a - -b ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -v | intersectBed -a - -b ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes.bed -wa -f 1 > ${prefix}.ext${extendlength}.unique.bed
      
      # calculate the number of deduplicated and unique mapped reads.
      unique_mapped_reads_num=$( wc -l ${prefix}.ext${extendlength}.unique.bed | awk '{print $1}' )
      echo "" >> ${prefix}_SE_KAS-seq_mapping_summary.txt
      echo "Number of unique mapped reads: $unique_mapped_reads_num" >> ${prefix}_SE_KAS-seq_mapping_summary.txt     

      genomeCoverageBed -bg -i ${prefix}.ext${extendlength}.unique.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.ext${extendlength}.unique.bg
      echo "done."
      echo ""
   fi 

   # transfer bam into bed.
   echo "Transfer ${prefix}_rmdup.bam into ${prefix}.bed with bamToBed."
   echo ""
   bamToBed -i ${prefix}_rmdup.bam > ${prefix}.bed
   intersectBed -a ${prefix}.bed -b ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -v > ${prefix}.filter.bed
   mv ${prefix}.filter.bed ${prefix}.bed
   echo "'bamToBed' done."
   echo ""

   # extend the single-end spKAS-seq deduplicated mapped reads to extendlength.
   echo "Extend the deduplicated reads in ${prefix}.bed to ${extendlength}."
   echo ""
   awk '$3-'${extendlength}'>0 {if ($6~"+") printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$2+'${extendlength}',$4,$5,$6); else if ($6~"-") printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$3-'${extendlength}',$3,$4,$5,$6)}' ${prefix}.bed | sortBed -i > ${prefix}.ext${extendlength}.bed
   echo "'extend' done."
   echo ""
   
   echo "Transfer ${prefix}.ext${extendlength}.bed into ${prefix}.ext${extendlength}.bg with genomeCoverageBed."
   echo ""
   genomeCoverageBed -bg -i ${prefix}.ext${extendlength}.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.ext${extendlength}.bg
   echo "'genomeCoverageBed' done."
   echo ""
     
   # Identify R-loops. 
   if [[ $rloops == "true" ]]; then
      echo "Identify R-loops ..."
      echo ""
        
      # filter out extended mapped reads that exceed the start or end of chromosome.
      intersectBed -a ${prefix}.ext${extendlength}.bed -b ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes.bed -wa -f 1 | sortBed -i > ${prefix}.ext${extendlength}.sort.bed
     
      # filter out chrM and unfinished random chromosome reads, separate spKAS-seq reads into minus and plus mapped reads.
      echo "Separate extended unique $read1 mapped reads into plus and minus and remove chrM."
      echo ""
      grep ^chrM -v ${prefix}.ext${extendlength}.sort.bed | grep + > ${prefix}.ext${extendlength}.plus.bed
      grep ^chrM -v ${prefix}.ext${extendlength}.sort.bed | grep + -v > ${prefix}.ext${extendlength}.minus.bed
      echo "done."
      echo ""
     
      # transfer minus or plus bed filed into bedGraph files.
      echo "Transfer ${prefix}.ext${extendlength}.plus.bed or .minus.bed to ${prefix}.ext${extendlength}.plus.bed or .minus.bg with genomeCoverageBed."
      echo ""
      genomeCoverageBed -bg -i ${prefix}.ext${extendlength}.plus.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.ext${extendlength}.plus.bg 
      genomeCoverageBed -bg -i ${prefix}.ext${extendlength}.minus.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.ext${extendlength}.minus.bg
      echo "done."
      echo ""
     
      # transfer minus or plus bedGraph files into bigWig files.
      echo "Transfer ${prefix}.ext${extendlength}.plus.bg and .minus.bg to .bigWig with bedGraphToBigWig."
      echo ""
      bedGraphToBigWig ${prefix}.ext${extendlength}.plus.bg ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes ${prefix}.ext${extendlength}.plus.bigWig
      bedGraphToBigWig ${prefix}.ext${extendlength}.minus.bg ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes ${prefix}.ext${extendlength}.minus.bigWig 
      echo "done."
      echo ""

      # generate the $binsize bp bins with $binsize/2 bp overlap and calculate the plus or minus mapped averaged spKAS-seq reads density.
      echo "Calculate plus and minus spKAS-seq density on bins of $binsize."
      echo ""
      distanceBetweenBins=$((-1*binsize/2))
      multiBigwigSummary bins -b ${prefix}.ext${extendlength}.minus.bigWig ${prefix}.ext${extendlength}.plus.bigWig --labels minus plus -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -bs $binsize -n $distanceBetweenBins -p $threads -out ${prefix}_plus_vs_minus.bins.rmbl.npz --outRawCounts ${prefix}_plus_vs_minus.bins.rmbl.tab 
     
      sed -i "s/nan/0/g" ${prefix}_plus_vs_minus.bins.rmbl.tab
      sed "1d" ${prefix}_plus_vs_minus.bins.rmbl.tab | awk '{printf("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n",$1,$2,$3,$4,$5,log(($5+0.1)/($4+0.1))/log(2))}' > ${prefix}_plus_vs_minus.bins.rmbl.bed
      echo "done."
      echo ""

      # peak calling for spKAS-seq data.
      echo "Call spKAS-seq peaks."
      echo ""
      sed -i '/^chrM/d' ${prefix}.ext${extendlength}.bed
 
      if [[ $assemblyid == "hg18" ]] || [[ $assemblyid == "hg19" ]] || [[ $assemblyid == "hg38" ]] || [[ $assemblyid == "hs1" ]] ;then
             assemblysize="2.7e9"
      elif [[ $assemblyid == "mm9" ]] || [[ $assemblyid == "mm10" ]] || [[ $assemblyid == "mm39" ]] ;then
             assemblysize="1.87e9"
      elif [[ $assemblyid == "ce10" ]] || [[ $assemblyid == "ce11" ]] ;then
             assemblysize="9e7"
      elif [[ $assemblyid == "dm3" ]] || [[ $assemblyid == "dm6" ]] ;then
             assemblysize="1.2e8"
      elif [[ $assemblyid == "rn6" ]] || [[ $assemblyid == "rn7" ]] ;then
             assemblysize="2.1e9"	     
      elif [[ $assemblyid == "danRer10" ]] || [[ $assemblyid == "danRer11" ]] ;then
	     assemblysize="9.5e8"
      fi

      macs2 callpeak -t ${prefix}.ext${extendlength}.bed -n ${prefix} --broad -g ${assemblysize} --broad-cutoff 0.01 -q 0.01 
      echo "spKAS-seq peak calling done."
      echo ""

      # filter bins overlap with spKAS-seq peaks.
      echo "Filter bins on spKAS-seq peaks."
      echo ""
      intersectBed -a ${prefix}_plus_vs_minus.bins.rmbl.bed -b ${prefix}_peaks.broadPeak -wa -f 0.5 | sort -u | sortBed -i > ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed
      echo "done."
      echo ""

      # identify bins with $foldchange different spKAS-seq reads density between plus and minus strands.
      echo "Identify R-loops bins based on the difference of spKAS-seq density on plus and minus strands."
      echo ""
      awk '$6>= '${foldchange}' {printf("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n",$1,$2,$3,$4,$5,$6)}' ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed > R-loop_plus_${prefix}.bins.rmbl.bed
      awk '$6<= -'${foldchange}' {printf("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n",$1,$2,$3,$4,$5,$6)}' ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed > R-loop_minus_${prefix}.bins.rmbl.bed
      echo "done."
      echo ""

      # merge R-loops bins into R-loops peaks.
      echo "Merge R-loops bins."
      echo ""
      mergeBed -i R-loop_plus_${prefix}.bins.rmbl.bed | awk '{printf("%s\t%d\t%d\t%s\n",$1,$2,$3,"+")}' > ${prefix}_R-loops.plus.bed
      mergeBed -i R-loop_minus_${prefix}.bins.rmbl.bed | awk '{printf("%s\t%d\t%d\t%s\n",$1,$2,$3,"-")}' > ${prefix}_R-loops.minus.bed
      echo "done."
      echo ""
     
      echo "Generate R-loops density at 50bp resolution."
      echo ""
      multiBigwigSummary bins -b ${prefix}.ext${extendlength}.minus.bigWig ${prefix}.ext${extendlength}.plus.bigWig --labels minus plus -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -bs 50 -p $threads -out ${prefix}_plus_vs_minus.50bp.bins.rmbl.npz --outRawCounts ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab 

      sed -i "s/nan/0/g" ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab
      sed "1d" ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab | awk '$5-$4>0 {printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$5-$4)}' > ${prefix}_R-loop.density.plus.bg
      sed "1d" ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab | awk '$4-$5>0 {printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$5-$4)}' > ${prefix}_R-loop.density.minus.bg
     
      cat ${prefix}_R-loop.density.plus.bg ${prefix}_R-loop.density.minus.bg | sortBed -i > ${prefix}_R-loop.density.bg
      echo "done."
      echo ""
#     bedGraphToBigWig ${prefix}_R-loop.density.bg ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes ${prefix}_R-loop.density.bigWig

      # clean up intermediate files generated during R-loops identification.
      rm -f ${prefix}.ext${extendlength}.sort.bed
      rm -f ${prefix}.ext${extendlength}.plus.bed
      rm -f ${prefix}.ext${extendlength}.minus.bed
      rm -f ${prefix}.ext${extendlength}.plus.bg
      rm -f ${prefix}.ext${extendlength}.minus.bg
      rm -f ${prefix}.ext${extendlength}.plus.bigWig
      rm -f ${prefix}.ext${extendlength}.minus.bigWig
      rm -f ${prefix}_plus_vs_minus.bins.rmbl.npz
      rm -f ${prefix}_plus_vs_minus.bins.rmbl.tab
      rm -f ${prefix}_plus_vs_minus.bins.rmbl.bed
      rm -f ${prefix}_model.r
      rm -f ${prefix}_peaks.broadPeak
      rm -f ${prefix}_peaks.gappedPeak
      rm -f ${prefix}_peaks.xls
      rm -f ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed
      rm -f R-loop_plus_${prefix}.bins.rmbl.bed
      rm -f R-loop_minus_${prefix}.bins.rmbl.bed
      # clean up intermediate files generated during R-loops density.
      rm -f ${prefix}_plus_vs_minus.50bp.bins.rmbl.npz
      rm -f ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab
      echo "Clean up R-loops related intermediate files. done."
      echo ""

      # move R-loops and R-loop density (bigWig files) into R-loops folder.
      cd ..
      mkdir -p R-loops
      cd R-loops
      mv ../${prefix}/${prefix}_R-loops.plus.bed ./
      mv ../${prefix}/${prefix}_R-loops.minus.bed ./
      mv ../${prefix}/${prefix}_R-loop.density.plus.bg ./
      mv ../${prefix}/${prefix}_R-loop.density.minus.bg ./
      mv ../${prefix}/${prefix}_R-loop.density.bg ./
#     mv ../${prefix}/${prefix}_R-loop.density.bigWig ./
      cd ..
      cd ${prefix}
      echo "Move R-loops identification output files into 'R-loops' folder. done."
      echo ""
      fi

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
   mv ../${prefix}/${prefix}.ext${extendlength}.*bg ./
   cd ..

   # move extended bed files into Bed_files.
   mkdir -p Bed_files
   cd Bed_files
   mv ../${prefix}/${prefix}.ext${extendlength}.*bed ./
   cd ..

   # move mapping summary files into Mapping_summary.
   mkdir -p Summary
   cd Summary
   mv ../${prefix}/${prefix}_SE_spKAS-seq_mapping_summary.txt ./
   cd ..
   echo "Move output files into folders. done."
   echo ""

elif [[ $paired_or_single_end == "paired" ]]; then

   # create $prefix folder and move read1 and read2 KAS-seq data into $prefix folder
   mkdir -p $prefix
   cd $prefix
   mv ../$read1 ./
   mv ../$read2 ./

   # output the numebr of spKAS-seq raw reads into mapping summary txt file.
   # test if $read1 is .gz compressed file.
   if [[ ${read1##*.} == gz ]] ;then
      reads_num=$( zcat $read1 | wc -l | awk '{print $1/4}' )
   else
      reads_num=$( wc -l $read1 | awk '{print $1/4}' )
   fi

   echo "Number of spKAS-seq reads pairs: $reads_num" >> ${prefix}_PE_spKAS-seq_mapping_summary.txt
   echo "" >> ${prefix}_PE_spKAS-seq_mapping_summary.txt

   # map paired-end spKAS-seq data to reference genome with bowtie2.
   echo "Map paired-end $read1 $read2 spKAS-seq data to $assemblyid genome with bowtie2..."
   echo ""
   bowtie2 -X 1000 -p $threads -x $indexpath -1 $read1 -2 $read2 -S ${prefix}.sam >> .${prefix}_PE_spKAS-seq_mapping_statistics.txt 2>&1
   mapping_ratios=$( grep "overall alignment rate" .${prefix}_PE_spKAS-seq_mapping_statistics.txt | awk '{print $1}' )
   mapped_reads_num=$( awk -v a=$reads_num -v b=$mapping_ratios 'BEGIN {printf("%d\n",a*b/100)}' )
   echo "Number of mapped reads. bowtie2: $mapped_reads_num" >> ${prefix}_PE_spKAS-seq_mapping_summary.txt
   echo "" >> ${prefix}_PE_spKAS-seq_mapping_summary.txt
   echo "Mapping ratios: $mapping_ratios" >> ${prefix}_PE_spKAS-seq_mapping_summary.txt
   echo "" >> ${prefix}_PE_spKAS-seq_mapping_summary.txt
   rm -f .${prefix}_PE_spKAS-seq_mapping_statistics.txt
   echo "bowtie2 mapping done."
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
   samtools view -h ${prefix}_rmdup.bam | ${SH_SCRIPT_DIR}/../src/SAMtoBED  -i - -o ${prefix}.bed -x -v >> .${prefix}_PE_spKAS-seq_fragment_length.txt 2>&1
   paired_alignments_num=$( grep "Paired alignments (fragments):" .${prefix}_PE_spKAS-seq_fragment_length.txt | awk '{printf("%d\n",$4/2)}' )
   unpaired_alignments_num=$( grep "Unpaired alignments:" .${prefix}_PE_spKAS-seq_fragment_length.txt | awk '{printf("%d\n",$3/2)}' )
   alignments_num=$(( paired_alignments_num + unpaired_alignments_num ))
   echo "Number of deduplicated mapped reads: $alignments_num" >> ${prefix}_PE_spKAS-seq_mapping_summary.txt
   echo "" >> ${prefix}_PE_spKAS-seq_mapping_summary.txt

   duplication_ratios=$( grep "Unknown Library" ${prefix}.PCR_duplicates | awk '{printf("%.2f\n",$9*100)}' | awk '{print $1"%"}' )
   echo "Duplication ratios. picard: $duplication_ratios" >> ${prefix}_PE_spKAS-seq_mapping_summary.txt
   echo "" >> ${prefix}_PE_spKAS-seq_mapping_summary.txt
   rm -f ${prefix}.PCR_duplicates

   if [[ $unique == "on" ]] ;then
      echo "Filter the unique mapped reads ..."
      echo ""
      samtools view -q 10 ${prefix}_rmdup.bam | ${SH_SCRIPT_DIR}/../src/SAMtoBED -i - -o ${prefix}.unique.bed -x -v >> /dev/null 2>&1
      intersectBed -a ${prefix}.unique.bed -b ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -v | sortBed -i > ${prefix}.unique.rmbl.bed
      mv ${prefix}.unique.rmbl.bed ${prefix}.unique.bed
      
      # calculate the number of deduplicated and unique mapped reads.
      unique_mapped_reads_num=$( wc -l ${prefix}.unique.bed | awk '{print $1}' )
      echo "Number of unique mapped reads: $unique_mapped_reads_num" >> ${prefix}_PE_KAS-seq_mapping_summary.txt
      echo "" ${prefix}_PE_KAS-seq_mapping_summary.txt

      genomeCoverageBed -bg -i ${prefix}.unique.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.unique.bg
      echo "done."
      echo ""
   fi

   sed -i "/^Warning/d" .${prefix}_PE_spKAS-seq_fragment_length.txt
   fragment_length=$( grep "Average fragment length:" .${prefix}_PE_spKAS-seq_fragment_length.txt | awk '{printf("%s\n",$4)}' )
   echo "Length of DNA fragments: $fragment_length" >> ${prefix}_PE_spKAS-seq_mapping_summary.txt
   echo "" >> ${prefix}_PE_spKAS-seq_mapping_summary.txt
   rm -f .${prefix}_PE_spKAS-seq_fragment_length.txt

   echo "'SAMtoBED' combines "properly paired" alignments done."
   echo "" 

   echo "Transfer ${prefix}.bed into ${prefix}.bg with genomeCoverageBed."
   echo ""
   grep "^chrM" -v ${prefix}.bed | intersectBed -a - -b ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -v > ${prefix}.rmbl.bed
   echo "done."
   echo ""

   echo "Attach strand information into combined alignments."
   echo ""
   samtools view -hbf 64 ${prefix}_rmdup.bam > ${prefix}_rmdup.R1.bam
   bamToBed -i ${prefix}_rmdup.R1.bam | grep "^chrM" -v | sed "s/\/1//g" > ${prefix}_rmdup.R1.bed
   intersectBed -a ${prefix}.rmbl.bed -b ${prefix}_rmdup.R1.bed -wa -wb -F 1 | awk '$4==$8 {printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$3,$8,$9,$10)}' | sortBed -i > ${prefix}.bed
   echo "done."
   echo ""

   if [[ $unique == "on" ]] ;then
      echo "Attach strand information into combined alignments for uniquely mapped reads ..."
      echo "done."
      intersectBed -a ${prefix}.unique.bed -b ${prefix}_rmdup.R1.bed -wa -wb -F 1 | awk '$4==$8 {printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$3,$8,$9,$10)}' | sortBed -i > ${prefix}.unique.6bed
      mv ${prefix}.unique.6bed ${prefix}.unique.bed
      echo "done."
      echo ""
   fi  

   genomeCoverageBed -bg -i ${prefix}.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.bg
   echo "'genomeCoverageBed' done."
   echo ""
     
   # Identify R-loops. 
   if [[ $rloops == "true" ]]; then
      echo "Identify R-loops ..."
      echo ""	     
     
      echo "Separate combined alignments into minus and plus."
      echo ""
      grep '+' ${prefix}.bed | sortBed -i > ${prefix}.plus.bed
      grep '+' -v ${prefix}.bed | sortBed -i > ${prefix}.minus.bed
      echo "done."
      echo ""

      # transfer minus or plus bed filed into bedGraph files.
      echo "Transfer ${prefix}.plus.bed or .minus.bed into .bg files."
      echo ""
      genomeCoverageBed -bg -i ${prefix}.plus.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.plus.bg
      genomeCoverageBed -bg -i ${prefix}.minus.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.minus.bg
      echo "done."
      echo ""

      # transfer minus or plus bedGraph files into bigWig files.
      echo "Transfer ${prefix}.minus.bg or .plus.bg files into .bigWig files."
      echo ""
      bedGraphToBigWig ${prefix}.plus.bg ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes ${prefix}.plus.bigWig
      bedGraphToBigWig ${prefix}.minus.bg ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes ${prefix}.minus.bigWig
      echo "done."
      echo ""

      # generate the $binsize bp bins with $binsize/2 bp overlap and calculate the plus or minus mapped averaged spKAS-seq reads density.
      echo "Calculate the plus or minus mapped spKAS-seq density on ${binsize}bp bins."
      echo ""
      distanceBetweenBins=$((-1*binsize/2))
      multiBigwigSummary bins -b ${prefix}.minus.bigWig ${prefix}.plus.bigWig --labels minus plus -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -bs $binsize -n $distanceBetweenBins -p $threads -out ${prefix}_plus_vs_minus.bins.rmbl.npz --outRawCounts ${prefix}_plus_vs_minus.bins.rmbl.tab
      sed "s/nan/0/g" ${prefix}_plus_vs_minus.bins.rmbl.tab | sed "1d" | awk '{printf("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n",$1,$2,$3,$4,$5,log(($5+0.1)/($4+0.1))/log(2))}' > ${prefix}_plus_vs_minus.bins.rmbl.bed
      echo "done."
      echo ""

      # peak calling for spKAS-seq data.
      echo "spKAS-seq peaking with macs2 to filter ${binsize}bp bins."
      echo ""
      if [[ $assemblyid == "hg18" ]] || [[ $assemblyid == "hg19" ]] || [[ $assemblyid == "hg38" ]] || [[ $assemblyid == "hs1" ]]  ;then
	     assemblysize="2.7e9"
      elif [[ $assemblyid == "mm9" ]] || [[ $assemblyid == "mm10" ]] || [[ $assemblyid == "mm39" ]] ;then
	     assemblysize="1.87e9"
      elif [[ $assemblyid == "ce10" ]] || [[ $assemblyid == "ce11" ]] ;then
	     assemblysize="9e7"
      elif [[ $assemblyid == "dm3" ]] || [[ $assemblyid == "dm6" ]] ;then     
	     assemblysize="1.2e8"
      elif [[ $assemblyid == "rn6" ]] || [[ $assemblyid == "rn7" ]] ;then
             assemblysize="2.1e9"
      elif [[ $assemblyid == "danRer10" ]] || [[ $assemblyid == "danRer11" ]] ;then
             assemblysize="9.5e8"
      fi

      macs2 callpeak -t ${prefix}.bed -n ${prefix} --broad -g $assemblysize --broad-cutoff 0.01 -q 0.01
      echo "done."
      echo ""

      # filter bins overlap with spKAS-seq peaks.
      echo "Filter bins on spKAS-seq peaks."
      echo ""
      intersectBed -a ${prefix}_plus_vs_minus.bins.rmbl.bed -b ${prefix}_peaks.broadPeak -wa -f 0.5 | sort -u | sortBed -i > ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed
      echo "done."
      echo ""

      # identify bins with $foldchange different spKAS-seq reads density between plus and minus strands.
      echo "Identify R-loop bins on plus or minus strands."
      echo ""
      awk '$6>='$foldchange' {printf("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n",$1,$2,$3,$4,$5,$6)}' ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed > R-loop_plus_${prefix}.bins.rmbl.bed
      awk '$6<=-'$foldchange' {printf("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n",$1,$2,$3,$4,$5,$6)}' ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed > R-loop_minus_${prefix}.bins.rmbl.bed
      echo "done."
      echo ""

      # merge R-loops bins into R-loops peaks.
      echo "Merge R-loop bins into R-loop enriched regions."
      echo ""
      sortBed -i R-loop_plus_${prefix}.bins.rmbl.bed | mergeBed -i - | awk '{printf("%s\t%d\t%d\t%s\n",$1,$2,$3,"+")}' > ${prefix}_R-loops.plus.bed
      sortBed -i R-loop_minus_${prefix}.bins.rmbl.bed | mergeBed -i - | awk '{printf("%s\t%d\t%d\t%s\n",$1,$2,$3,"-")}' > ${prefix}_R-loops.minus.bed
      echo "done."
      echo ""

      echo "Generate the R-loop density at 50bp resolution ..."
      echo ""
      multiBigwigSummary bins -b ${prefix}.minus.bigWig ${prefix}.plus.bigWig --labels minus plus -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -bs 50 -p $threads -out ${prefix}_plus_vs_minus.50bp.bins.rmbl.npz --outRawCounts ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab

      sed -i "s/nan/0/g" ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab
      sed "1d" ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab | awk '$5-$4>0 {printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$5-$4)}' > ${prefix}_R-loop.density.plus.bg
      sed "1d" ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab | awk '$4-$5>0 {printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$5-$4)}' > ${prefix}_R-loop.density.minus.bg

      cat ${prefix}_R-loop.density.plus.bg ${prefix}_R-loop.density.minus.bg | sortBed -i > ${prefix}_R-loop.density.bg
      echo "done."
      echo ""
#     bedGraphToBigWig ${prefix}_R-loop.density.bg ${SH_SCRIPT_DIR}/../chrom_size/${assembly}.chrom.sizes ${prefix}_R-loop.density.bigWig

      # clean up intermediate files generated during R-loops identification.
      rm -f ${prefix}.plus.bed
      rm -f ${prefix}.minus.bed
      rm -f ${prefix}.plus.bg
      rm -f ${prefix}.minus.bg
      rm -f ${prefix}.plus.bigWig
      rm -f ${prefix}.minus.bigWig
      rm -f ${prefix}_plus_vs_minus.bins.rmbl.npz
      rm -f ${prefix}_plus_vs_minus.bins.rmbl.tab
      rm -f ${prefix}_plus_vs_minus.bins.rmbl.bed
      rm -f ${prefix}_model.r
      rm -f ${prefix}_peaks.broadPeak
      rm -f ${prefix}_peaks.gappedPeak
      rm -f ${prefix}_peaks.xls
      rm -f ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed
      rm -f R-loop_plus_${prefix}.bins.rmbl.bed
      rm -f R-loop_minus_${prefix}.bins.rmbl.bed
      # clean up intermediate files generated during R-loops density.
      rm -rf ${prefix}_plus_vs_minus.50bp.bins.rmbl.npz
      rm -rf ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab
      echo "Clean up R-loops related intermediate files. done."
      echo ""

      # move R-loops and R-loops signal density (bigWig file) into R-loops directory.
      cd ..
      mkdir -p R-loops
      cd R-loops
      mv ../${prefix}/${prefix}_R-loops.plus.bed ./
      mv ../${prefix}/${prefix}_R-loops.minus.bed ./
      mv ../${prefix}/${prefix}_R-loop.density.plus.bg ./
      mv ../${prefix}/${prefix}_R-loop.density.minus.bg ./
      mv ../${prefix}/${prefix}_R-loop.density.bg ./
#     mv ../${prefix}/${prefix}_R-loop.density.bigWig ./
      cd ..
      cd ${prefix}
      echo "Move R-loops identification output files into 'R-loops' folder. done."
      echo ""
   fi

   # clean up the .sam, sorted .bam and unextended .bed files.i
   rm -f ${prefix}.sam
   rm -f ${prefix}_rmdup.bam.bai
   rm -f ${prefix}_rmdup.bam
   rm -f ${prefix}.rmbl.bed
   rm -f ${prefix}_rmdup.R1.bam
   rm -f ${prefix}_rmdup.R1.bed
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
   mv ../${prefix}/${prefix}.*bg ./
   cd ..

   # move extended bed files into Bed_files.
   mkdir -p Bed_files
   cd Bed_files
   mv ../${prefix}/${prefix}.*bed ./
   cd ..

   # move mapping summary files into Mapping_summary.
   mkdir -p Summary
   cd Summary
   mv ../${prefix}/${prefix}_PE_spKAS-seq_mapping_summary.txt ./
   cd ..
   echo "Move output files into folders. done."
   echo ""

fi 

echo "'KAS-Analyzer spKAS-seq' run successfully."
