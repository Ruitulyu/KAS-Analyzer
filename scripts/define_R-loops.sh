#!/bin/bash
# 'KAS-pipe2 R-loop' was developed by Ruitu Lyu on 12-22-2021.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 R-loop [ -h/--help ] [ -t threads ] [ -o prefix ] [ -s assembly id ] [ -p peaks ] [ -b bin size ] [ -f fold change ] [ -l labels ] [ -n Input ] [ -k spKAS-seq ]"
exampleHelp="Example: nohup KAS-pipe2 R-loop -o KAS-seq_R-loops -t 10 -s mm10 -l labels.txt -k KAS-seq.txt &"
threadsHelp="-t [threads]: please specify the number of threads used for R-loops identification. DEFAULT: 1."
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 R-loop' output files. Default: basename of txt files containing spKAS-seq data."
assemblyidHelp="-s [assembly id]: please specify the genome assembly id of spKAS-seq data. -s [assembly id]. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11."
peaksHelp="-p [peaks]: please specify the spKAS-seq peaks file. if not specified, spKAS-seq peaks will be called by using macs2 without spKAS-seq Input. OPTIONAL."
binsizeHelp="-b [bin size]: please specify the size of bins used to identify R-loops. Default: 500."
foldchangeHelp="-f [fold change]: please specify the fold change cutoff of spKAS-seq reads difference between plus and minus strands used for R-loops identification. DEFAULT: 2."
labelsHelp="-l [labels]: please input the text file containing the labels spKAS-seq data that used to identify R-loops. Default: basename of KAS-seq files.
Example:
rep1
rep2
rep3                        ---labels.txt"
InputHelp="-n [Input]: please input the text file containing Input bed files for R-loops identification. OPTIONAL.
Example:
Input_rep1.bed
Input_rep2.bed
Input_rep3.bed              ---Input.txt"
KASseqHelp="-k [KAS-seq]: please input the text file containing spKAS-seq bed files for R-loops identification. The order and number of spKAS-seq data should be the consistent with the labels file in -l [labels]. REQUIRED.
Example:
spKAS-seq_rep1.bed
spKAS-seq_rep2.bed
spKAS-seq_rep3.bed          ---KAS-seq.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 R-loops' shell script is applied to identify R-loops from multiples pKAS-seq data."

printHelpAndExit() {
    echo -e ""
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
    echo -e "$peaksHelp"
    echo -e ""
    echo -e "$binsizeHelp"
    echo -e ""
    echo -e "$foldchangeHelp"
    echo -e ""
    echo -e "$labelsHelp"
    echo -e ""
    echo -e "$InputHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 R-loops' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

while getopts 'ht:o:s:p:b:f:l:n:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
	t) threads=$OPTARG ;;
        o) prefix=$OPTARG ;;
	s) assemblyid=$OPTARG ;;
	p) peaks=$OPTARG ;;
        b) binsize=$OPTARG ;;
	f) foldchange=$OPTARG ;;
        l) labels=$OPTARG ;;
	n) Input=$OPTARG ;;
        k) spKASseq=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $assemblyid ;then
   echo ""      
   echo "Please specify the reference genome assembly id of (sp)KAS-seq data. -s [assembly id]. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melan
ogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11."
   echo ""
   exit 1
fi

if test -z $spKASseq ;then
   echo ""
   echo "Please input the txt file containing spKAS-seq bed files for R-loops identification. REQUIRED. -k [spKAS-seq]"
   echo ""
   exit 1
fi

# get the genome size for peaks calling based on assembly id with macs2.
# genome size. e.g. human(hs): 2.7e9; mouse(mm): 1.87e9; C.elegans(ce): 9e7; fruitfly(dm): 1.2e8; rat(rn): 2.5e9; zebrafish(danRer): 1e9.
if [[ $assemblyid == "hg18" ]] || [[ $assemblyid == "hg19" ]] || [[ $assemblyid == "hg38" ]] ;then
   genomesize="2.7e9"
elif [[ $assemblyid == "mm9" ]] || [[ $assemblyid == "mm10" ]] || [[ $assemblyid == "mm39" ]] ;then     
   genomesize="1.87e9"
elif [[ $assemblyid == "ce10" ]] || [[ $assemblyid == "ce11" ]] ;then   
   genomesize="9e7"
elif [[ $assemblyid == "dm3" ]] || [[ $assemblyid == "dm6" ]] ;then     
   genomesize="1.2e8"
elif [[ $assemblyid == "rn6" ]] || [[ $assemblyid == "rn7" ]] ;then     
   genomesize="2.5e9"
elif [[ $assemblyid == "danRer10" ]] || [[ $assemblyid == "danRer11" ]] ;then   
   genomesize="1e9"
else 
   echo ""
   echo "Error: unsupported assembly id: $assemblyid. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11."
   echo ""
   exit 0
fi

# setup parameters of default options.
if test -z $prefix ;then
   prefix=$( basename ${KASseq} .txt )
fi

if test -z $binsize ;then
   binsize=500
fi

if test -z $threads ;then
   threads=1
fi

if test -z $foldchange ;then
   foldchange=2
fi

# get the number of spKAS-seq samples.
number_of_samples=$( awk 'END {print NR}' $spKASseq )

# get the labels if labels.txt is not provided.
if test -z $labels ;then
   for ((i=1; i<=${number_of_samples}; i++))
   do
   sample_selected=$(sed -n ''$i'p' $spKASseq)
   label_basename=$(basename ${sample_selected} .bed)
   echo $label_basename >> ${prefix}.labels_basename.txt
   done
   labels="${prefix}.labels_basename.txt"

else
   number_of_labels=$( awk 'END {print NR}' $labels )
   if [[ $number_of_labels != $number_of_samples ]] ;then
      echo ""
      echo "Error:the number of labels isn't consistent with the number of spKAS-seq samples!"
      echo ""
      exit -1
   fi
fi

# get the absolute path of 'KAS-pipe2 R-loop' shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

if [ $number_of_samples -eq 1 ]; then

   sample_selected=$(sed -n '1p' $spKASseq)
   
   echo "Only one spKAS-seq data: $sample_selected is provided, 'KAS-pipe2 R-loop' will not perform statistical analysis for R-loops identification."
   echo ""
   echo "Separate $sample_selected spKAS-seq mapped reads into minus and plus strands."
   echo ""
   grep ^chrM -v $sample_selected | grep '+'  | sortBed -i > ${prefix}.plus.bed
   grep ^chrM -v $sample_selected | grep '+' -v | sortBed -i > ${prefix}.minus.bed
   echo "done."
   echo ""	

   # transfer minus or plus bed filed into bedGraph files.
   echo "Converting ${prefix}.plus.bed or .minus.bed into .bg files ..."
   echo ""
   genomeCoverageBed -bg -i ${prefix}.plus.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.plus.bg
   genomeCoverageBed -bg -i ${prefix}.minus.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${prefix}.minus.bg
   echo "done."
   echo ""

   # transfer minus or plus bedGraph files into bigWig files.
   echo "Converting ${prefix}.minus.bg or .plus.bg files into .bigWig files ..."
   echo ""
   bedGraphToBigWig ${prefix}.plus.bg ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes ${prefix}.plus.bigWig
   bedGraphToBigWig ${prefix}.minus.bg ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes ${prefix}.minus.bigWig
   echo "done."
   echo ""

   # generate the $binsize bp bins with $binsize/2 bp overlap and calculate the plus or minus mapped averaged spKAS-seq reads density.
   echo "Calculate the plus or minus mapped spKAS-seq reads density on ${binsize}bp bins."
   echo ""
   distanceBetweenBins=$((-1*binsize/2))
   multiBigwigSummary bins -b ${prefix}.minus.bigWig ${prefix}.plus.bigWig --labels minus plus -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -bs $binsize -n $distanceBetweenBins -p $threads -out ${prefix}_plus_vs_minus.bins.rmbl.npz --outRawCounts ${prefix}_plus_vs_minus.bins.rmbl.tab > /dev/null 2>&1
   sed "s/nan/0/g" ${prefix}_plus_vs_minus.bins.rmbl.tab | sed "1d" | awk '{printf("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n",$1,$2,$3,$4,$5,log(($5+0.1)/($4+0.1))/log(2))}' > ${prefix}_plus_vs_minus.bins.rmbl.bed
   echo "done."
   echo ""

   if test -z $peaks ;then
      if test -z $Input ;then
         echo "spKAS-seq peaks were not provided, spKAS-seq peaks will be called for $sample_selected using macs2 ..."
         echo ""
         spKAS_basename=$(basename ${sample_selected} .bed)
         macs2 callpeak -t $sample_selected -n $spKAS_basename --broad -g $genomesize --broad-cutoff 0.01 -q 0.01
         echo "done."
         echo ""

         echo "Filter bins on spKAS-seq peaks."
         echo ""
         intersectBed -a ${prefix}_plus_vs_minus.bins.rmbl.bed -b ${spKAS_basename}_peaks.broadPeak -wa -f 0.5 | sort -u | sortBed -i > ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed
         echo "done."
         echo ""

         rm -f ${spKAS_basename}_model.r
         rm -f ${spKAS_basename}_peaks.broadPeak
         rm -f ${spKAS_basename}_peaks.gappedPeak
         rm -f ${spKAS_basename}_peaks.xls   
         
      elif test -n "$Input" ;then
         Input_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $Input)
         echo "spKAS-seq peaks were not provided, spKAS-seq peaks will be called for $sample_selected vs $Input_list using macs2 ..."
         echo ""
         spKAS_basename=$(basename ${sample_selected} .bed)
         macs2 callpeak -t $sample_selected -c $Input_list -n $spKAS_basename --broad -g $genomesize --broad-cutoff 0.01 -q 0.01
         echo "done."
         echo ""

         echo "Filter bins on spKAS-seq peaks."
         echo ""
         intersectBed -a ${prefix}_plus_vs_minus.bins.rmbl.bed -b ${spKAS_basename}_peaks.broadPeak -wa -f 0.5 | sort -u | sortBed -i > ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed
         echo "done."
         echo ""

         rm -f ${spKAS_basename}_model.r
         rm -f ${spKAS_basename}_peaks.broadPeak
         rm -f ${spKAS_basename}_peaks.gappedPeak
         rm -f ${spKAS_basename}_peaks.xls
      fi

   else 
      echo "$peaks was provided. spKAS-seq peaks will not be recalled."
      echo ""

      echo "Filter bins on spKAS-seq peaks."
      echo ""
      intersectBed -a ${prefix}_plus_vs_minus.bins.rmbl.bed -b $peaks -wa -f 0.5 | sort -u | sortBed -i > ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed
      echo "done."
      echo ""

   fi

   # identify bins with $foldchange different spKAS-seq reads density between plus and minus strands.
   echo "Identify R-loop bins on plus or minus strands with $foldchange fold spKAS-seq density change."
   echo ""
   awk -v x=$foldchange '$6>=x {printf("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n",$1,$2,$3,$4,$5,$6)}' ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed > R-loop_plus_${prefix}.bins.rmbl.bed
   awk -v x=$foldchange '$6<=-x {printf("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n",$1,$2,$3,$4,$5,$6)}' ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed > R-loop_minus_${prefix}.bins.rmbl.bed
   echo "done."
   echo ""

   # merge R-loops bins into R-loops peaks.
   echo "Merge R-loop bins into R-loop enriched regions."
   echo ""
   sortBed -i R-loop_plus_${prefix}.bins.rmbl.bed | mergeBed -i - | awk '{printf("%s\t%d\t%d\t%s\n",$1,$2,$3,"+")}' > ${prefix}_R-loops.plus.bed
   sortBed -i R-loop_minus_${prefix}.bins.rmbl.bed | mergeBed -i - | awk '{printf("%s\t%d\t%d\t%s\n",$1,$2,$3,"-")}' > ${prefix}_R-loops.minus.bed
   echo "done."
   echo ""

   echo "Generate R-loops density at 50bp resolution."
   echo ""
   multiBigwigSummary bins -b ${prefix}.minus.bigWig ${prefix}.plus.bigWig --labels minus plus -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -bs 50 -p $threads -out ${prefix}_plus_vs_minus.50bp.bins.rmbl.npz --outRawCounts ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab > /dev/null 2>&1

   sed "s/nan/0/g" ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab | sed "1d" | awk '$5-$4>0 {printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$5-$4)}' > ${prefix}_R-loop.density.plus.bg
   sed "s/nan/0/g" ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab | sed "1d" | awk '$4-$5>0 {printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$5-$4)}' > ${prefix}_R-loop.density.minus.bg
   cat ${prefix}_R-loop.density.plus.bg ${prefix}_R-loop.density.minus.bg | sortBed -i > ${prefix}_R-loop.density.bg
   echo "done."
   echo ""
#  bedGraphToBigWig ${prefix}_R-loop.density.bg ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes ${prefix}_R-loop.density.bigWig

   echo "Clean up intermediate files."
   echo ""
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
   rm -f ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed
   rm -f R-loop_plus_${prefix}.bins.rmbl.bed
   rm -f R-loop_minus_${prefix}.bins.rmbl.bed
   # clean up intermediate files generated during R-loops density.
   rm -f ${prefix}_plus_vs_minus.50bp.bins.rmbl.npz
   rm -f ${prefix}_plus_vs_minus.50bp.bins.rmbl.tab
   rm -f ${prefix}.labels_basename.txt
   echo "done."
   echo ""

   # move R-loops and R-loops signal density (bigWig file) into R-loops directory.
   mkdir -p R-loops
   cd R-loops
   mv ../${prefix}_R-loops.plus.bed ./
   mv ../${prefix}_R-loops.minus.bed ./
   mv ../${prefix}_R-loop.density.plus.bg ./
   mv ../${prefix}_R-loop.density.minus.bg ./
   mv ../${prefix}_R-loop.density.bg ./
   # mv ../${prefix}/${prefix}_R-loop.density.bigWig ./
   cd ..
   echo "Move R-loops identification output files into 'R-loops' folder. done."
   echo ""

elif [ $number_of_samples -gt 1 ]; then

   echo "$number_of_samples spKAS-seq data was provided, 'KAS-pipe2 R-loop' will perform statistical analysis for R-loops identification."
   echo ""
   
   # generate the annotation file with header.
   echo -e "\tcondition" > annotation.header.txt

   for ((i=1; i<=${number_of_samples}; i++))
   do
   spKAS_selected=$( sed -n ''$i'p' $spKASseq )
   label_selected=$( sed -n ''$i'p' $labels )
   echo "Separate ${spKAS_selected} into minus and plus strands."
   echo ""
   grep ^chrM -v $spKAS_selected | grep + > ${label_selected}.plus.bed
   grep ^chrM -v $spKAS_selected | grep + -v > ${label_selected}.minus.bed
   echo "done."
   echo ""

   # convert minus or plus bed files into bedGraph files.
   echo "Convert ${label_selected}.plus.bed or .minus.bed into .bg files."
   echo ""
   genomeCoverageBed -bg -i ${label_selected}.plus.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${label_selected}.plus.bg
   genomeCoverageBed -bg -i ${label_selected}.minus.bed -g ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes > ${label_selected}.minus.bg
   echo "done."
   echo ""

   # convert minus or plus bedGraph files into bigWig files.
   echo "Convert ${label_selected}.plus.bg or .minus.bg files into .bigWig files."
   echo ""
   bedGraphToBigWig ${label_selected}.plus.bg ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes ${label_selected}.plus.bigWig
   bedGraphToBigWig ${label_selected}.minus.bg ${SH_SCRIPT_DIR}/../chrom_size/${assemblyid}.chrom.sizes ${label_selected}.minus.bigWig
   echo "done."
   echo ""

   echo "Generate the R-loop density at 50bp resolution ..."
   echo ""
   multiBigwigSummary bins -b ${label_selected}.minus.bigWig ${label_selected}.plus.bigWig --labels minus plus -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -bs 50 -p $threads -out ${label_selected}_plus_vs_minus.50bp.bins.rmbl.npz --outRawCounts ${label_selected}_plus_vs_minus.50bp.bins.rmbl.tab > /dev/null 2>&1
   sed "s/nan/0/g" ${label_selected}_plus_vs_minus.50bp.bins.rmbl.tab | sed "1d" | awk '$5-$4>0 {printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$5-$4)}' > ${label_selected}_R-loop.density.plus.bg
   sed "s/nan/0/g" ${label_selected}_plus_vs_minus.50bp.bins.rmbl.tab | sed "1d" | awk '$4-$5>0 {printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$5-$4)}' > ${label_selected}_R-loop.density.minus.bg
   cat ${label_selected}_R-loop.density.plus.bg ${label_selected}_R-loop.density.minus.bg | sortBed -i > ${label_selected}_R-loop.density.bg
   echo "done."
   echo ""

   echo "Clean up ${spKAS_selected} intermediate files."
   echo ""
   rm -f ${label_selected}.plus.bed
   rm -f ${label_selected}.minus.bed
   rm -f ${label_selected}.plus.bg
   rm -f ${label_selected}.minus.bg
   rm -f ${label_selected}_plus_vs_minus.50bp.bins.rmbl.npz
   rm -f ${label_selected}_plus_vs_minus.50bp.bins.rmbl.tab
   echo "done."
   echo ""

   echo "${label_selected}.plus.bigWig" >> ${prefix}.spKAS-seq.plus.txt
   echo "${label_selected}.minus.bigWig" >> ${prefix}.spKAS-seq.minus.txt

   echo "${label_selected}.plus" >> ${prefix}.labels.plus.txt
   echo "${label_selected}.minus" >> ${prefix}.labels.minus.txt
    
   # input the spKAS-seq labels into annatation file.
   echo -e "${label_selected}.plus\tplus" >> annotation.plus.txt
   echo -e "${label_selected}.minus\tminus" >> annotation.minus.txt

   mkdir -p R-loops
   cd R-loops
   mv ../${label_selected}_R-loop.density.plus.bg ./
   mv ../${label_selected}_R-loop.density.minus.bg ./
   mv ../${label_selected}_R-loop.density.bg ./
   cd ..
   echo "Move ${spKAS_selected} R-loops density output files into 'R-loops' folder. done."
   echo ""

   done

   # generate the $binsize bp bins with $binsize/2 bp overlap and calculate the plus or minus mapped averaged spKAS-seq reads density.
   spKAS_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $spKASseq)
   spKAS_plus_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' ${prefix}.spKAS-seq.plus.txt)
   spKAS_minus_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' ${prefix}.spKAS-seq.minus.txt)
   labels_plus_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' ${prefix}.labels.plus.txt)
   labels_minus_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' ${prefix}.labels.minus.txt)
   
   # generate the labels files.
   cat ${prefix}.labels.minus.txt ${prefix}.labels.plus.txt > labels.plus_vs_minus.txt 
   cat annotation.header.txt annotation.minus.txt annotation.plus.txt > annotation.txt

   distanceBetweenBins=$((-1*binsize/2))

   echo "Calculate the plus or minus mapped spKAS-seq density on ${binsize}bp bins with ${distanceBetweenBins}bp overlap."
   echo ""
   multiBigwigSummary bins -b ${spKAS_minus_list} ${spKAS_plus_list} --labels ${labels_minus_list} ${labels_plus_list} -bl ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -bs $binsize -n $distanceBetweenBins -p $threads -out ${prefix}_plus_vs_minus.bins.rmbl.npz --outRawCounts ${prefix}_plus_vs_minus.bins.rmbl.tab > /dev/null 2>&1

   sed "s/nan/0/g" ${prefix}_plus_vs_minus.bins.rmbl.tab | sed "1d" > ${prefix}_plus_vs_minus.bins.rmbl.bed
   echo "done."
   echo ""
    
   if test -z $peaks ;then
      
      if test -z $Input ;then 
      echo "spKAS-seq peaks were not provided, spKAS-seq peaks will be called for ${spKAS_list} using macs2 ..."
      echo ""
      macs2 callpeak -t ${spKAS_list} -n ${prefix} --broad -g $genomesize --broad-cutoff 0.01 -q 0.01
      echo "done."
      echo ""

      # filter bins overlap with spKAS-seq peaks.
      echo "Filter ${binsize}bp bins on spKAS-seq peaks."
      echo ""
      intersectBed -a ${prefix}_plus_vs_minus.bins.rmbl.bed -b ${prefix}_peaks.broadPeak -wa -f 0.5 | sort -u | sortBed -i > ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed
      echo "done."
      echo ""

      rm -f ${prefix}_model.r
      rm -f ${prefix}_peaks.broadPeak
      rm -f ${prefix}_peaks.gappedPeak
      rm -f ${prefix}_peaks.xls
   
      elif test -n "$Input" ;then
      Input_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $Input)
      echo "spKAS-seq peaks were not provided, spKAS-seq peaks will be called for ${spKAS_list} vs ${Input_list} using macs2 ..."
      echo ""
      macs2 callpeak -t ${spKAS_list} -c ${Input_list} -n ${prefix} --broad -g $genomesize --broad-cutoff 0.01 -q 0.01
      echo "done."
      echo ""

      # filter bins overlap with spKAS-seq peaks.
      echo "Filter ${binsize}bp bins on spKAS-seq peaks."
      echo ""
      intersectBed -a ${prefix}_plus_vs_minus.bins.rmbl.bed -b ${prefix}_peaks.broadPeak -wa -f 0.5 | sort -u | sortBed -i > ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed
      echo "done."
      echo ""

      rm -f ${prefix}_model.r
      rm -f ${prefix}_peaks.broadPeak
      rm -f ${prefix}_peaks.gappedPeak
      rm -f ${prefix}_peaks.xls

      fi
   	 
   else
      echo "$peaks was provided. spKAS-seq peaks will not be recalled."
      echo ""

      # filter bins overlap with spKAS-seq peaks.
      echo "Filter ${binsize}bp bins on spKAS-seq peaks."
      echo ""
      intersectBed -a ${prefix}_plus_vs_minus.bins.rmbl.bed -b $peaks -wa -f 0.5 | sort -u | sortBed -i > ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed
      echo "done."
      echo ""
   fi

   echo "Generate plus vs minus spKAS-seq matrix ..."
   echo ""
   cut -f1,2,3 --complement ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed | awk '{for(i=1;i<=NF;i++){printf "%d\t", $i*10}; printf "\n"}' > ${prefix}_plus_vs_minus.bins.rmbl.peaks.matrix
   awk '{printf("%s\n",$1"-"$2"-"$3)}' ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed > ${prefix}_plus_vs_minus.bins.rmbl.peaks.rowname
   paste ${prefix}_plus_vs_minus.bins.rmbl.peaks.rowname ${prefix}_plus_vs_minus.bins.rmbl.peaks.matrix > ${prefix}_plus_vs_minus.bins.rmbl.peaks.without_header.txt

   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf "\t" a[i,j];print ""}}' labels.plus_vs_minus.txt > header.txt

   cat header.txt ${prefix}_plus_vs_minus.bins.rmbl.peaks.without_header.txt > ${prefix}_plus_vs_minus.bins.rmbl.peaks.txt
   echo "done."
   echo ""

   echo "Clean up the intermediate files."
   rm -f ${prefix}.spKAS-seq.plus.txt
   rm -f ${prefix}.spKAS-seq.minus.txt
   rm -f ${prefix}.labels.plus.txt
   rm -f ${prefix}.labels.minus.txt
   rm -f ${prefix}.labels_basename.txt
   rm -f *plus.bigWig
   rm -f *minus.bigWig
   rm -f ${prefix}_plus_vs_minus.bins.rmbl.npz
   rm -f ${prefix}_plus_vs_minus.bins.rmbl.tab
   rm -f ${prefix}_plus_vs_minus.bins.rmbl.bed 
   rm -f ${prefix}_plus_vs_minus.bins.rmbl.peaks.bed
   rm -f ${prefix}_plus_vs_minus.bins.rmbl.peaks.matrix
   rm -f ${prefix}_plus_vs_minus.bins.rmbl.peaks.rowname
   rm -f ${prefix}_plus_vs_minus.bins.rmbl.peaks.without_header.txt
   rm -f header.txt
   rm -f annotation.header.txt 
   rm -f annotation.minus.txt 
   rm -f annotation.plus.txt
   echo "done."
   echo ""


   echo "Perform differential spKAS-seq analysis between Watson and Crick strands with DEseq2, and identify genome-wide R-loops."
   echo ""
   Rscript --vanilla ${SH_SCRIPT_DIR}/../R/Identify_R-loops.R ${prefix}_plus_vs_minus.bins.rmbl.peaks.txt annotation.txt
   echo "done."
   echo ""

   sed "s/\"//g" spKAS-seq_plus_vs_minus_DESeq2_R-loops_output.csv | sed "s/\,/\t/g" > ${prefix}_DESeq2_R-loops_output.txt
   sed "s/\"//g" spKAS-seq_plus_vs_minus_DESeq2_R-loops_Fold1.5_padj0.05_output.csv | sed "s/\,/\t/g" > ${prefix}_DESeq2_R-loops_Fold1.5_padj0.05_output.txt

   rm -f spKAS-seq_plus_vs_minus_DESeq2_R-loops_output.csv
   rm -f spKAS-seq_plus_vs_minus_DESeq2_R-loops_Fold1.5_padj0.05_output.csv
#  rm -f ${prefix}_plus_vs_minus.bins.rmbl.peaks.txt
   rm -f annotation.txt
   rm -f labels.plus_vs_minus.txt

fi

echo "'KAS-pipe2 R-loop' run successfully!"
