#!/bin/bash
# 'KAS-pipe2 correlation' was developed by Ruitu Lyu on 12-10-2021.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-pipe2 correlation [ -h/--help ] [ -m correlation method ] [ -t threads ] [ -s assembly id ] [ -r regions ] [ -f peaks file ] [ -p plot types ] [ -o prefix ] [ -l labels ] [ -k KAS-seq ]"
exampleHelp="Example:
On peaks:             
nohup KAS-pipe2 correlation -m pearson -t 10 -s hg19 -r peaks -f KAS-seq_peaks.bed -p heatmap -o KAS-seq -l labels.txt -k KAS-seq.txt &
On bins:
nohup KAS-pipe2 correlation -m pearson -t 10 -s hg19 -r bin -p heatmap -o KAS-seq -l labels.txt -k KAS-seq.txt &"
methodsHelp="-m [correlation method]: please specify the methods to calculate correlation coefficients. e.g. pearson, kendall or spearman. DEFAULT: pearson."
threadsHelp="-t [threads]: please specify the number of threads. DEFAULT: 1."
assemblyidHelp="-s [assembly id]: please specify the genome assembly id of your (sp)KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
regionsHelp="-r [regions]: please specify the region types to calculate the (sp)KAS-seq density matrix. e.g. bin or peak. DEFAULT: bin."
peakfileHelp="-f [peaks file]: please input the merged peaks list file. Note: only valid when '-r peaks' is specified. REQUIRED." 
plotsHelp="-p [plot types]: please specify the plot types to generate correlation plot. DEFAULT: scatterplot."
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 correlation' output files. REQUIRED."
labelsHelp="-l [labels]: please input the text file containing the labels of (sp)KAS-seq data that show in heatmap or scatterplot. REQUIRED.
Example:
WT_rep1
WT_rep2
WT_rep3
KO_rep1
KO_rep2   
KO_rep2                       ---labels.txt"
KASseqHelp="-k [KAS-seq]: please input the text file containing the indexed bam files of (sp)KAS-seq data that used to calculate correlation coefficient and generate correlation plots. REQUIRED.
Example:
KAS-seq_WT.rep1.bigWig
KAS-seq_WT.rep2.bigWig
KAS-seq_WT.rep3.bigWig
KAS-seq_KO.rep1.bigWig
KAS-seq_KO.rep2.bigWig
KAS-seq_KO.rep3.bigWig        ---KAS-seq_data.txt"
helpHelp="-h: print this help and exit.
Note: The 'KAS-pipe2 correlation' shell script is applied to calculate correlation coefficients and generate correlation plots between replicates or (sp)KAS-seq data of different conditions."

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$methodsHelp"
    echo -e ""
    echo -e "$threadsHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$regionsHelp"
    echo -e ""
    echo -e "$peakfileHelp"
    echo -e ""
    echo -e "$plotsHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$labelsHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 correlation' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'hm:t:s:r:f:p:o:l:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        m) methods=$OPTARG ;;
        t) threads=$OPTARG ;;
        s) assemblyid=$OPTARG ;;
        r) regions=$OPTARG ;;
	f) peakfile=$OPTARG ;;
        p) plots=$OPTARG ;;
        o) prefix=$OPTARG ;;
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
    echo "Error: unsupported assembly id: $assemblyid "
    echo ""
    exit 1
fi

if test -z $prefix ;then
   echo ""
   echo "Please specify the prefix (basename) of 'KAS-pipe2 correlation' output files. -o [prefix]"
   echo ""
   exit 1
fi

if test -z $labels ;then
   echo ""
   echo "Please input the text file containing the labels of (sp)KAS-seq data that show in heatmap or scatterplot. -l [labels]"
   echo ""
   exit 1
fi

if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing the indexed bam files of (sp)KAS-seq data to calculate correlation coefficients and generate correlation plots. -k [KAS-seq]"
   echo ""
   exit 1
fi

# test if peaks file is provided by user.
if [[ $regions == "peak" ]] && test -z $peakfile ;then
   echo ""
   echo "Please input the peaks file if you want to generated generate density matrix on KAS-seq peaks. Note: only valid when '-r peaks' is specified. -f [peaks.bed]"
   echo ""
   exit 1
fi
   	

# setup the default options.
if test -z $methods ;then
   methods="pearson"
fi

if [[ $methods != "pearson" ]] && [[ $methods != "kendall" ]] && [[ $methods != "spearman" ]] ;then
   echo ""
   echo "Error: unsupported correlation methods: $methods "
   echo ""
   exit -1
fi 

if test -z $threads ;then
   threads=1
fi

if test -z $regions ;then
   regions="bin"
fi

if test -z $plots ;then
   plots="scatterplot"
fi

# test if the number of labels is same to the number of samples.
number_of_samples=$(awk 'END {print NR}' $KASseq )
number_of_labels=$(awk 'END {print NR}' $labels )

# test if number of labels is consistent with number of samples.
if [[ ${number_of_labels} != ${number_of_samples} ]]
then
   echo ""
   echo "Error: the number of labels isn't consistent with the number of samples."
   echo ""
   exit -1
fi

# supported regions types.
if [[ $regions != "bin" ]] && [[ $regions != "peak" ]]
   then
   echo ""
   echo "Error: unsupported types of regions: $regions. e.g. bin or peak."
   echo ""
   exit -1
fi

# supported plots types.
if [[ $plots != "scatterplot" ]] && [[ $plots != "heatmap" ]]
   then
   echo	""   
   echo "Error: unsupported types of plot: $plots. e.g. scatterplot or heatmap."
   echo ""
   exit -1
fi

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# get the list of KAS-seq samples and labels.
KASseq_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $KASseq)
labels_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $labels)

# generate the matrix of KAS-seq data on 1kb bins
if [[ $regions == "bin" ]]; then

   echo "Generate the matrix of KAS-seq data on 1kb bins ..."
   echo ""
   multiBigwigSummary bins -b $KASseq_list --labels $labels_list --binSize 1000 -p $threads --blackListFileName ${SH_SCRIPT_DIR}/../blacklist/${assemblyid}-blacklist.bed -out ${prefix}_on_${regions}.npz --outRawCounts ${prefix}_on_${regions}.tab
   echo "done."
   echo ""

   sed "s/nan/0/g" ${prefix}_on_${regions}.tab | sed "1d" > ${prefix}_on_${regions}.bed 

   # calculate the average value of every single row in the table.
   echo "Calculate the average KAS-seq density on 1kb bins ..."
   echo ""
   awk 'BEGIN{if(NR>0) a[NR]=0}{if(NR>0) for(i=4; i<=NF; i++) a[NR]+=$i}END{for(j in a) print a[j]/(NF-3) }' ${prefix}_on_${regions}.bed > ${prefix}_on_${regions}.average
   echo "done."
   echo ""

   # filter the bins with averaged KAS-seq RPKM lower than 2.
   echo "Filter matrix of KAS-seq density on 1kb bins ..."
   echo ""
   paste ${prefix}_on_${regions}.average ${prefix}_on_${regions}.bed | sort -k 1 -n -r | sed '1,100d' | awk '$1>=5 {print $0}' > ${prefix}_on_${regions}.filter.bed
   cut -f1,2,3,4 --complement ${prefix}_on_${regions}.filter.bed | awk '{for(i=1;i<=NF;i++){printf "%.2f\t", log($i+1)/log(2)}; printf "\n"}' > ${prefix}_on_${regions}.matrix
   awk '{printf("%s\n",$2"-"$3"-"$4)}' ${prefix}_on_${regions}.filter.bed > ${prefix}_on_${regions}.rowname
   echo "done."
   echo ""

   echo "Generate the final matrix of KAS-seq density on 1kb bins ..."
   echo ""
   paste ${prefix}_on_${regions}.rowname ${prefix}_on_${regions}.matrix > ${prefix}_on_${regions}.without_header.txt
   echo -e "" > ${prefix}.header1.txt
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels > ${prefix}.header2.txt
   paste ${prefix}.header1.txt ${prefix}.header2.txt > ${prefix}.header.txt
   cat ${prefix}.header.txt ${prefix}_on_${regions}.without_header.txt > ${prefix}_on_${regions}.txt
   echo "done."
   echo ""

   # clean up the intermediate files.
   rm -f ${prefix}_on_${regions}.npz
   rm -f ${prefix}_on_${regions}.tab
   rm -f ${prefix}_on_${regions}.bed
   rm -f ${prefix}_on_${regions}.average
   rm -f ${prefix}_on_${regions}.filter.bed
   rm -f ${prefix}_on_${regions}.matrix
   rm -f ${prefix}_on_${regions}.rowname
   rm -f ${prefix}_on_${regions}.without_header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}.header.txt
   echo "Intermediate files clean up, done."

   # generate the matrix of KAS-seq data on peaks
elif [[ $regions == "peak" ]]; then
   echo "Generate the matrix of KAS-seq data on $peakfile ..."
   echo ""
   peakfile_basename=$(basename ${peakfile} .bed)
   multiBigwigSummary BED-file --bwfiles $KASseq_list --BED $peakfile --labels $labels_list -p $threads -out ${prefix}_on_${peakfile_basename}.npz --outRawCounts ${prefix}_on_${peakfile_basename}.tab
   sed "s/nan/0/g" ${prefix}_on_${peakfile_basename}.tab | sed "1d" > ${prefix}_on_${peakfile_basename}.bed
   echo "done."
   echo ""
   
   echo "Generate the final matrix of (sp)KAS-seq density on $peakfile ..."
   echo ""
   cut -f1,2,3 --complement ${prefix}_on_${peakfile_basename}.bed | awk '{for(i=1;i<=NF;i++){printf "%.2f\t", $i}; printf "\n"}' > ${prefix}_on_${peakfile_basename}.matrix
   awk '{printf("%s\n",$1"-"$2"-"$3)}' ${prefix}_on_${peakfile_basename}.bed > ${prefix}_on_${peakfile_basename}.rowname
   paste ${prefix}_on_${peakfile_basename}.rowname ${prefix}_on_${peakfile_basename}.matrix > ${prefix}_on_${regions}.without_header.txt
   echo -e "" > ${prefix}.header1.txt
   awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] "\t";print ""}}' $labels > ${prefix}.header2.txt
   paste ${prefix}.header1.txt ${prefix}.header2.txt > ${prefix}.header.txt
   cat ${prefix}.header.txt ${prefix}_on_${regions}.without_header.txt > ${prefix}_on_${regions}.txt
   echo "done."
   echo ""

   rm -f ${prefix}_on_${peakfile_basename}.npz
   rm -f ${prefix}_on_${peakfile_basename}.tab
   rm -f ${prefix}_on_${peakfile_basename}.bed
   rm -f ${prefix}_on_${peakfile_basename}.matrix
   rm -f ${prefix}_on_${peakfile_basename}.rowname
   rm -f ${prefix}_on_${regions}.without_header.txt
   rm -f ${prefix}.header1.txt
   rm -f ${prefix}.header2.txt
   rm -f ${prefix}.header.txt
   echo "Intermediate files clean up done."
   echo ""

fi

if [[ $plots == "heatmap" ]]; then
   echo "Plot correlation heatmap ..."
   echo ""
   Rscript --vanilla ${SH_SCRIPT_DIR}/../R/Plotcorr_heatmap.R ${prefix}_on_${regions}.txt

   mv KAS-seq_corr_heatmap.png ${prefix}_KAS-seq_corr_heatmap.png
   mv KAS-seq_corr_heatmap.svg ${prefix}_KAS-seq_corr_heatmap.svg
   mv KAS-seq_corr_circle.png ${prefix}_KAS-seq_corr_circle.png
   mv KAS-seq_corr_circle.svg ${prefix}_KAS-seq_corr_circle.svg
   mv KAS-seq_correlation.csv ${prefix}_KAS-seq_correlation.csv
   mv KAS-seq_correlation_pvalue.csv ${prefix}_KAS-seq_correlation_pvalue.csv
   echo "done."
   echo ""

elif [[ $plots == "scatterplot" ]]; then
   if [ $number_of_samples -eq 2 ]; then
   
   echo "Plot correlation scatterplot using ${methods} method ..."
   echo ""   
   Rscript --vanilla ${SH_SCRIPT_DIR}/../R/Plotcorr_scatterplot.R ${prefix}_on_${regions}.txt ${methods}

   mv KAS-seq_corr_scatterplot.png ${prefix}_KAS-seq_corr_scatterplot.png
   mv KAS-seq_corr_scatterplot.svg ${prefix}_KAS-seq_corr_scatterplot.svg
   echo "done."
   echo ""
   
   elif [ $number_of_samples -gt 2 ]; then
   echo "Plot correlation scatter matrix ..."
   echo ""   
   Rscript --vanilla ${SH_SCRIPT_DIR}/../R/Plotcorr_scatterplot_matrix.R ${prefix}_on_${regions}.txt
   
   mv KAS-seq_corr_scatterplot_matrix.png ${prefix}_KAS-seq_corr_scatterplot_matrix.png
   mv KAS-seq_corr_scatterplot_matrix.svg ${prefix}_KAS-seq_corr_scatterplot_matrix.svg
   echo "done."
   echo ""

   elif [ $number_of_samples -eq 1 ]; then
   echo ""
   echo "please make sure the number of samples for correlation analysis must be over 2."
   echo ""
   printHelpAndExit
   fi
fi

echo "'KAS-pipe2 correlation' run successfully!"
