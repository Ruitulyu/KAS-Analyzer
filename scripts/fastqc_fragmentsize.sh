#!/bin/bash
# 'KAS-pipe2 fragmentsize' was developed by Ruitu Lyu on 12-14-2021.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 fragmentsize [ -h/--help ] [ -o prefix ] [ -l labels ] [ -k KAS-seq ]"
exampleHelp="Example: nohup KAS-pipe2 fragmentsize -o KAS-seq_fragmentsize -l labels.txt -k KAS-seq.txt &"
prefixHelp="-o [prefix]: please specify the prefix (basename) of 'KAS-pipe2 fragmentsize' output files. REQUIRED."
labelsHelp="-l [labels.txt]: please input the text file containing the labels of paired-end KAS-seq or spKAS-seq data that show in fragment size plot. Default: basename of (sp)KAS-seq bed files.
Example:
WT.rep1
WT.rep2
WT.rep3
WT.rep4                        ---labels.txt"
KASseqHelp="-k [KAS-seq.txt]: please input the text file containing bed files (uniquely mapped reads from 'KAS-pipe2 (sp)KAS-seq'), which are used to calcuate fragment size of DNA fragments. The order and number of (sp)KAS-seq bed files should be the consistent with the labels in labels.txt file. REQUIRED.
Example:
KAS-seq_WT_PE.rep1.bed
KAS-seq_WT_PE.rep2.bed
KAS-seq_WT_PE.rep3.bed
KAS-seq_WT_PE.rep4.bed         ---KAS-seq.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 fragmentsize' shell script is applied to calculate and plot fragment size of (sp)KAS-seq data. Note: it only works for paired-end (sp)KAS-seq data."

printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
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

# if no parameters was provided, 'KAS-pipe2 fragmentsize' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'ho:l:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        o) prefix=$OPTARG ;;
        l) labels=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $prefix ;then
   echo ""
   echo "Please specify the prefix (basename) of 'KAS-pipe2 fragmentsize' output files. -o [prefix]"
   echo ""
   exit -1
fi

if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing bed files (uniquely mapped reads from 'KAS-pipe2 (sp)KAS-seq'), which are used to calcuate fragment size of DNA fragments. -k [KAS-seq]"
   echo ""
   exit -1
fi

# get the path of shell scripts.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# get the number of samples.
number_of_samples=$(awk 'END {print NR}' $KASseq )

# setup the default parameter of $labels.
if test -n "$labels" ;then
   number_of_labels=$(awk 'END {print NR}' $labels )
   if [[ $number_of_labels != $number_of_samples ]] ;then 
      echo "" 
      echo "Error: the number of labels isn't consistent with the number of samples." 
      echo ""
      exit -1   
   fi
fi

# if labels.txt is not provided, use the basename of KAS-seq bed files.
if test -z $labels ;then
for ((i=1; i<=${number_of_samples}; i++))
do
sample_selected=$(sed -n ''$i'p' $KASseq)
label_basename=$(basename ${sample_selected} .bed)
echo $label_basename >> ${prefix}.labels_basename.txt
done
labels = "${prefix}.labels_basename.txt"
fi

# generate the fragment size txt file.
echo -e "Samples\tLength" > ${prefix}_fragmentsize.txt

for ((i=1; i<=${number_of_samples}; i++))       
do
sample_selected=$(sed -n ''$i'p' $KASseq)
label_selected=$(sed -n ''$i'p' $labels)
echo "Calculate the size of DNA fragments of $sample_selected ..."
echo ""
awk '{printf("%d\n",$3-$2)}' $sample_selected | shuf | head -n 200000 | awk -v x=$label_selected '{printf("%s\t%d\n",x,$1)}' >> ${prefix}_fragmentsize.txt
echo "done."
echo ""
done

rm -f ${prefix}.labels_basename.txt

echo "Plot the size of DNA fragments of paired-end KAS-seq data."
Rscript --vanilla ${SH_SCRIPT_DIR}/../R/Fragmentsize_plot.R ${prefix}_fragmentsize.txt
echo "done."
echo ""

mv KAS-seq_fragment_size_density_plot.png ${prefix}_fragment_size_density_plot.png
mv KAS-seq_fragment_size_density_plot.svg ${prefix}_fragment_size_density_plot.svg

# rm -f ${prefix}_fragmentsize.txt

echo "'KAS-pipe2 fragmentsize' run successfully!"
