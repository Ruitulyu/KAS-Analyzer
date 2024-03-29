#!/bin/bash
# 'KAS-Analyzer statistics' was developed by Ruitu Lyu on 12-21-2021.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-Analyzer statistics [ -h/--help ] [ -o prefix ] [ -l labels ] [ -s summary folder ]"
exampleHelp="Example: nohup KAS-Analyzer statistics -o KAS-seq_statistics -l labels.txt -s summary.txt &"
prefixHelp="-o [prefix]: please input the prefix (basename) of 'KAS-Analyzer statistics' output files. REQUIRED."
labelsHelp="-l [labels]: please input the txt file containing labels of (sp)KAS-seq summary file generated by 'KAS-Analyzer (sp)KAS-seq'. DEFAULT: basename of summary file.
KAS-seq.rep1
KAS-seq.rep2
KAS-seq.rep3
KAS-seq.rep4        ---labels.txt"
summaryHelp="-s [summary folder or summary.txt]: please input the absolute path of folder with summary files or the txt file containing the summary files. e.g. -s summary.txt or /absolute path of summary folder. REQUIRED.
WT.rep1.KAS-seq_mapping_summary.txt                     or  /absolute path/Summary/
WT.rep2.KAS-seq_mapping_summary.txt
KO.rep1.KAS-seq_mapping_summary.txt
KO.rep2.KAS-seq_mapping_summary.txt     ---summary.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-Analyzer statistics' shell script is applied to generate the table containing (sp)KAS-seq mapping statistics."

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
    echo -e "$summaryHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer statistics' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit 0
fi

# get the value of options.
while getopts 'ho:l:s:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        o) prefix=$OPTARG ;;
        l) labels=$OPTARG ;;
        s) summary=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $prefix ;then
   echo ""
   echo "Please specify the prefix (basename) of 'KAS-Analyzer statistics' output files. -o [prefix]"
   echo ""
   exit -1
fi

if test -z $summary ;then
   echo ""
   echo "Please specify the absolute path of folder with summary files or the txt file containing the summary files. e.g. -s summary.txt or /absolute path/Summary."
   echo ""
   exit -1
fi

# setup parameters of default options.
if test -z $labels ;then

   if [[ ${summary##*.} == txt ]] ;then
      number_of_summary=$( awk 'END {print NR}' $summary )

      for ((i=1; i<=${number_of_summary}; i++))
      do
      summary_selected=$( sed -n ''$i'p' $summary )
      summary_basename=$( basename ${summary_selected} .txt )
      echo $summary_basename >> ${prefix}.summary_basename.txt
      done

      labels="${prefix}.summary_basename.txt"

   else 

      cd $summary	    

      for file in ./*mapping_summary.txt
      do
      summary_basename=$( basename $file .txt )
      echo $summary_basename >> ${prefix}.summary_basename.txt
      echo $file >> ${prefix}.summary_files.txt
      done

      number_of_summary=$(awk 'END {print NR}' ${prefix}.summary_files.txt )     
      labels="${prefix}.summary_basename.txt"
      summary="${prefix}.summary_files.txt"
   fi

else 
   if [[ ${summary##*.} == txt ]] ;then
      # the the numebr of summary files.	   
      number_of_summary=$(awk 'END {print NR}' $summary )
   else 
     
      cd $summary

      for file in ./*mapping_summary.txt
      do
      echo $file >> ${prefix}.summary_files.txt
      number_of_summary=$(awk 'END {print NR}' ${prefix}.summary_files.txt )
      summary="${prefix}.summary_files.txt"    
      done
    
   fi
fi

number_of_summary=$( awk 'END {print NR}' $summary )

# Input the header of output summary file.
echo -e "Samples\tClean_reads\tMapped_reads\tDeduplicated_reads\tMapping_ratios\tDuplication_ratios" > ${prefix}_summary.txt

for ((j=1; j<=${number_of_summary}; j++))
do
Summary_selected=$( sed -n ''$j'p' $summary )
Basename_selected=$( sed -n ''$j'p' $labels )

Clean_reads=$( grep "KAS-seq reads:" $Summary_selected | awk '{print $5}' )
Mapped_reads=$( grep "Number of mapped reads." $Summary_selected | awk '{print $6}' )
Deduplicated_reads=$( grep "Number of deduplicated mapped reads" $Summary_selected | awk '{print $6}' )
Mapping_ratios=$( grep "Mapping ratios:" $Summary_selected | awk '{print $3}' )
Duplication_ratios=$( grep "Duplication ratios." $Summary_selected | awk '{print $4}' )

echo -e "$Basename_selected\t$Clean_reads\t$Mapped_reads\t$Deduplicated_reads\t$Mapping_ratios\t$Duplication_ratios" >> ${prefix}_summary.txt
done

rm -f ${prefix}.summary_basename.txt
rm -f ${prefix}.summary_files.txt

echo "'KAS-Analyzer statistics' run successfully!"
