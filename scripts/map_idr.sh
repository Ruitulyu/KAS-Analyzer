#!/bin/bash
  
# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-Analyzer idr [ -h ] [ -1 rep1 peaklist ] [ -2 rep2 peaklist ] [ -c peak caller] [ -t file type ] [ -r rank ] [ -o output file ]"
exampleHelp="Example:
     narrowPeak (macs2):
     nohup KAS-Analyzer idr -1 KAS-seq.narrow_peakslist.rep1.bed -2 KAS-seq.narrow_peakslist.rep2.bed -c macs2 -t narrowPeak -r signal.value -o KAS-seq.narrow_peakslist &

     broadPeak (epic2):
     nohup KAS-Analyzer idr -1 KAS-seq.broad_peakslist.rep1.bed -2 KAS-seq.broad_peakslist.rep2.bed -c epic2 -t broadPeak -r signal.value -o KAS-seq.broad_peakslist &"
rep1peakslistHelp="-1 [rep1 peaklist]: please provide the KAS-seq peaks file of the first replicate. REQUIRED"
rep2peakslistHelp="-2 [rep2 peaklist]: please provide the KAS-seq peaks file of the second replicate. REQUIRED"
peakscallerHelp="-c [peaks caller]: please specify the peaks caller used to generate the peak list. macs2 or epic2. REQUIRED"
filetypeHelp="-t: please specify the output file type, narrowPeak or broadPeak. DEFAULT: broadPeak"
rankHelp="-r [rank]: please specify the type of signals to rank the peaks in your provided lists. signal.value, p.value or q.value. DEFAULT: signal.value"
outputHelp="-o [output file]: please specify the file where the output should be written. REQUIRED"
helpHelp="-h: print this help and exit.
Note: The 'KAS-Analyzer idr' shell script mainly invoke the idr (Irreproducible Discovery Rate) framework is a uniﬁed approach to measure the reproducibility of ﬁndings identiﬁed from KAS-seq replicate experiments and provide highly stable thresholds based on reproducibility, please refer to https://github.com/nboley/idr for more information."

# print help function.
printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$rep1peakslistHelp"
    echo -e ""
    echo -e "$rep2peakslistHelp"
    echo -e ""
    echo -e "$peakscallerHelp"
    echo -e ""
    echo -e "$filetypeHelp"
    echo -e ""
    echo -e "$rankHelp"
    echo -e ""
    echo -e "$outputHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-Analyzer idr' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'h1:2:c:t:r:o:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        1) peaksrep1=$OPTARG ;;
        2) peaksrep2=$OPTARG ;;
	c) peakscaller=$OPTARG ;;
        t) filetype=$OPTARG ;;
        r) rank=$OPTARG ;;
        o) outputfile=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# required options.
if test -z $peaksrep1 ;then
   echo ""
   echo "Please provide the KAS-seq peaks file of the first replicate. -1 [rep1 peaklist]"
   echo ""
   exit 1
fi

if test -z $peaksrep2 ;then
   echo ""
   echo "Please provide the KAS-seq peaks file of the second replicate. -2 [rep2 peaklist]"
   echo ""
   exit 1
fi

if test -z $peakscaller ;then
   echo ""
   echo "Please specify the peaks caller used to generate the peak list. macs2 or epic2. -c [peaks caller]"
   echo ""
   exit 1
fi

if test -z $outputfile ;then
   echo ""
   echo "Please specify the file where the output should be written. -o [output file]"
   echo ""
   exit 1
fi

# setup the default options.

if test -z $filetype ;then
   filetype="broadPeak"
fi

if test -z $rank ;then
   rank="signal.value"
fi

# test if idr (Irreproducible Discovery Rate) framework was installed.
if ! type idr > /dev/null 2>&1 ;then
   echo "idr was not installed or not export to the \$PATH'"
   echo ""
   echo "Install idr framework with 'conda install -c bioconda idr' or 'mamba install idr'."
   echo ""
   exit 1
fi

# test unsuported peaks caller. macs2 and epic2

if [[ $peakscaller != "macs2" ]] && [[ $peakscaller != "epic2" ]] ;then
    echo ""
    echo "Error: unsupported peaks caller: $peakscaller. e.g. epic2 or macs2."
    echo ""
    exit -1
fi

# test unsuported output file type. narrowPeak,broadPeak,bed
if test -n "$filetype" && [[ $filetype != "narrowPeak" ]] && [[ $filetype != "broadPeak" ]] && [[ $filetype != "bed" ]] ;then
    echo ""
    echo "Error: unsupported output file type: $filetype. e.g. narrowPeak, broadPeak or bed. DEFAULT: broadPeak."
    echo ""
    exit -1
fi

# test unsuported signal types to rank peaks. signal.value, p.value or q.value.
if test -n "$rank" && [[ $rank != "signal.value" ]] && [[ $rank != "p.value" ]] ;then
    echo ""
    echo "Error: unsupported signal types to rank peaks: $rank. e.g. signal.value or p.value. DEFAULT: signal.value."
    echo ""
    exit -1
fi

# get the path of 'KAS-Analyzer idr' shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# # sort the peak list using $rank, signal.value or p.value.

if [[ $peakscaller == "epic2" ]] ;then
	if [[ $rank == "signal.value" ]] ;then
		
		sed "1d" $peaksrep1 | awk '{printf("%s\t%d\t%d\t%s\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%d\n",$1,$2,$3,"epic2_peaks"NR,log($4)/log(10)*(-10),".",$10,log($4)/log(10)*(-1),log($9)/log(10)*(-1),$7)}' | sed "s/+inf/1000/g" | sort -k7,7nr > ${peaksrep1}_macs2_10_columns.sorted.bed
		sed "1d" $peaksrep2 | awk '{printf("%s\t%d\t%d\t%s\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%d\n",$1,$2,$3,"epic2_peaks"NR,log($4)/log(10)*(-10),".",$10,log($4)/log(10)*(-1),log($9)/log(10)*(-1),$7)}' | sed "s/+inf/1000/g" | sort -k7,7nr > ${peaksrep2}_macs2_10_columns.sorted.bed

	elif [[ $rank == "p.value" ]] ;then

		sed "1d" $peaksrep1 | awk '{printf("%s\t%d\t%d\t%s\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%d\n",$1,$2,$3,"epic2_peaks"NR,log($4)/log(10)*(-10),".",$10,log($4)/log(10)*(-1),log($9)/log(10)*(-1),$7)}' | sed "s/+inf/1000/g" | sort -k8,8nr > ${peaksrep1}_macs2_10_columns.sorted.bed
		sed "1d" $peaksrep2 | awk '{printf("%s\t%d\t%d\t%s\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%d\n",$1,$2,$3,"epic2_peaks"NR,log($4)/log(10)*(-10),".",$10,log($4)/log(10)*(-1),log($9)/log(10)*(-1),$7)}' | sed "s/+inf/1000/g" | sort -k8,8nr > ${peaksrep2}_macs2_10_columns.sorted.bed
	
	fi	

elif [[ $peakscaller == "macs2" ]] ;then

	if [[ $rank == "signal.value" ]] ;then
		echo "hahahahaha"
		
		sort -k7,7nr $peaksrep1 > ${peaksrep1}_macs2_10_columns.sorted.bed
		sort -k7,7nr $peaksrep2 > ${peaksrep2}_macs2_10_columns.sorted.bed

	elif [[ $rank == "p.value" ]] ;then
         
		sort -k8,8nr $peaksrep1 > ${peaksrep1}_macs2_10_columns.sorted.bed
		sort -k7,7nr $peaksrep2 > ${peaksrep2}_macs2_10_columns.sorted.bed
	fi

fi
  
# Utilize the Irreproducible Discovery Rate (IDR) framework to assess the reproducibility of KAS-seq peaks identified from replicate KAS-seq experiments.

idr --samples ${peaksrep1}_macs2_10_columns.sorted.bed ${peaksrep2}_macs2_10_columns.sorted.bed --input-file-type narrowPeak --output-file-type $filetype --rank $rank --output-file ${outputfile}_idr --plot --log-output-file ${outputfile}_idr.log

rm -rf ${peaksrep1}_macs2_10_columns.sorted.bed
rm -rf ${peaksrep2}_macs2_10_columns.sorted.bed

echo "'KAS-Analyzer idr' run successfully."

