#!/bin/bash
# 'KAS-Analyzer UCSC' was developed by Ruitu Lyu on 12-10-2021.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-Analyzer UCSC [ -h/--help ] [ -k KAS-seq ] [ -l UCSC track ] [ -c track colors ]"
exampleHelp="Example: nohup KAS-Analyzer UCSC -k KAS-seq_data.txt -l UCSC_track_names.txt &"
KASseqHelp="-k [KAS-seq]: please input the text file containing the bedGraph files generated from 'KAS-Analyzer KAS-seq'. REQUIRED.
Example:
KAS-seq_WT.rep1.nor.bg
KAS-seq_WT.rep2.nor.bg
KAS-seq_KO.rep1.nor.bg 
KAS-seq_KO.rep2.nor.bg    ---KAS-seq_data.txt"
tracknameHelp="-l [UCSC track]: please input the text file containing the track names of KAS-seq or spKAS-seq data that you want to visualize on UCSC genome browser. REQUIRED.
Example: 
KAS-seq_WT.rep1 
KAS-seq_WT.rep2
KAS-seq_KO.rep1
KAS-seq_KO.rep2           ---UCSC_track_names.txt"
trackcolorHelp=" -c [track colors]: please input the text file containing the R,G,B colors list for KAS-seq or spKAS-seq data that you want to visualize on UCSC genome browser. Default: black (0,0,0).
Example:
255,102,102
255,178,102
102,255,102
102,178,255               ---track_colors.txt

For more colors options, please refer to https://www.w3schools.com/colors/colors_rgb.asp or https://www.rapidtables.com/web/color/RGB_Color.html.
"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-Analyzer UCSC' shell script is used to generate files for uploading into UCSC genome browser."

# print help function.
printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$tracknameHelp"
    echo -e ""
    echo -e "$trackcolorHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters or option '--help' was provided, 'KAS-Analyzer UCSC' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then	
   printHelpAndExit
fi

# get the value of options.
while getopts 'hk:l:c:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        k) KASseq=$OPTARG ;;
        l) trackname=$OPTARG ;;
	c) trackcolor=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# Required options.
if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing the normalized bedGraph files that you want to upload into UCSC genome browser. REQUIRED. -k [KAS-seq_data.txt]"
   echo ""
   printHelpAndExit 0
fi

if test -z $trackname ;then
   echo ""
   echo "Please input the text file containing the track names of KAS-seq or spKAS-seq data that you want to visualize on UCSC genome browser. REQUIRED. -n [UCSC_track_names.txt]"
   echo ""
   printHelpAndExit 0
fi

# get the number of KAS-seq samples.
number_of_samples=$( awk 'END {print NR}' $KASseq )

# get the number of track names.
number_of_tracknames=$( awk 'END {print NR}' $trackname )

# Test if the number of track names is consistent with the number of KAS-seq samples.
if [[ ${number_of_tracknames} != ${number_of_samples} ]]
   then
   echo ""
   echo "Error: the number of track names isn't consistent with the number of KAS-seq or spKAS-seq samples."
   echo ""
   exit
fi

# setup the default option and parameters.
if test -z $trackcolor ;then
   for ((i=1; i<=${number_of_samples}; i++))
   do
   echo "0,0,0" >> .track_colors.txt
   done
   trackcolor=".track_colors.txt"
else
   number_of_colors=$( awk 'END {print NR}' $trackcolor )
       if [[ ${number_of_tracknames} != ${number_of_colors} ]] ;then
       echo ""
       echo "Error: the number of track colors isn't consistent with the number of KAS-seq or spKAS-seq samples and track names."
       echo ""
       exit
       fi
fi

# Test if the number of track colors is consistent with the number of KAS-seq samples or tracks.
# if [[ ${number_of_tracknames} != ${number_of_colors} ]]
# then
#   echo " "
#   echo "Error: the number of track colors isn't consistent with the number of KAS-seq or spKAS-seq samples and track names."
#   exit
# fi

# generate files for uploading into UCSC genome browser.
for ((j=1; j<=${number_of_samples}; j++))
do
sample_selected=$( sed -n ''$j'p' $KASseq )
trackname_selected=$( sed -n ''$j'p' $trackname)
KASseq_basename=$(basename ${sample_selected} .bg)
color_selected=$( sed -n ''$j'p' $trackcolor)
echo "Processing $sample_selected ..."
echo ""
echo "track type=bedGraph name=\"${trackname_selected}.$(date +%Y-%m-%d)\" description=\"${trackname_selected}.$(date +%Y-%m-%d)\" visibility=full color=${color_selected}" > .${trackname_selected}.track
cat .${trackname_selected}.track $sample_selected | head -n 50000000 > ${KASseq_basename}.UCSC.bg
gzip ${KASseq_basename}.UCSC.bg
rm -f .${trackname_selected}.track
echo "done."
echo ""
done

rm -f .track_colors.txt

echo "'KAS-Analyzer UCSC' run successfully!"
