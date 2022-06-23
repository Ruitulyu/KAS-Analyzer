#!/bin/bash
# 'KAS-pipe2 heatmap' was developed by Ruitu Lyu on 12-14-2021.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 heatmap [ -h/--help ] [ -t threads ] [ -e length ] [ -s assembly id ] [ -q ] [ -u samples using ] [ -m maximum value ] [ -o prefix ] [ -r regions ] [ -p peaks ] [ -l labels ] [ -c colors ] [ -k KAS-seq ]"
exampleHelp="Example: 
Genomic features(genebody, TSS or TES):
nohup KAS-pipe2 heatmap -t 10 -s hg19 -o KAS-seq_heatmap -r genebody -q -c Reds -l labels.txt -k KAS-seq.txt &
Custom regions(peaks. e.g. enhancers.bed):
nohup KAS-pipe2 heatmap -t 10 -o KAS-seq_heatmap -r peaks -q -p KAS-seq_peaks.bed -c Reds,Reds,Blues,Blues -l labels.txt -k KAS-seq.txt &"
threadsHelp="-t [threads]: please specify the number of threads used for generating (sp)KAS-seq heatmap plot. DEFAULT: 1."
assemblyidHelp="-s [assemblyid]: please specify the genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED only for 'genomic features' mode."
sortHelp="-q: please specify to sort the regions. DEFAULT: off."
sortsampleusingHelp="-u [sample using]: please specify the samples list to sort the regions. e.g. -u 1,2,3"
lengthHelp="-e [length]: please specify the distance upstream of the start site of the regions defined in the region file. If the regions are genebody, this would be the distance upstream of the transcription start site. DEFAULT: 3000."
maximumvaluesHelp="-m [maximum value]: please specify the maximum value of the heatmap intensities. e.g. -m 15,20,60. Note: the number of values need to be consistent with KAS-seq samples."
prefixHelp="-o [KAS-seq_heatmap]: please specify the prefix (basename) of 'KAS-pipe2 heatmap' output files. REQUIRED."
regionsHelp="-r [regions]: please specify the regions types to generate heatmap plot. e.g. Genomic features: genebody, TSS or TES; Custom regions: peaks. REQUIRED."
peaksHelp="-p [peakslist]: please input the peak list file. REQUIRED only for 'custom regions' mode."
colorsHelp="-c [colors]: please specify the colors for (sp)KAS-seq data in heatmap plot. Note: you can specify only one colors (e.g. Reds) or specify a color list (e.g. Reds,Greens,Blues,Purples), which needs to be consistent with the number of KAS-seq bigWig files. REQUIRED. Note: please refer to http://matplotlib.org/users/colormaps.html to get the list of valid colors names."
labelsHelp="-l [labels.txt]: please input the text file containing the labels of (sp)KAS-seq data to generate metagene profile. Default: basename of (sp)KAS-seq bigWig files.
Example:
WT_rep1
WT.rep2
KO.rep1
KO.rep2                        ---labels.txt"
KASseqHelp="-k [KAS-seq.txt]: please input the text file containing normalized (sp)KAS-seq bigWig files, which can be generated with 'KAS-pipe2 normalize' and 'KAS-pipe2 bedGraphToBigWig' shell scripts. The order and number of (sp)KAS-seq bigWig files should be the consistent with the labels file when specified. REQUIRED.
Example:
KAS-seq_WT_rep1.nor.bigWig
KAS-seq_WT_rep2.nor.bigWig
KAS-seq_KO_rep1.nor.bigWig
KAS-seq_KO_rep2.nor.bigWig     ---KAS-seq.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 heatmap' shell script is applied to generate heatmap plots for (sp)KAS-seq data on genomic features( genebody, TSS or TES) or provided custom regions. 'KAS-pipe2 heatmap' shell script mainly invoke deeptools 'computeMatrix' and 'plotHeatmap', please refer to https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html for more information."

# print help function.
printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$threadsHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$sortHelp"
    echo -e ""
    echo -e "$sortsampleusingHelp"
    echo -e ""
    echo -e "$lengthHelp"
    echo -e ""
    echo -e "$maximumvaluesHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$regionsHelp"
    echo -e ""
    echo -e "$peaksHelp"
    echo -e ""
    echo -e "$colorsHelp"
    echo -e ""
    echo -e "$labelsHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 heatmap' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'ht:s:qu:e:m:o:r:p:l:c:k:' opt; do
    case $opt in
        h) printHelpAndExit;;
        t) threads=$OPTARG ;;
        s) assemblyid=$OPTARG ;;
	q) sort="on" ;;
	u) sortsampleusing=$OPTARG ;;
	e) length=$OPTARG ;;
	m) maximumvalues=$OPTARG ;;
        o) prefix=$OPTARG ;;
        r) regions=$OPTARG ;;
        p) peakslist=$OPTARG ;;
        l) labels=$OPTARG ;;
        c) colors=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        ?) printHelpAndExit;;
    esac
done

# test if deeptools was installed.
if ! type deeptools > /dev/null 2>&1 ;then
   echo "deeptools was not installed or not export to the \$PATH'"
   echo ""
   echo "Install deeptools with 'conda install -c bioconda deeptools' or refer the official website of 'deeptools'."
   echo ""
   exit 1
fi

# Required options.
if test -z $prefix ;then
   echo ""
   echo "Please input the prefix (basename) of 'KAS-pipe2 heatmap' output files. -o [prefix]"
   echo ""
   exit -1
fi

if test -z $regions ;then
   echo ""
   echo "Please specify the region types for generating heatmap. e.g. Genomic features: genebody, TSS or TES; Custom regions: peaks."
   echo ""
   exit -1

elif [[ $regions == "genebody" ]] || [[ $regions == "TSS" ]] || [[ $regions == "TES" ]] ;then
   
   if test -z $assemblyid ;then 
   echo ""
   echo "please specify the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebr
a fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. -s [assembly id]; REQUIRED only for 'genomic features' mode."
   echo ""
   exit -1
   
   elif [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]]; then
   echo ""
   echo "Error: unsupported assembly id: $assemblyid." 
   echo ""
   exit -1
   fi 

elif [[ $regions == "peaks" ]] ;then 
   if test -z $peakslist ;then
   echo ""
   echo "please input the peak list file. REQUIRED only for 'custom regions' mode. -p [peakslist]"
   echo ""
   exit -1

   elif test -n "$peakslist" ;then
   peakslist_basename=$( basename ${peakslist} .bed )	   
   echo ""
   fi
else 
echo ""
echo "Error: unsupported regions types: $regions. Please specify the region types for generating heatmap. e.g. Genomic features: genebody, TSS or TES; Custom regions: peaks." 
echo ""
exit -1

fi

if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing normalized (sp)KAS-seq bigWig files. -k [KAS-seq]"
   echo ""
   exit -1
fi

number_of_samples=$( awk 'END {print NR}' $KASseq )

if test -n "$labels" ;then
number_of_labels=$( awk 'END {print NR}' $labels )
   if [[ ${number_of_samples} != ${number_of_labels} ]] ;then
      echo ""
      echo "Error:the number of labels isn't consistent with the number of samples!"
      echo ""
      exit -1
   fi
fi 

if test -z $colors ;then 
   echo ""
   echo "please specify the color list for (sp)KAS-seq data in heatmap plot. REQUIRED. -c [colors]"
   echo ""
   exit -1
else
   echo $colors > .colors.txt
   sed -i "s/\,/ /g" .colors.txt
   colors_list=$( cat .colors.txt )
   number_of_colors=$( awk 'END {print NF}' .colors.txt )
   if [[ $number_of_colors -ge 2 ]] && [[ $number_of_colors != $number_of_samples ]] ;then
       echo ""
       echo "Error:the number of colors isn't consistent with the number of samples!"
       echo ""
       exit
   fi    
rm -f .colors.txt
fi

# setup the maximum and minimum heatmap density value.
if test -n "$maximumvalues" ;then
   echo $maximumvalues > .maximumvalues.txt 	
   sed -i "s/\,/ /g" .maximumvalues.txt
   maximumvalues_list=$( cat .maximumvalues.txt )
   number_of_maximumvalues=$( awk 'END {print NF}' .maximumvalues.txt )
     if [[ $number_of_maximumvalues -ge 2 ]] && [[ $number_of_maximumvalues != $number_of_samples ]] ;then
	 echo ""
	 echo "Error:the number of maximum values isn't consistent with the number of samples!"   
	 echo ""
	 exit
     fi     
   minimumvalues_list=$( echo "0" | awk '{for(j=1;j<='$number_of_samples';j++) $j=$1;print}' )

   rm -f .maximumvalues.txt
fi 


if test -n "$sortsampleusing" ;then
   echo $sortsampleusing > .sortsampleusing.txt
   sed -i "s/\,/ /g" .sortsampleusing.txt
   sortsampleusing_list=$( cat .sortsampleusing.txt )
   rm -f .sortsampleusing.txt
fi

# setup the default parameters.

# setup the basename of bigWig files as default labels.
if test -z $labels ;then
   for ((i=1; i<=${number_of_samples}; i++))
   do
   sample_selected=$(sed -n ''$i'p' $KASseq)
   label_basename=$(basename ${sample_selected} .bigWig)
   echo $label_basename >> .labels_basename.txt
   done
   labels=".labels_basename.txt"
fi

# setup the default number of threads as 1.
if test -z $threads ;then
   threads=1
fi

if test -z $length ;then
   length=3000
fi

# default status of sort is off. 
if test -z $sort ;then
   sort="off"
fi

# get the sample and labels list. 
samples_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $KASseq)
labels_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $labels)
rm -f .labels_basename.txt

# get the absolute path of 'KAS-pipe2 heatmap' shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

if [[ $regions == "genebody" ]] && [[ $sort == "on" ]] ;then
   # generate the matrix file on genebody.
   echo "Compute (sp)KAS-seq matrix of $samples_list on ${assemblyid}_Refseq_${regions} ..."
   echo ""
   computeMatrix scale-regions -R ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.bed -S $samples_list -b $length -a $length --regionBodyLength 6000 --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz > /dev/null 2>&1
   echo "done."
   echo ""

   if test -z $sortsampleusing && test -z $maximumvalues ;then
      echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with all (sp)KAS-seq data ..."	
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
      echo "done."
      echo ""

      echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with all (sp)KAS-seq data ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
      echo "done."
      echo ""
        
   elif test -z $sortsampleusing && test -n "$maximumvalues" ;then
      echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with all (sp)KAS-seq data ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
      echo "done."
      echo ""

      echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with all (sp)KAS-seq data ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
      echo "done."
      echo ""
        
   elif test -n "$sortsampleusing" && test -z $maximumvalues ;then
      echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
      echo "done."
      echo ""

      echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
      echo "done."
      echo ""

   elif test -n "$sortsampleusing" && test -n "$maximumvalues" ;then
      echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
      echo "done."
      echo ""

      echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
      echo "done."
      echo ""

   fi

# rm -f {prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz
elif [[ $regions == "genebody" ]] && [[ $sort == "off" ]] ;then
# generate the matrix file on genebody.
   echo "Compute (sp)KAS-seq matrix of $samples_list on ${assemblyid}_Refseq_${regions} ..."
   echo ""
   computeMatrix scale-regions -R ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.bed -S $samples_list -b $length -a $length --regionBodyLength 6000 --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz > /dev/null 2>&1
   echo "done."
   echo ""
   
   if test -z $maximumvalues ;then
      echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' in original ${assemblyid} ${regions} order ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
      echo "done."
      echo ""

      echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' in original ${assemblyid} ${regions} order ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
      echo "done."
      echo ""

   elif test -n "$maximumvalues" ;then
      echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' in original ${assemblyid} ${regions} order ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
      echo "done."
      echo ""

      echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' in original ${assemblyid} ${regions} order ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
      echo "done."
      echo ""
        
   fi

# rm -f ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz
elif [[ $regions == "TSS" ]] && [[ $sort == "on" ]] ;then
   # generate the matrix file on TSS.	
   echo "Compute (sp)KAS-seq matrix of $samples_list on ${assemblyid}_Refseq_${regions} ..."
   echo ""
   computeMatrix reference-point --referencePoint TSS -R ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.bed -S $samples_list -b $length -a $length --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz > /dev/null 2>&1
   echo "done."
   echo ""
		
   if test -z $sortsampleusing && test -z $maximumvalues ;then
      echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with all (sp)KAS-seq data ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
      echo "done."
      echo ""
        
      echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with all (sp)KAS-seq data ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
      echo "done."
      echo ""

   elif test -z $sortsampleusing && test -n "$maximumvalues" ;then
      echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with all (sp)KAS-seq data ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
      echo "done."
      echo ""

      echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with all (sp)KAS-seq data ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
      echo "done."
      echo ""

   elif test -n "$sortsampleusing" && test -z $maximumvalues ;then
      echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
      echo "done."
      echo ""

      echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
      echo "done."
      echo ""
        
      elif test -n "$sortsampleusing" && test -n "$maximumvalues" ;then
      echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
      echo ""
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
      echo "done."
      echo ""

      echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
      plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
      echo "done."
      echo ""

      fi

# rm -f ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz	
elif [[ $regions == "TSS" ]] && [[ $sort == "off" ]] ;then	
   # generate the matrix file on TSS.
   echo "Compute (sp)KAS-seq matrix of $samples_list on ${assemblyid}_Refseq_${regions} ..."
   echo ""
   computeMatrix reference-point --referencePoint TSS -R ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.bed -S $samples_list -b $length -a $length --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz > /dev/null 2>&1
   echo "done."
   echo ""
        
   if test -z $maximumvalues ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' in original ${assemblyid} ${regions} order ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""

       echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' in original ${assemblyid} ${regions} order ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""

    elif test -n "$maximumvalues" ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' in original ${assemblyid} ${regions} order ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""

       echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' in original ${assemblyid} ${regions} order ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""

   fi
# rm -f ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz

elif [[ $regions == "TES" ]] && [[ $sort == "on" ]] ;then
   # generate the matrix file on TES.        
   echo "Compute (sp)KAS-seq matrix of $samples_list on ${assemblyid}_Refseq_${regions} ..."	
   echo ""
   computeMatrix reference-point --referencePoint TES -R ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.bed -S $samples_list -b $length -a $length --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz > /dev/null 2>&1
   echo "done."
   echo ""	
	
   if test -z $sortsampleusing && test -z $maximumvalues ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with all (sp)KAS-seq data ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""
        
       echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with all (sp)KAS-seq data ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""

   elif test -z $sortsampleusing && test -n "$maximumvalues" ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with all (sp)KAS-seq data ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""	

       echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with all (sp)KAS-seq data ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo "" 
        
   elif test -n "$sortsampleusing" && test -z $maximumvalues ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples 
$sortsampleusing_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""

       echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples 
$sortsampleusing_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""

       elif test -n "$sortsampleusing" && test -n "$maximumvalues" ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""

       echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""
        
   fi
# rm -f ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz

elif [[ $regions == "TES" ]] && [[ $sort == "off" ]] ;then     
   # generate the matrix file on TES.
   echo "Compute (sp)KAS-seq matrix of $samples_list on ${assemblyid}_Refseq_${regions} ..."
   echo ""
   computeMatrix reference-point --referencePoint TES -R ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.bed -S $samples_list -b $length -a $length --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz > /dev/null 2>&1
   echo "done."
   echo ""
   if test -z $maximumvalues ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' in original ${assemblyid} ${regions} order ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""

       echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' in original ${assemblyid} ${regions} order ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""

   elif test -n "$maximumvalues" ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' in original ${assemblyid} ${regions} order ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""

       echo "Generate the svg format heatmap plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' in original ${assemblyid} ${regions} order ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""
        
   fi

# rm -f ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz	
elif [[ $regions == "peaks" ]] && [[ $sort == "on" ]] ;then
   # generate the matrix file on user provided peak list.	
   echo "Compute (sp)KAS-seq matrix of $samples_list on ${peakslist} ..."
   echo ""
   computeMatrix reference-point --referencePoint center -R $peakslist -S $samples_list -b $length -a $length --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${peakslist_basename}.matrix.gz > /dev/null 2>&1
   echo "done."
   echo ""

   if test -z $sortsampleusing && test -z $maximumvalues ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${peakslist_basename}.matrix.gz', and sort with all (sp)KAS-seq data ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${peakslist_basename}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" -out ${prefix}_on_${peakslist_basename}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""

       echo "Generate the svg format heatmap plot with '${prefix}_on_${peakslist_basename}.matrix.gz', and sort with all (sp)KAS-seq data ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${peakslist_basename}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" -out ${prefix}_on_${peakslist_basename}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""

   elif test -z $sortsampleusing && test -n "$maximumvalues" ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${peakslist_basename}.matrix.gz', and sort with all (sp)KAS-seq data ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${peakslist_basename}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${peakslist_basename}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""

       echo "Generate the svg format heatmap plot with '${prefix}_on_${peakslist_basename}.matrix.gz', and sort with all (sp)KAS-seq data ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${peakslist_basename}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${peakslist_basename}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""

   elif test -n "$sortsampleusing" && test -z $maximumvalues ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${peakslist_basename}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${peakslist_basename}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list -out ${prefix}_on_${peakslist_basename}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""

       echo "Generate the svg format heatmap plot with '${prefix}_on_${peakslist_basename}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${peakslist_basename}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list -out ${prefix}_on_${peakslist_basename}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""

   elif test -n "$sortsampleusing" && test -n "$maximumvalues" ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${peakslist_basename}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${peakslist_basename}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${peakslist_basename}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""

       echo "Generate the svg format heatmap plot with '${prefix}_on_${peakslist_basename}.matrix.gz', and sort with (sp)KAS-seq data number: $sortsampleusing_list ..."
       plotHeatmap -m ${prefix}_on_${peakslist_basename}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortUsingSamples $sortsampleusing_list --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${peakslist_basename}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""

   fi
# rm -f ${prefix}_on_${peakslist_basename}.matrix.gz	
	
elif [[ $regions == "peaks" ]] && [[ $sort == "off" ]] ;then	
   # generate the matrix file on user provided peak list.
   echo "Compute (sp)KAS-seq matrix of $samples_list on ${peakslist} ..."
   echo ""
   computeMatrix reference-point --referencePoint center -R $peakslist -S $samples_list -b $length -a $length --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${peakslist_basename}.matrix.gz > /dev/null 2>&1
   echo "done."
   echo ""
        
   if test -z $maximumvalues ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${peakslist_basename}.matrix.gz' in original ${peakslist} order ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${peakslist_basename}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no -out ${prefix}_on_${peakslist_basename}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""

       echo "Generate the svg format heatmap plot with '${prefix}_on_${peakslist_basename}.matrix.gz' in original ${peakslist} order ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${peakslist_basename}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no -out ${prefix}_on_${peakslist_basename}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""

   elif test -n "$maximumvalues" ;then
       echo "Generate the png format heatmap plot with '${prefix}_on_${peakslist_basename}.matrix.gz' in original ${peakslist} order ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${peakslist_basename}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${peakslist_basename}.matrix_heatmap.png --plotFileFormat png
       echo "done."
       echo ""

       echo "Generate the svg format heatmap plot with '${prefix}_on_${peakslist_basename}.matrix.gz' in original ${peakslist} order ..."
       echo ""
       plotHeatmap -m ${prefix}_on_${peakslist_basename}.matrix.gz --colorMap $colors_list --boxAroundHeatmaps no --whatToShow "heatmap and colorbar" --sortRegions no --zMin $minimumvalues_list --zMax $maximumvalues_list -out ${prefix}_on_${peakslist_basename}.matrix_heatmap.svg --plotFileFormat svg
       echo "done."
       echo ""

   fi

# rm -f ${prefix}_on_${peakslist_basename}.matrix.gz
fi

echo "'KAS-pipe2 heatmap' run successfully!"
