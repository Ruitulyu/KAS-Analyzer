#!/bin/bash
# 'KAS-pipe2 profile' was developed by Ruitu Lyu on 12-14-2021.

# Stop on error
set -e

# help arguments
usageHelp="Usage: KAS-pipe2 profile [ -h/--help ] [ -t threads ] [ -s assembly id ] [ -e length ] [ -o prefix ] [ -r regions ] [ -p peaks file ] [ -f peaks files list ] [ -l labels ] [ -c colors ] [ -k KAS-seq ]"
exampleHelp="Example: 
Genomic features(genebody, TSS or TES):
nohup KAS-pipe2 profile -t 10 -s hg19 -o KAS-seq_genebody -r genebody -c red,blue,green,purple -l labels.txt -k KAS-seq.txt &

Custom regions(peaks. e.g. enhancers.bed):
nohup KAS-pipe2 profile -t 10 -o KAS-seq_peaks -r peaks -p KAS-seq_peaks.bed -c red,blue,green,purple -l labels.txt -k KAS-seq.txt &

KAS-seq signal on different regions (-f [peaks list] must be specified):
nohup KAS-pipe2 profile -t 10 -o KAS-seq_different_clusters -r peakslist -f peaks_cluster1.bed,peaks_cluster2.bed,peaks_cluster3.bed -c red,green,purple -l labels.txt -k KAS-seq.txt &"
threadsHelp="-t [threads]: please input the number of threads used for generating (sp)KAS-seq metagene profile plot. DEFAULT: 1."
assemblyidHelp="-s [assemblyid]: please specify the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. REQUIRED only for 'genomic features' mode."
lengthHelp="-e [length]: please specify the distance upstream of the start site of the regions defined in the region file. If the regions are genebody, this would be the distance upstream of the transcription start site. DEFAULT: 3000."
prefixHelp="-o [KAS-seq_profile]: please input the prefix (basename), which will be used to generate the name of 'KAS-pipe2 metageneprofile' output files. REQUIRED."
regionsHelp="-r [regions]: please specify the regions types for generating metagene profile plot. e.g. Genomic features: genebody, TSS or TES; Custom regions: peaks or peakslist. REQUIRED."
peaksHelp="-p [peaks file]: please specify the peak file. REQUIRED only for 'custom regions: peaks' mode."
peakslistHelp="-f [peaks list]: please specify the peak files list. e.g. peaks_cluster1.bed,peaks_cluster2.bed,peaks_cluster3.bed. REQUIRED only for 'custom regions: peakslist' mode. Note: only one KAS-seq bigWig file is needed."
colorsHelp="-c [colors]: please specify the color list for (sp)KAS-seq data in metagene profile plot. Note: the number of colors in the profile plot needs to be consistent with the number of KAS-seq bigWig files. REQUIRED. Note: the list of valid color names https://matplotlib.org/examples/color/named_colors.html."
labelsHelp="-l [labels.txt]: please input the text file containing the labels of (sp)KAS-seq data or peaks files (-f need to be specified) that used for generating metagene profile. DEFAULT: basename of (sp)KAS-seq bigWig files or peaks files.
Example:
WT_rep1          Cluster1
WT.rep2          Cluster2
KO.rep1          Cluster3
KO.rep2   or                        ---labels.txt"
KASseqHelp="-k [KAS-seq.txt]: please input the text file containing normalized (sp)KAS-seq bigWig files, which can be generated with 'KAS-pipe2 normalize' and 'KAS-pipe2 bedGraphToBigWig' shell scripts. The order and number of (sp)KAS-seq bigWig files should be the consistent with the labels file when provided. REQUIRED.
Example:
KAS-seq_WT_rep1.nor.bigWig
KAS-seq_WT_rep2.nor.bigWig
KAS-seq_KO_rep1.nor.bigWig
KAS-seq_KO_rep2.nor.bigWig          ---KAS-seq.txt"
helpHelp="-h/--help: print this help and exit.
Note: The 'KAS-pipe2 profile' shell script is applied to generate metagene profile for (sp)KAS-seq data on genomic features (genebody, TSS or TES) or provided custom regions. 'KAS-pipe2 metageneprofile' shell script mainly invoke deeptools 'computeMatrix' and 'plotProfile', please refer to https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html for more information."

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
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$lengthHelp"
    echo -e ""
    echo -e "$regionsHelp"
    echo -e ""
    echo -e "$peaksHelp"
    echo -e ""
    echo -e "$peakslistHelp"
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

# if no parameters was provided, 'KAS-pipe2 metageneprofile' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'ht:s:o:e:r:p:f:l:c:k:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        t) threads=$OPTARG ;;
	s) assemblyid=$OPTARG ;;
        o) prefix=$OPTARG ;;
	e) length=$OPTARG ;;
        r) regions=$OPTARG ;;
        p) peaks=$OPTARG ;;
	f) peakslist=$OPTARG ;;
        l) labels=$OPTARG ;;
        c) colors=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        ?) printHelpAndExit 0;;
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
   echo "Please input the prefix (basename) of 'KAS-pipe2 metageneprofile' output files. -o [prefix]"
   echo ""
   exit -1
fi

if test -z $regions ;then
   echo ""
   echo "Please specify the region types for generating metagene profile plot. e.g. Genomic features: genebody, TSS or TES; Custom regions: peaks or peakslist. -r [regions]."
   echo ""
   exit -1

elif [[ $regions == "genebody" ]] || [[ $regions == "TSS" ]] || [[ $regions == "TES" ]] ;then
   
   if test -z $assemblyid ;then 
      echo ""
      echo "please specify the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. -s [assembly id]; REQUIRED only for 'genomic features' mode."
      echo ""
      exit -1

   elif [[ $assemblyid != "hg18" ]] && [[ $assemblyid != "hg19" ]] && [[ $assemblyid != "hg38" ]] && [[ $assemblyid != "mm9" ]] && [[ $assemblyid != "mm10" ]] && [[ $assemblyid != "mm39" ]] && [[ $assemblyid != "dm3" ]] && [[ $assemblyid != "dm6" ]] && [[ $assemblyid != "rn6" ]] && [[ $assemblyid != "rn7" ]] && [[ $assemblyid != "ce10" ]] && [[ $assemblyid != "ce11" ]] && [[ $assemblyid != "danRer10" ]] && [[ $assemblyid != "danRer11" ]]; then
      echo ""
      echo "Error: unsupported assembly id: $assemblyid. supported assembly id: Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11." 
      echo ""
      exit -1
   fi 

elif [[ $regions == "peaks" ]] ;then

   if test -z $peaks && test -z $peakslist ;then
      echo ""
      echo "Please specify the peaks file or peaks files list. -p [peaks] or -f [peaks file list]"
      echo ""
      exit -1
   fi

else 
   echo ""
   echo "Error: unsupported regions types: $regions. please specify the regions types for generating metagene profile plot. e.g. Genomic features: genebody, TSS or TES; Custom regions: peaks. REQUIRED. -r [regions]." 
   echo ""
   exit -1
fi 

if test -z $KASseq ;then
   echo ""
   echo "Please input the text file containing normalized (sp)KAS-seq bigWig files. -k [KAS-seq]"
   echo ""
   exit -1
fi

# get the number of samples.
number_of_samples=$( awk 'END {print NR}' $KASseq )

if test -z $colors ;then 
   echo ""
   echo "please specify the color list for (sp)KAS-seq data in metagene profile plot. REQUIRED. -c [colors]"
   echo ""
   exit -1
fi

echo $colors > ${prefix}.colors.txt
sed -i "s/\,/ /g" ${prefix}.colors.txt
colors_list=$( cat ${prefix}.colors.txt )
number_of_colors=$( awk 'END {print NF}' ${prefix}.colors.txt )
rm -f ${prefix}.colors.txt

if test -n "$peakslist" ;then
   echo $peakslist > ${prefix}.peakslist.txt
   sed -i "s/\,/ /g" ${prefix}.peakslist.txt
   peaks_list=$( cat ${prefix}.peakslist.txt )
   number_of_peaksfiles=$( awk 'END {print NF}' ${prefix}.peakslist.txt )
fi

# test if number of colors is consistent with the number of samples or peaks files.
if test -z $peakslist && [[ ${number_of_samples} != ${number_of_colors} ]] ;then
   echo ""
   echo "Error: the number of colors isn't consistent with the number of samples!"
   echo ""
   exit -1

elif test -n "$peakslist" && [[ ${number_of_peaksfiles} != ${number_of_colors} ]] ;then
   echo ""
   echo "Error: the number of colors isn't consistent with the number of peaks files!"
   echo ""
   exit -1
fi

# setup the default parameters.

# setup the basename of bigWig files or peak files as default labels.
if test -z $labels ;then
   if test -z $peakslist ;then
      for ((i=1; i<=${number_of_samples}; i++))
      do
      sample_selected=$( sed -n ''$i'p' $KASseq )
      label_basename=$( basename ${sample_selected} .bigWig )
      echo $label_basename >> ${prefix}.labels_basename.txt
      done
      labels="${sample_selected}.labels_basename.txt"
     
   elif test -n "$peakslist" ;then
      for ((j=1; j<=${number_of_peaksfiles}; j++))
      do
      peaks_selected=$( awk -v a=$j '{print $a}' ${prefix}.peakslist.txt )	     
      label_basename=$( basename ${peaks_selected} .bed )
      echo $label_basename >> ${prefix}.labels_basename.txt
      done
      labels="${prefix}.labels_basename.txt"
      rm -f ${prefix}.peakslist.txt
   fi 
else 
   echo ""
fi


# setup the default number of threads as 1.
if test -z $threads ;then
   threads=1
fi

if test -z $length ;then
   length=3000
fi

# get the sample and labels list. 
samples_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $KASseq)
labels_list=$(awk '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] " ";print ""}}' $labels)

rm -f ${prefix}.labels_basename.txt

# get the absolute path of 'KAS-pipe2 metageneprofile' shell script.
SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

if [[ $regions == "genebody" ]] ;then
   if test -z $peaks && test -z $peakslist ;then
      # generate the matrix file on genebody.
      echo "Compute (sp)KAS-seq matrix of $samples_list on ${assemblyid}_Refseq_${regions} ..."
      echo ""
      computeMatrix scale-regions -R ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.bed -S $samples_list -b $length -a $length --regionBodyLength 6000 --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz 
      echo "done."
      echo ""

      # generate the metagene profile plot.
      echo "Generate the png format metagene profile plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' ..."
      echo ""
      plotProfile -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz -out ${prefix}_on_${assemblyid}_Refseq_${regions}_profile.png --samplesLabel $labels_list --colors $colors_list --plotFileFormat png  --perGroup --plotHeight 7 --plotWidth 9
      echo "done."
      echo ""

      echo "Generate the svg format metagene profile plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' ..."
      echo ""
      plotProfile -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz -out ${prefix}_on_${assemblyid}_Refseq_${regions}_profile.svg --samplesLabel $labels_list --colors $colors_list --plotFileFormat svg  --perGroup --plotHeight 7 --plotWidth 9
      echo "done."
      echo ""
      # rm -f ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz

   elif test -n "$peaks" ;then
      # get the basename of peak file.
      peaks_basename=$( basename ${peaks} .bed )

      # generate the matrix file on genebody.
      echo "Compute (sp)KAS-seq matrix of $samples_list on user provided regions $peaks ..."
      echo ""
      computeMatrix scale-regions -R ${peaks} -S $samples_list -b $length -a $length --regionBodyLength 6000 --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${peaks_basename}.matrix.gz
      echo "done."
      echo ""

      # generate the metagene profile plot.
      echo "Generate the png format metagene profile plot with '${prefix}_on_${peaks_basename}.matrix.gz' ..."
      echo ""
      plotProfile -m ${prefix}_on_${peaks_basename}.matrix.gz -out ${prefix}_on_${peaks_basename}_profile.png --samplesLabel $labels_list --colors $colors_list --plotFileFormat png  --perGroup --plotHeight 7 --plotWidth 9
      echo "done."
      echo ""
  
      echo "Generate the svg format metagene profile plot with '${prefix}_on_${peaks_basename}.matrix.gz' ..."
      echo ""
      plotProfile -m ${prefix}_on_${peaks_basename}.matrix.gz -out ${prefix}_on_${peaks_basename}_profile.svg --samplesLabel $labels_list --colors $colors_list --plotFileFormat svg  --perGroup --plotHeight 7 --plotWidth 9
      echo "done."
      echo ""

   elif test -n "$peakslist" ;then
      # generate the matrix file on genebody	      
      peakslist_basename=$( basename ${peakslist} .bed )
      
      echo "Compute (sp)KAS-seq matrix of $samples_list on $peaks_list ..."
      echo ""
      computeMatrix scale-regions -R $peaks_list -S $samples_list -b $length -a $length --regionBodyLength 6000 --skipZeros --samplesLabel KAS-seq --missingDataAsZero -p $threads -o ${prefix}_on_${peakslist_basename}.matrix.gz
      echo "done."
      echo ""

      # generate the metagene profile plot.
      echo "Generate the png format metagene profile plot with '${prefix}_on_${peakslist_basename}.matrix.gz' ..."
      echo ""
      plotProfile -m ${prefix}_on_${peakslist_basename}.matrix.gz -out ${prefix}_on_${peakslist_basename}.matrix_profile.png --regionsLabel $labels_list --colors $colors_list --plotFileFormat png --plotHeight 7 --plotWidth 9
      echo "done."
      echo ""

      echo "Generate the svg format metagene profile plot with '${prefix}_on_${peakslist_basename}.matrix.gz' ..."
      echo ""
      plotProfile -m ${prefix}_on_${peakslist_basename}.matrix.gz -out ${prefix}_on_${peakslist_basename}.matrix_profile.svg --regionsLabel $labels_list --colors $colors_list --plotFileFormat svg --plotHeight 7 --plotWidth 9
      echo "done."
      echo ""      

   fi
      	      
elif [[ $regions == "TSS" ]] ;then

   # generate the matrix file on TSS.
   echo "Compute (sp)KAS-seq matrix of $samples_list on ${assemblyid}_Refseq_${regions} ..."
   echo ""
   computeMatrix reference-point --referencePoint TSS -R ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.bed -S $samples_list -b $length -a $length --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz
   echo "done."
   echo ""

   # generate the metagene profile plot.
   echo "Generate the png format metagene profile plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' ..."
   echo ""
   plotProfile -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz -out ${prefix}_on_${assemblyid}_Refseq_${regions}_profile.png --samplesLabel $labels_list --colors $colors_list --plotFileFormat png  --perGroup --plotHeight 7 --plotWidth 9
   echo "done."
   echo ""

   echo "Generate the svg format metagene profile plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' ..."
   echo ""
   plotProfile -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz -out ${prefix}_on_${assemblyid}_Refseq_${regions}_profile.svg --samplesLabel $labels_list --colors $colors_list --plotFileFormat svg  --perGroup --plotHeight 7 --plotWidth 9
   echo "done."
   echo ""

   # rm -f ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz

elif [[ $regions == "TES" ]] ;then

   # generate the matrix file on TES.
   echo "Compute (sp)KAS-seq matrix of $samples_list on ${assemblyid}_Refseq_${regions} ..."
   echo ""
   computeMatrix reference-point --referencePoint TES -R ${SH_SCRIPT_DIR}/../annotation/${assemblyid}/${assemblyid}_Refseq.bed -S $samples_list -b $length -a $length --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz
   echo "done."
   echo ""

   # generate the metagene profile plot.
   echo "Generate the png format metagene profile plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' ..."
   echo ""
   plotProfile -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz -out ${prefix}_on_${assemblyid}_Refseq_${regions}_profile.png --samplesLabel $labels_list --colors $colors_list --plotFileFormat png  --perGroup --plotHeight 7 --plotWidth 9
   echo "done."
   echo ""

   echo "Generate the svg format metagene profile plot with '${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz' ..."
   echo ""
   plotProfile -m ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz -out ${prefix}_on_${assemblyid}_Refseq_${regions}_profile.svg --samplesLabel $labels_list --colors $colors_list --plotFileFormat svg  --perGroup --plotHeight 7 --plotWidth 9 
   echo "done."
   echo ""

   # rm -f ${prefix}_on_${assemblyid}_Refseq_${regions}.matrix.gz

elif [[ $regions == "peaks" ]] ;then
   if test -n "$peaks" ;then
       # get the basename of peak file.
       peaks_basename=$( basename ${peaks} .bed )
       # generate the matrix file on peaks.
       echo "Compute (sp)KAS-seq matrix of $samples_list on $peaks ..."
       echo ""
       computeMatrix reference-point --referencePoint center -R $peaks -S $samples_list -b $length -a $length --skipZeros --samplesLabel $labels_list --missingDataAsZero -p $threads -o ${prefix}_on_${peaks_basename}.matrix.gz   
       echo "done."
       echo ""

       # generate the metagene profile plot.
       echo "Generate the png format metagene profile plot with '${prefix}_on_${peaks_basename}.matrix.gz' ..."
       echo ""
       plotProfile -m ${prefix}_on_${peaks_basename}.matrix.gz -out ${prefix}_on_${peaks_basename}_profile.png --samplesLabel $labels_list --colors $colors_list --plotFileFormat png  --perGroup --plotHeight 7 --plotWidth 9
       echo "done."
       echo ""

       echo "Generate the svg format metagene profile plot with '${prefix}_on_${peaks_basename}.matrix.gz' ..."
       echo ""
       plotProfile -m ${prefix}_on_${peaks_basename}.matrix.gz -out ${prefix}_on_${peaks_basename}_profile.svg --samplesLabel $labels_list --colors $colors_list --plotFileFormat svg  --perGroup --plotHeight 7 --plotWidth 9
       echo "done."
       echo ""
       # rm -f ${prefix}_on_${peakslist_basename}.matrix.gz

   elif test -n "$peakslist" ;then
       peakslist_basename=$( basename ${peakslist} .bed )

       echo "Compute (sp)KAS-seq matrix of $samples_list on $peaks_list ..."
       echo ""
       computeMatrix reference-point --referencePoint center -R $peaks_list -S $samples_list -b $length -a $length --skipZeros --missingDataAsZero -p $threads -o ${prefix}_on_${peakslist_basename}.matrix.gz
       echo "done."
       echo ""
 
       echo "Generate the png format metagene profile plot with '${prefix}_on_${peakslist_basename}.matrix.gz' ..."
       echo ""
       plotProfile -m ${prefix}_on_${peakslist_basename}.matrix.gz -out ${prefix}_on_${peakslist_basename}_profile.png --regionsLabel $labels_list --colors $colors_list --plotFileFormat png --plotHeight 7 --plotWidth 9
       echo "done."
       echo ""

       echo "Generate the svg format metagene profile plot with '${prefix}_on_${peakslist_basename}.matrix.gz' ..."
       echo ""
       plotProfile -m ${prefix}_on_${peakslist_basename}.matrix.gz -out ${prefix}_on_${peakslist_basename}_profile.svg --regionsLabel $labels_list --colors $colors_list --plotFileFormat svg --plotHeight 7 --plotWidth 9
       echo "done."
       echo ""
       # rm -f ${prefix}_on_${peakslist_basename}.matrix.gz

   else 
       echo ""
       echo "Please input the regions files. -p [peaks] or -f [peaks files list]"	
       echo ""
       exit -1
   fi
fi

echo "'KAS-pipe2 metageneprofile' run successfully!"
