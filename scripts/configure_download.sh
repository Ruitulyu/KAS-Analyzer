#!/bin/bash
# 'KAS-pipe2 download' was developped by Ruitu Lyu on 12-01-2021.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-pipe2 download [ -l ] [ -h ] [ -a aligner ] [ -g assembly id ] [ -d directory to save index of aligner ]"
exampleHelp="Example: KAS-pipe2 download -a bowtie2 -g hg19 -d /Software/reference_genome/ "
listHelp="-l list all of the available aligner index for reference genomes."
alignerHelp="-a [aligner]: aligner name you want to use, e.g. bowtie, bowtie2 or bwa. REQUIRED."
assemblyidHelp="-g [assembly id]: reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9,mm10,mm39; Fruitfly: dm3, dm6; Rat: rn6, rn7; C.elegans: ce10, ce11; Zebra fish: danRer10, danRer11. REQUIRED."
dirHelp="-d [directory to save index of aligner]: directory to save downloaded reference genome index. REQUIRED."
helpHelp="-h\-help: print this help and exit."


printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "" 
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$listHelp"
    echo -e ""
    echo -e "$alignerHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$dirHelp"
    echo -e ""
    echo "$helpHelp"
    echo -e ""
    exit -1
}

list() {
    echo "
    aligner:
    bowtie2 or bwa.
    
    Human:
    hg18, hg19, hg38.
    
    Mouse:
    mm9, mm10, mm39.
    
    Fruit fly:
    dm3, dm6.
    
    Rat:
    rn6, rn7.
    
    C. elegans:
    ce10, ce11.
    
    Zebra fish:
    danRer10, danRer11.
    "
    exit -1
}

# if no parameters or '--help' was provided, KAS-pipe2 download will print the help.
if [[ $# == 0 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'lha:g:d:' opt; do
    case $opt in
	l) list 0;;
	h) printHelpAndExit 0;;
	a) aligner=$OPTARG ;;
	g) assemblyid=$OPTARG ;;
	d) indexdir=$OPTARG ;;
	?) printHelpAndExit ;;
    esac
done

# Required options.
if test -z $aligner ;then
   echo "please input aligner: bowtie, bowtie2 or bwa. -a [aligner]"
   echo ""
   printHelpAndExit
fi 

if test -z $assemblyid ;then
   echo "please specify the reference assembly id: Human: hg18, hg19, hg38; Mouse: mm9,mm10,mm39; Fruitfly: dm3, dm6; Rat: rn6, rn7; C.elegans: ce10, ce11; Zebra fish: danRer10, danRer11. -g [assemblyid] "
   echo ""
   printHelpAndExit
fi

if test -z $indexdir ;then
   echo "please input directory to save downloaded reference genome index. -d [indexdir]"
   echo ""
   printHelpAndExit
fi

# select the right download link.
if [ "$aligner" == "bowtie2" ]; then
   # download human bowtie2 index.
   if [[ "$assemblyid" == "hg18" ]] ;then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/hg18_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "hg19" ]] ;then      
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/hg19_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "hg38" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/hg38_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   # download mouse bowtie2 index.   
   elif [[ "$assemblyid" == "mm9" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/mm9_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "mm10" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/mm10_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "mm39" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/mm39_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   # download fruit fly bowtie2 index.
   elif [[ "$assemblyid" == "dm3" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/dm3_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "dm6" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/dm6_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   # download rat bowtie2 index.
   elif [[ "$assemblyid" == "rn6" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/rn6_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "dn7" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/rn7_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   # download C. elegans bowtie index.
   elif [[ "$assemblyid" == "ce10" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/ce10_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "ce11" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/ce11_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   # download zebra fish bowtie2 index.
   elif [[ "$assemblyid" == "danRer10" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/danRer10_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "danRer11" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/danRer11_Bowtie2Index/
      echo "$assemblyid $aligner index download successfully!"

   # test the unavailable assembly id.
   else
      echo ""	   
      echo "Unsupported assembly id: $assemblyid. Please specify the reference assembly id: Human: hg18, hg19, hg38; Mouse: mm9,mm10,mm39; Fruitfly: dm3, dm6; Rat: rn6, rn7; C.elegans: ce10, ce11; Zebra fish: danRer10, danRer11. -s [assembly id]"
      echo ""
      printHelpAndExit
   fi

elif [ "$aligner" == "bwa" ]; then

   # download human bwa index.
   if [[ "$assemblyid" == "hg18" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/hg18_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "hg19" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/hg19_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "hg38" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/hg38_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   # download mouse bwa index.
   elif [[ "$assemblyid" == "mm9" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/mm9_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "mm10" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/mm10_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "mm39" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/mm39_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   # download fruit fly bwa index.
   elif [[ "$assemblyid" == "dm3" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/dm3_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "dm6" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/dm6_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   # download rat bwa index.
   elif [[ "$assemblyid" == "rn6" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/rn6_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "dn7" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/rn7_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   # download C. elegans bwa index.
   elif [[ "$assemblyid" == "ce10" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/ce10_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "ce11" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/ce11_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   # download zebra fish bwa index.
   elif [[ "$assemblyid" == "danRer10" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/danRer10_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   elif [[ "$assemblyid" == "danRer11" ]]; then
      wget -r -np -nH --cut-dirs=1 --reject "index.html*" -P $indexdir http://128.135.119.125/helab_genome_index/danRer11_BWAIndex/
      echo "$assemblyid $aligner index download successfully!"

   # test the unavailable assembly id.
   else
      echo ""	   
      echo "Unsupported assembly id: $assemblyid. Please specify the reference assembly id: Human: hg18, hg19, hg38; Mouse: mm9,mm10,mm39; Fruitfly: dm3, dm6; Rat: rn6, rn7; C.elegans: ce10, ce11; Zebra fish: danRer10, danRer11. -g [assembly id]"
      echo ""
      printHelpAndExit
   fi
else
   echo ""	
   echo "Unsupported aligner: $aligner. e.g. bowtie2 or bwa. -a [aligner]"
   echo ""
   printHelpAndExit

fi
