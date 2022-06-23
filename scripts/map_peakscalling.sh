#!/bin/bash
# 'KAS-pipe2 peakscalling' was developed by Ruitu Lyu on 12-09-2021.

# Stop on error
set -e

## Read arguments                                                     
usageHelp="Usage: KAS-pipe2 peakscalling [ -h ] [ -m peaks caller ] [ -k KAS-seq ] [ -c Control ] [ -b mode ] [ -o prefix ] [ -p pvalue or qvalue ] [ -g assembly id ]."
exampleHelp="Example: nohup KAS-pipe2 peakscalling -t KAS-seq.rep1.bed,KAS-seq.rep2.bed -c Control_Input.rep1.bed,Control_Input.rep2.bed -o KAS-seq -g hg19 &"
peakscallerHelp="-m [peaks caller]: please input the peaks caller (macs14, macs2) that you want to use for KAS-seq peaks calling. DEFAULT: macs2."
KASseqHelp="-k [KAS-seq]: please input the KAS-seq bed or bam files. e.g. KAS-seq.rep1.bed,KAS-seq.rep2.bed or KAS-seq.rep1.bam,KAS-seq.rep2.bam. REQUIRED."
controlHelp="-c [Control]: please input the KAS-seq control bed or bam files. e.g. KAS-seq_Input.rep1.bed,KAS-seq_Input.rep2.bed or KAS-seq_Input.rep1.bam,KAS-seq_Input.rep2.bam. OPTIONAL."
modeHelp="-b [mode]: specify macs2 to perferm KAS-seq peaks calling with 'broad' or 'sharp' mode. -b option only works for macs2. DEFAULT: broad."
prefixHelp="-o [prefix]: please input the prefix (basename), which will be used to generate the name of 'KAS-pipe2 peakscalling' output files. REQUIRED."
cutoffHelp="-p [pvalue or qvalue]: please input the pvalue or qvalue for KAS-seq peaks calling with macs14 or macs2. DEFAULT: macs14: 1e-7; macs2: 0.01"
assemblyidHelp="-g [assembly id]: please specify the reference genome assembly id of KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED."
# genome size. e.g. human(hs): 2.7e9; mouse(mm): 1.87e9; C.elegans(ce): 9e7; fruitfly(dm): 1.2e8; rat(rn): 2.5e9; zebrafish(danRer): 1e9.
helpHelp="-h: print this help and exit.
Note: This shell script mainly invoke macs14 or macs2 for calling (sp)KAS-seq data peaks, please google their official websites for more information."

# print help function.
printHelpAndExit() {
    echo -e ""
    echo -e "$usageHelp"
    echo -e ""
    echo -e "$exampleHelp"
    echo -e ""
    echo -e "$peakscallerHelp"
    echo -e ""
    echo -e "$KASseqHelp"
    echo -e ""
    echo -e "$controlHelp"
    echo -e ""
    echo -e "$modeHelp"
    echo -e ""
    echo -e "$prefixHelp"
    echo -e ""
    echo -e "$cutoffHelp"
    echo -e ""
    echo -e "$assemblyidHelp"
    echo -e ""
    echo -e "$helpHelp"
    echo -e ""
    exit -1
}

# if no parameters was provided, 'KAS-pipe2 peakscalling' will print the help.
if [[ $# == 1 ]] || [[ $1 == "--help" ]] || [[ $1 == "-help" ]] ;then
    printHelpAndExit
fi

# get the value of options.
while getopts 'hm:k:c:b:o:p:g:' opt; do
    case $opt in
        h) printHelpAndExit 0;;
        m) peakscaller=$OPTARG ;;
        k) KASseq=$OPTARG ;;
        c) control=$OPTARG ;;
	b) mode=$OPTARG ;;
        o) prefix=$OPTARG ;;
        p) cutoff=$OPTARG ;;
        g) assemblyid=$OPTARG ;;
        ?) printHelpAndExit 0;;
    esac
done

# check if macs2 was installed.
if ! type macs2 > /dev/null 2>&1 ;then
   echo "macs2 was not installed or not export to the \$PATH'"
   echo ""
   echo "Install macs2 with 'conda install -c bioconda macs2' or refer the official website of 'macs2'."
   echo ""
   exit 1
fi

# Required options.
if test -z $KASseq ;then
   echo ""
   echo "Please input the KAS-seq file. e.g. KAS-seq.bed or KAS-seq.bam. -t [KAS-seq]"
   echo ""
   exit -1
fi

if test -z $prefix ;then
   echo ""
   echo "Please input the prefix (basename) of peaks files. -o [prefix]"
   echo ""
   exit -1
fi

if test -z $assemblyid ;then
   echo ""	
   echo "Please specify the reference genome assembly id of KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. -g [assembly id]"
   echo ""
   exit -1
fi

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
echo "Error: unsupported assembly id: $assemblyid."
echo ""
exit -1
fi


# setup the default options.
if test -z $peakscaller ;then
   peakscaller="macs2"
fi

if test -z $mode && [[ $peakscaller == "macs2" ]] ;then 
    mode="broad" 
fi 


if test -n "$mode" && [[ $mode != "broad" ]] && [[ $mode != "sharp" ]] ;then
    echo ""
    echo "Error: unsupported peak calling mode using macs2: $mode. e.g. broad or sharp. Default: broad."
    echo ""
    exit -1
fi

# test unsuported peak caller.
if [[ $peakscaller != "macs2" ]] && [[ $peakscaller != "macs14" ]] ;then
   echo ""
   echo "Error: unsupported peak caller: $peakscaller. e.g. macs14 or macs2. Default: macs2."
   echo " "
   exit -1
fi

# setup the default pvalue or qvalue.
if test -z $cutoff ;then
   if [[ $peakscaller == "macs2" ]] ;then
      cutoff="0.01"
   elif [[ $peakscaller == "macs14" ]] ;then
      cutoff="1e-7"
   fi	
fi

# call KAS-seq peaks without control using macs2
if test -z $control && [[ $peakscaller == "macs2" ]] && [[ $mode == "broad" ]] ;then	

   echo $KASseq > .KASseq.txt
   sed -i "s/\,/ /g" .KASseq.txt
   KASseq_list=$(cat .KASseq.txt)
   echo "Call (sp)KAS-seq peaks for $KASseq_list at broad mode using macs2 ..." 
   echo ""
   macs2 callpeak -t $KASseq_list -n $prefix --broad -g $genomesize --broad-cutoff $cutoff -q $cutoff 2> ${prefix}_output_${peakscaller}.log
   rm -f .KASseq.txt
   echo "done."
   echo ""

elif test -z $control && [[ $peakscaller == "macs2" ]] && [[ $mode == "sharp" ]] ;then
   echo $KASseq > .KASseq.txt
   sed -i "s/\,/ /g" .KASseq.txt
   KASseq_list=$(cat .KASseq.txt)
   echo "Call (sp)KAS-seq peaks for $KASseq_list at sharp mode using macs2 ..." 
   echo ""
   macs2 callpeak -t $KASseq_list -n $prefix -g $genomesize -q $cutoff 2> ${prefix}_output_${peakscaller}.log
   rm -f .KASseq.txt
   echo "done."
   echo ""

# call KAS-seq peaks with control using macs2
elif test -n "$control" && [[ $peakscaller == "macs2" ]] && [[ $mode == "broad" ]] ;then
   echo $KASseq > .KASseq.txt
   echo $control > .control.txt
   sed -i "s/\,/ /g" .KASseq.txt
   sed -i "s/\,/ /g" .control.txt
   KASseq_list=$(cat .KASseq.txt)
   control_list=$(cat .control.txt)
   echo "Call (sp)KAS-seq peaks for $KASseq_list compaired to $control_list at broad mode using macs2 ..." 
   echo ""
   macs2 callpeak -t $KASseq_list -c $control_list -n $prefix --broad -g $genomesize --broad-cutoff $cutoff -q $cutoff 2> ${prefix}_output_${peakscaller}.log
   rm -f .KASseq.txt
   rm -f .control.txt
   echo "done."
   echo ""

elif test -n "$control" && [[ $peakscaller == "macs2" ]] && [[ $mode == "sharp" ]] ;then	
   echo $KASseq > .KASseq.txt
   echo $control > .control.txt
   sed -i "s/\,/ /g" .KASseq.txt
   sed -i "s/\,/ /g" .control.txt
   KASseq_list=$(cat .KASseq.txt)
   control_list=$(cat .control.txt)
   echo "Call (sp)KAS-seq peaks for $KASseq_list compared to $control_list at sharp mode using macs2 ..."
   echo ""
   macs2 callpeak -t $KASseq_list -c $control_list -n $prefix -g $genomesize -q $cutoff 2> ${prefix}_output_${peakscaller}.log
   rm -f .KASseq.txt
   rm -f .control.txt
   echo "done."
   echo ""

# call KAS-seq peaks without control using macs14
elif test -z $control && [[ $peakscaller == "macs14" ]] ;then
   echo "Call (sp)KAS-seq peaks for $KASseq_list using macs2 ..."
   echo ""
   macs14 -t $KASseq -n $prefix -p $cutoff 2> ${prefix}_output_${peakscaller}.log
   echo "done."
   echo ""

# call KAS-seq peaks with control using macs14	
elif test -n "$control" && [[ $peakscaller == "macs14" ]] ;then
   echo "Call (sp)KAS-seq peaks for $KASseq compared to $control using macs2 ..."
   echo ""
   macs14 -t $KASseq -c $control -n $prefix -p $cutoff 2> ${prefix}_output_${peakscaller}.log
   echo "done."
   echo ""

else
printHelpAndExit 0

fi

echo "'KAS-pipe2 peakscalling' run successfully!"
