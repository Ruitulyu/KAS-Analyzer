# KAS-pipe2: a user-friendly toolkit for exploring KAS-seq data

## Author: Ruitu Lyu, Chemistry department, the University of Chicago 
---------------------------------------------------------------------------
<script type="text/javascript" src="//rf.revolvermaps.com/0/0/6.js?i=5f2zpl0mkwd&amp;m=7&amp;c=e63100&amp;cr1=ffffff&amp;f=arial&amp;l=0&amp;bv=90&amp;lx=-420&amp;ly=420&amp;hi=20&amp;he=7&amp;hc=a8ddff&amp;rs=80" async="async"></script>
![Image](https://raw.githubusercontent.com/Ruitulyu/KAS-pipe2/main/image/KAS-pipe2.jpg)

## Introduction

KAS-pipe2 is a collection of command line tools specifically developped for exploring KAS-seq or strand-specific (sp)KAS-seq data, including the basic processing tools for quality control, reference genome index, raw reads mapping, and heatmaps and summary plots. KAS-pipe2 also includes many novel features and completely new frame design compared to KAS-pipe. e.g. time-courese(TC) KAS-seq differential analysis, R-loop identification (only for spKAS-seq), ss-enhancers, motif, index calculation, termination length and so on.

KAS-pipe2 is still on active development and the source code is hosted on [GitHub](https://github.com/Ruitulyu/KAS-pipe2).

## Installation

**Install by cloning KAS-pipe2 git repository on github:**

You can install KAS-pipe2 on command line (linux/mac) by cloning git repository on github:

	$ git clone https://github.com/Ruitulyu/KAS-pipe2.git
	$ cd KAS-pipe2
	$ bash ./setup.sh
	
	# If anaconda or miniconda was not installed on your system.
	$ KAS-pipe2 install -conda
	
	# Install conda 'KAS-pipe2' environment. 
	$ KAS-pipe2 install -KAS-pipe2
	
	# Activate conda 'KAS-pipe2' environment.
	$ KAS-pipe2 activate
	
## Overview
---------------------------------------------------------------------------

KAS-pipe2 is a collection of command line tools for KAS-seq or strand specific KAS-seq(spKAS-seq) data analysis.
```
Version:   v2.0
About:     KAS-pipe2 is mainly developed by Ruitu Lyu, postdoc fellow in Prof. Chuan He's group at the University of Chicago.
Docs:      https://ruitulyu.github.io/KAS-pipe2/
Code:      https://github.com/Ruitulyu/KAS-pipe2
Mail:      https://github.com/Ruitulyu/KAS-pipe2/discussions
```
### Usage:
KAS-pipe2 sub-command [options]

### The KAS-pipe2 sub-commands include:
----------------------------------------
#### Configure
```	
   download        Downlaod the index of reference genome for aligners (bowtie, bowtie2, bwa, star).
   build           Build the index of reference genome for aligners (bowtie, bowtie2, bwa, star).
   install         Install and check the 'KAS-pipe2' conda environment; check the installation of tools.
   uninstall       Uninstall the 'KAS-pipe2' conda environment.
   activate        Activate the 'KAS-pipe2' conda environment.
   deactivate      Deactivate the 'KAS-pipe2' conda environment.
```	
#### Fastqc
```
   fastqc          Generate basic quality control metrics for KAS-seq data.
   readsnum        Calcuate the reads number of raw fastq files.
   statistics      Calculate the mapping statistics for KAS-seq data.
   FRiP            Calculate the fraction of reads in peaks.
   fragmentsize    Measure the fragment size of paired-end KAS-seq data.
   correlation     Calculate the correlation coefficient and pvalue, generate correlation plot for replicates of KAS-seq data.
   saturation      Perform saturation analysis for KAS-seq data.
   complexity      Calculate the complexity metric for (sp)KAS-seq data, including the PCR Bottlenecking Coefficient and Non-Redundant Fraction (NRF).
   genomicdist     Visualize the genomic distribution for KAS-seq peaks (table and plot).
   fingerprint     Plot fingerprint for KAS-seq data.
```
#### Mapping
```
   trim            Trim adapter and low quality sequence, perform quality control for raw KAS-seq data.
   KAS-seq         Align KAS-seq data to the reference genome, deduplicate mapped reads, and generate several files with maped reads (bam, bed and bedGraph).
   spKAS-seq       Align strand specific KAS-seq (spKAS-seq) data. Note: we strongly recommend paired-end sequencing for spKAS-seq data to accurately measure the fragment size.
   peakscalling    Call broad or sharp peaks for KAS-seq data.
   normalize       Normalize KAS-seq data with bedGraph density files.
   bedGraphToBigWig Transfer normalized bedGraph file to bigWig file.
   targetgenes     Define target or associated genes (promoter, genebody, terminator or gene) of KAS-seq peaks, R-loops or enhancers loci."
   UCSC            Generate bedGraph files ready for submitting to UCSC genome browser.
```

#### Plot
```
   profile         Generate metagene profile for KAS-seq data (normalized bigWig files are needed).
   heatmap         Generate heatmap for KAS-seq data (normalized bigWig files are needed).
```

#### Differential KAS-seq analysis
```
   KASexpre        Calculate normalized KAS-seq expression levels on promoter, genebody, genes or custom regions. 
   diff            Perform differential KAS-seq analysis on promoter, genebody, gene, bin or custom regions.        
   TC              Perform 'case-only' or 'case-control' differential time course(TC) analysis for (sp)KAS-seq data.
   PCA             Perform and plot PCA analysis for (sp)KAS-seq data.
```
#### R-loops
```
   R-loop          Identify R-loops regions with spKAS-seq data.
```   
  
#### Single-stranded enhancers identification
```
   ss_enhancer     Identify the single stranded (ss) enhancers.
   motif           Identify enriched TF binding motifs on ss_enhancers.
```   
  
#### Termination length
```
   termilength     Calculate the transcription termination length of protein coding genes.
```   
 
#### KAS-seq index
```
   index           Calculate the pausing or termination index.
```   

#### General help
```
   --help          Print this help menu.
   --version       Print the version of KAS-pipe2 you are using.
   --contact       Feature requests, bugs, mailing lists, etc.
```   

## Tutorial
---------------------------------------------------------------------------

### Configure sub-commands:
---------------------------------------------------------------------------
#### download
#### usage: KAS-pipe2 download [ -l ] [ -h ] [ -a aligner ] [ -g assembly id ] [ -d directory to save index of aligner ]
```  
Example: KAS-pipe2 download -a bowtie2 -g hg19 -d /Software/reference_genome/ 

-l list all of the available aligner index for reference genomes.

-a [aligner]: aligner name you want to use, e.g. bowtie, bowtie2 or bwa. REQUIRED.

-g [assembly id]: reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9,mm10,mm39; Fruitfly: dm3, dm6; Rat: rn6, rn7; C.elegans: ce10, ce11; Zebra fish: danRer10, danRer11. REQUIRED.

-d [directory to save index of aligner]: directory to save downloaded reference genome index. REQUIRED.

-h\-help: print this help and exit.
```   

#### build
#### usage: KAS-pipe2 build [ -h ] [ -a aligner ] [ -g genome fasta ] [ -p index prefix ] [ -t threads ] [ -d index dir ]
```  
Example: KAS-pipe2 build -a bowtie2 -g ./genome.fa -p hg19 -t 10 -d /Software/hg19_Bowtie2Index/ 

-a [aligner]: please specify aligner name you want to use, e.g. bowtie, bowtie2 or bwa. REQUIRED.

-g [genome fasta]: please input the path of reference genome fasta file. REQUIRED.

-t [threads]: please specify the number of threads. Default: 1.

-p [index prefix]: please input the prefix (basename) of the aligners' (bowtie, bowtie2 or bwa) reference genome index. Default: basename of fasta file .

-d [index dir]: directory to save newly built reference genome index. REQUIRED.

-h\-help: print this help and exit.
```  

#### install
#### usage: KAS-pipe2 install [ -h\--help ] [ -conda ] [ -check ] [ -t tools ] [ -KAS-pipe2 ]
``` 
Example: KAS-pipe2 install or KAS-pipe2 install -check

-conda: check the installation or install anaconda in your computer.

-check: list the installed and uninstalled tools.

-t [tools]: check the installation of specific tool, if not, will install automaticially.

-h\-help: print this help and exit.
Note: this subcommand is used to install conda environment and specific tool that needed in KAS-pipe2.
``` 

#### uninstall
#### usage: KAS-pipe2 uninstall 
This subcommand is used to uninstall KAS-pipe2 conda environment.

#### activate
#### usage: KAS-pipe2 activate
This subcommand is used to activate KAS-pipe2 conda environment.

#### inactivate
#### usage: KAS-pipe2 deactivate
This subcommand is used to deactivate KAS-pipe2 conda environment.

### Fastqc sub-commands:
---------------------------------------------------------------------------
#### fastqc
#### usage: KAS-pipe2 fastqc [ -h/--help ] [ -t threads ] [ -c contaminants ] [ -o output dir ] [ -k KAS-seq ] 
``` 
Example: nohup KAS-pipe2 fastqc -t 10 -k KAS-seq.rep1.fastq.gz,KAS-seq.rep2.fastq.gz,KAS-seq.rep3.fastq.gz &

-t [threads]: please input number of threads to be used for quality control check. Default: 1.

-c [contaminants]: please specify a file which contains the list of contaminants (format: name[tab]sequence) to screen overrepresented sequences against. Default: no.

-o [output dir]: please specify the output directory with output files.

-k [KAS-seq]: please input the KAS-seq data that you want to know the quality control, like sequencing quality, duplicates, contaminants (adapter sequence).

-h/--help: print this help and exit.
Note: The 'KAS-pipe2 fastqc' shell script is applied to check quality control and identify a potential type of problem in your KAS-seq data in non-interactive mode. It mainly invoke FASTQC, please refer to the FASTQC official website for more information.
``` 

#### readsnum
#### usage: KAS-pipe2 readsnum [ -h/--help ] [ -o prefix ] [ -f format ] 
``` 
Example: nohup KAS-pipe2 readsnum -p KAS-seq_reads_num -f fastq.gz &

-o [prefix]: please specify the prefix (basename) of 'KAS-pipe2 readsnum' output files. REQUIRED.

-f [format]: please specify the format of raw reads data. e.g. fastq, fq, fastq.gz, fasta, fa or fa.gz. REQUIRED.

-h/--help: print this help and exit.
Note: The 'KAS-pipe2 readsnum' shell script is applied to calculate the reads number of raw sequencing files.
``` 

#### readsnum
#### usage: KAS-pipe2 statistics [ -h/--help ] [ -o prefix ] [ -l labels ] [ -s summary folder ]
```
Example: nohup KAS-pipe2 statistics -o KAS-seq_statistics -l labels.txt -s summary.txt &

-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 statistics' output files. REQUIRED.

-l [labels]: please input the txt file containing labels of (sp)KAS-seq summary file generated by 'KAS-pipe2 (sp)KAS-seq'. Default: basename of summary file.
KAS-seq.rep1
KAS-seq.rep2
KAS-seq.rep3
KAS-seq.rep4        ---labels.txt

-s [summary folder or summary.txt]: please input the absolute path of folder with summary files or the txt file containing the summary files. e.g. -s summary.txt or /absolute path of summary folder. REQUIRED.
WT.rep1.KAS-seq_mapping_summary.txt                     or  /absolute path/Summary/
WT.rep2.KAS-seq_mapping_summary.txt
KO.rep1.KAS-seq_mapping_summary.txt
KO.rep2.KAS-seq_mapping_summary.txt     ---summary.txt

-h/--help: print this help and exit.
Note: The 'KAS-pipe2 statistics' shell script is applied to generate the table containing (sp)KAS-seq mapping statistics.
```

#### FRiP

#### Usage: KAS-pipe2 FRiP [ -h/--help ] [ -o prefix ] [ -p peaks ] [ -l labels ] [ -k KAS-seq ]
```
Example: nohup KAS-pipe2 FRiP -o KAS-seq_FRiP -p peaks_files.txt -l labels.txt -k KAS-seq.txt &

-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 FRiP' output files. REQUIRED.

-p [peaks]: please input the text file containing the peaks files. REQUIRED.
Example:
KAS-seq_WT_rep1_peaks.bed
KAS-seq_WT_rep2_peaks.bed
KAS-seq_KO_rep1_peaks.bed
KAS-seq_KO_rep2_peaks.bed      ---peaks_files.txt

-l [labels]: please input the text file containing the labels of KAS-seq or spKAS-seq data that show in fraction of reads in peaks (FRiP) plot. Default: basename of KAS-seq files.
Example:
WT_rep1
WT.rep2
KO.rep1
KO.rep2                        ---labels.txt

-k [KAS-seq]: please input the text file containing bed files (uniquely mapped reads used for 'KAS-pipe2 peakcalling'), which are used to calcuate fraction of reads in peaks (FRiP) score. The order and number of (sp)KAS-seq data should be the consistent with the labels file. REQUIRED.
Example:
KAS-seq_WT_rep1.bed
KAS-seq_WT_rep2.bed
KAS-seq_KO_rep1.bed
KAS-seq_KO_rep2.bed            ---KAS-seq.txt

-h/--help: print this help and exit.
Note: The 'KAS-pipe2 FRiP' shell script is applied to calculate and plot fraction of reads in peaks (FRiP) scores.
```

#### fragmentsize
#### Usage: KAS-pipe2 fragmentsize [ -h/--help ] [ -o prefix ] [ -l labels ] [ -k KAS-seq ]
```
Example: nohup KAS-pipe2 fragmentsize -o KAS-seq_fragmentsize -l labels.txt -k KAS-seq.txt &

-o [prefix]: please specify the prefix (basename) of 'KAS-pipe2 fragmentsize' output files. REQUIRED.

-l [labels.txt]: please input the text file containing the labels of paired-end KAS-seq or spKAS-seq data that show in fragment size plot. Default: basename of (sp)KAS-seq bed files.
Example:
WT.rep1
WT.rep2
WT.rep3
WT.rep4                        ---labels.txt

-k [KAS-seq.txt]: please input the text file containing bed files (uniquely mapped reads from 'KAS-pipe2 (sp)KAS-seq'), which are used to calcuate fragment size of DNA fragments. The order and number of (sp)KAS-seq bed files should be the consistent with the labels in labels.txt file. REQUIRED.
Example:
KAS-seq_WT_PE.rep1.bed
KAS-seq_WT_PE.rep2.bed
KAS-seq_WT_PE.rep3.bed
KAS-seq_WT_PE.rep4.bed         ---KAS-seq.txt

-h/--help: print this help and exit.
Note: The 'KAS-pipe2 fragmentsize' shell script is applied to calculate and plot fragment size of (sp)KAS-seq data. Note: it only works for paired-end (sp)KAS-seq data.
```

### Map sub-commands:
---------------------------------------------------------------------------

#### trim
#### Usage: KAS-pipe2 trim [ -h ] [ -a adapter ] [ -t threads ] [ -f ] [ -q quality ] [ -l length ] [ -1 read1 ] [ -2 read2 ]
```
Example:
       Single-end:
       nohup KAS-pipe2 trim -a illumina -t 10 -1 KAS-seq.fastq.gz &
       Paired-end:
       nohup KAS-pipe2 trim -a illumina -t 10 -1 KAS-seq.R1.fastq.gz -2 KAS-seq.R2.fastq.gz &

-a [adapter types]: adapter sequence to be trimmed. e.g. illumina, nextera or small_rna. Hint: most of the NGS data used Illumina adapter. If not specified explicitly. KAS-pipe2 trim will auto-detect.

-t [threads]: number of threads to be used for trimming. Default: 1.

-f: instruct 'KAS-pipe2 trim' to check quality control before trimming. Default: off.

-q [quality]: trim low-quality ends from reads in addition to adapter removal. Default Phred score(ASCII+33): 20.

-l [length]: discard reads that became shorter than length INT bp because of either quality or adapter trimming. Default: 30.

-1 [read1]: please input single-end KAS-seq raw fastq file or read 1 of paired-end KAS-seq raw fastq files. REQUIRED.

-2 [read2]: please input read2 of paired-end KAS-seq raw fastq files.

-h: print this help and exit.
Note: The 'KAS-pipe2 trim' shell script mainly invoke the trim-galore, please refer to http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/ for more information.
```

#### KAS-seq 
#### Usage: KAS-pipe2 KAS-seq [ -h ] [ -a aligner ] [ -t threads ] [ -i index path ] [ -u ] [ -e extend length ] [ -o prefix ] [ -s assembly id ] [ -1 read1 ] [ -2 read2 ]
```
Example:
       Single-end:
       nohup KAS-pipe2 KAS-seq -a bowtie2 -t 10 -i /absolute path/hg19_Bowtie2Index/hg19 -e 150 -o KAS-seq -s hg19 -1 KAS-seq.trim.fastq.gz &
       Paired-end:
       nohup KAS-pipe2 KAS-seq -a bowtie2 -t 10 -i /absolute path/hg19_Bowtie2Index/hg19 -o KAS-seq -s hg19 -1 KAS-seq.trim.R1.fastq.gz -2 KAS-seq.trim.R2.fastq.gz &
       Note: Bowtie2 Index example: /absolute path/Bowtie2Index/hg19; bwa Index example: /absolute path/BWAIndex/hg19.fa

-a [aligner]: please specify the aligner (bowtie2 or bwa) you want to use to map KAS-seq data. DEFAULT: bowtie2.

-t [threads]: please input the number of threads used for KAS-seq data mapping. DEFAULT: 1.

-i [index path]: please input the absolute path of reference genome index for aligner. REQUIRED.

-u: please specify to filter the unique mapped reads. DEFAULT: off.

-e [extendlengthHelp]: please input the extend length for single-end KAS-seq data. DEFAULT: 150.

-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 KAS-seq' output files. REQUIRED.

-s [assembly id]: please specify the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. REQUIRED.

-1 [read1]: please input trimmed single-end KAS-seq fastq file or read 1 of paired-end KAS-seq raw fastq files; Compressed .fastq.gz is accepted. REQUIRED.

-2 [read2]: please input trimmed read2 of paired-end KAS-seq raw fastq files. Compressed .fastq.gz is accepted.

-h: print this help and exit.
Note: The 'KAS-pipe2 KAS-seq' shell script mainly invoke the specific aligner (bowtie, bowtie2 or bwa) for KAS-seq data mapping, please refer to their official websites for more information.
```

#### spKAS-seq
#### Usage: KAS-pipe2 spKAS-seq [ -h ] [ -t threads ] [ -i index path ] [ -u ] [ -r ] [ -f fold change ] [ -b bin size ] [ -e extend length ] [ -o prefix ] [ -s assembly id ] [ -1 read1 ] [ -2 read2 ]
```
Note: we strongly recommend paired-end sequencing for strand specific KAS-seq (spKAS-seq) data to accurately measure the fragments size.

Example:
       Single-end:
       nohup KAS-pipe2 spKAS-seq -t 10 -i /absolute path/hg19_Bowtie2Index/hg19 -o spKAS-seq -r -s hg19 -1 spKAS-seq.trim.R1.fastq.gz &
       Paired-end:
       nohup KAS-pipe2 spKAS-seq -t 10 -i /absolute path/hg19_Bowtie2Index/hg19 -o spKAS-seq -r -s hg19 -1 spKAS-seq.trim.R1.fastq.gz -2 spKAS-seq.trim.R2.fastq.gz &
       Note: Bowtie2 Index example: /absolute path/Bowtie2Index/hg19.

-t [threads]: please specify the number of threads used for spKAS-seq data mapping. DEFAULT: 1.

-i [index path]: please input the absolute path of reference genome index for aligner. REQUIRED.

-u: please specify to filter the unique mapped reads. DEFAULT: off.

-r: please specify if identify R-loops regions with spKAS-seq data. DEFAULT: off.

-f [fold change cutoff]: please specify the fold change cutoff of spKAS-seq reads difference between plus and minus strands used for R-loops identification. DEFAULT: 2.

-b [bin size]: please specify the size of bins used to identify R-loops. DEFAULT: 500.

-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 spKAS-seq' output files. REQUIRED.

-s [assembly id]: please input the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. REQUIRED.

-1 [read1]: please input trimmed single-end spKAS-seq fastq file or read1 of paired-end spKAS-seq fastq files; compressed .fastq.gz is accepted. REQUIRED.

-2 [read2]: please input trimmed read2 of paired-end spKAS-seq raw fastq files. compressed fastq.gz is accepted.

-h: print this help and exit.
Note: The 'KAS-pipe2 spKAS-seq' shell script mainly invoke bowtie2 for spKAS-seq data mapping and R-loops identification, please refer to their official websites for more information.
```

#### peakscalling
#### Usage: KAS-pipe2 peakscalling [ -h ] [ -m peaks caller ] [ -t KAS-seq ] [ -c Control ] [ -b mode ] [ -o prefix ] [ -p pvalue or qvalue ] [ -g assembly id ].
```
Example: nohup KAS-pipe2 peakscalling -t KAS-seq.rep1.bed,KAS-seq.rep2.bed -c Control_Input.rep1.bed,Control_Input.rep2.bed -o KAS-seq -g hg19 &

-m [peaks caller]: please input the peaks caller (macs14, macs2) that you want to use for KAS-seq peaks calling. DEFAULT: macs2.

-t [KAS-seq]: please input the KAS-seq bed or bam files. e.g. KAS-seq.rep1.bed,KAS-seq.rep2.bed or KAS-seq.rep1.bam,KAS-seq.rep2.bam. REQUIRED.

-c [Control]: please input the KAS-seq control bed or bam files. e.g. KAS-seq_Input.rep1.bed,KAS-seq_Input.rep2.bed or KAS-seq_Input.rep1.bam,KAS-seq_Input.rep2.bam. OPTIONAL.

-b [mode]: specify macs2 to perferm KAS-seq peaks calling with 'broad' or 'sharp' mode. -b option only works for macs2. DEFAULT: broad.

-o [prefix]: please input the prefix (basename), which will be used to generate the name of 'KAS-pipe2 peakscalling' output files. REQUIRED.

-p [pvalue or qvalue]: please input the pvalue or qvalue for KAS-seq peaks calling with macs14 or macs2. Default: macs14: 1e-7; macs2: 0.01

-g [assembly id]: please specify the reference genome assembly id of KAS-seq data. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED.

-h: print this help and exit.
Note: This shell script mainly invoke macs14 or macs2 for calling (sp)KAS-seq data peaks, please google their official websites for more information.
```

#### normalize
#### Usage: KAS-pipe2 normalize [ -h/--help ] [ -k KAS-seq ] [ -r ratios ]
```
Example: nohup KAS-pipe2 normalize -k KAS-seq_data.txt -r ratios.txt &

-k [KAS-seq_data.txt]: please input the text file containing the bedGraph files generated from 'KAS-pipe2 (sp)KAS-seq'. REQUIRED.
Example:
KAS-seq_WT.rep1.bg
KAS-seq_WT.rep2.bg
KAS-seq_KO.rep1.bg 
KAS-seq_KO.rep2.bg   ---KAS-seq_data.txt

-r [ratios.txt]: please input the text file containing ratios that used to normalize KAS-seq data, which can be calculated based on mapped reads number or SpikeIn reads. The order and number of ratios should be the consistent with KAS-seq bedGraph files. REQUIRED.
Example:
1.10
1.20
1.30
1.23                 ---ratios.txt

-h/--help: print this help and exit.
Note: The 'KAS-pipe2 normalize' shell script is applied to normalize spKAS-seq or KAS-seq data.
```

#### bedGraphToBigWig
#### Usage: KAS-pipe2 bedGraphToBigWig [ -h/--help ] [ -k KAS-seq ] [ -s assembly id ]
```
Example: nohup KAS-pipe2 bedGraphToBigWig -k KAS-seq_data.txt -s hg19 &

-k [KAS-seq]: please input the text file containing the bedGraph files generated from 'KAS-pipe2 KAS-seq'. REQUIRED.
Example:
KAS-seq_WT.rep1.nor.bg
KAS-seq_WT.rep2.nor.bg
KAS-seq_KO.rep1.nor.bg 
KAS-seq_KO.rep2.nor.bg    ---KAS-seq_data.txt

-s [assembly id]: please input the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the normalized KAS-seq bedGraph files. REQUIRED.

-h/--help: print this help and exit.
Note: The 'KAS-pipe2 bedGraphToBigWig' shell script is applied to transfer (sp)KAS-seq bedGraph to bigWig files.
```

#### targetgenes
#### Usage: KAS-pipe2 targetgenes [ -h/--help ] [ -o prefix ] [ -s assembly id ] [ -f features ] [ -l length ] [ -p peaks ]
```
Example: 
Gene features:
nohup KAS-pipe2 targetgenes -o KAS-seq_peaks_target_genes -s mm10 -f promoter -p KAS-seq_peaks.bed &
Associated genes of enhancers:
nohup KAS-pipe2 targetgenes -o KAS-seq_ss_enhancers_asso_genes -s mm10 -f enhancer -l 50000 -p KAS-seq_ss_enhancers.bed &

-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 targetgenes' output files. DEFAULT: basename of peaks file.

-s [assembly id]: please specify the genome assembly id of KAS-seq peaks. -s [assembly id]. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED.

-f [features]: please specify the gene feagures used to define target genes. e.g. promoter, genebody, terminator, gene or enhancer. REQUIRED.

-l [length]: please specify the length cutoff (length to TSS) to define enhancer's associated genes. DEFAULT: 50000.

-p [peaks]: please specify the KAS-seq peaks, R-loops or enhancer to define their target genes. REQUIRED.

-h/--help: print this help and exit.
Note: The 'KAS-pipe2 targetgenes' shell script is applied to define target or associated genes (promoter, genebody, terminator or gene) of KAS-seq peaks, R-loops or enhancers loci.
```

#### UCSC
#### Usage: KAS-pipe2 UCSC [ -h/--help ] [ -k KAS-seq ] [ -n UCSC track ] [ -c track colors ]
```
Example: nohup KAS-pipe2 UCSC -k KAS-seq_data.txt -n UCSC_track_names.txt &

-k [KAS-seq]: please input the text file containing the bedGraph files generated from 'KAS-pipe2 KAS-seq'. REQUIRED.
Example:
KAS-seq_WT.rep1.nor.bg
KAS-seq_WT.rep2.nor.bg
KAS-seq_KO.rep1.nor.bg 
KAS-seq_KO.rep2.nor.bg    ---KAS-seq_data.txt

-n [UCSC track]: please input the text file containing the track names of KAS-seq or spKAS-seq data that you want to visualize on UCSC genome browser. REQUIRED.
Example: 
KAS-seq_WT.rep1 
KAS-seq_WT.rep2
KAS-seq_KO.rep1
KAS-seq_KO.rep2           ---UCSC_track_names.txt

 -c [track colors]: please input the text file containing the R,G,B colors list for KAS-seq or spKAS-seq data that you want to visualize on UCSC genome browser. Default: black (0,0,0).
Example:
255,102,102
255,178,102
102,255,102
102,178,255               ---track_colors.txt

For more colors options, please refer to https://www.w3schools.com/colors/colors_rgb.asp or https://www.rapidtables.com/web/color/RGB_Color.html.


-h/--help: print this help and exit.
Note: The 'KAS-pipe2 UCSC' shell script is used to generate files for uploading into UCSC genome browser.
```

### Summary Plot sub-commands:
---------------------------------------------------------------------------

#### profile
#### Usage: KAS-pipe2 profile [ -h/--help ] [ -t threads ] [ -s assembly id ] [ -e length ] [ -o prefix ] [ -r regions ] [ -p peaks file ] [ -f peaks files list ] [ -l labels ] [ -c colors ] [ -k KAS-seq ]
```
Example: 
Genomic features(genebody, TSS or TES):
nohup KAS-pipe2 profile -t 10 -s hg19 -o KAS-seq_genebody -r genebody -c red,blue,green,purple -l labels.txt -k KAS-seq.txt &

Custom regions(peaks. e.g. enhancers.bed):
nohup KAS-pipe2 profile -t 10 -o KAS-seq_peaks -r peaks -p KAS-seq_peaks.bed -c red,blue,green,purple -l labels.txt -k KAS-seq.txt &

KAS-seq signal on different regions (-f [peaks list] must be specified):
nohup KAS-pipe2 profile -t 10 -o KAS-seq_different_clusters -r peakslist -f peaks_cluster1.bed,peaks_cluster2.bed,peaks_cluster3.bed -c red,green,purple -l labels.txt -k KAS-seq.txt &

-t [threads]: please input the number of threads used for generating (sp)KAS-seq metagene profile plot. DEFAULT: 1.

-s [assemblyid]: please specify the reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. Note: the assembly id need to be consistent with the reference genome index. REQUIRED only for 'genomic features' mode.

-o [KAS-seq_profile]: please input the prefix (basename), which will be used to generate the name of 'KAS-pipe2 metageneprofile' output files. REQUIRED.

-e [length]: please specify the distance upstream of the start site of the regions defined in the region file. If the regions are genebody, this would be the distance upstream of the transcription start site. DEFAULT: 3000.

-r [regions]: please specify the regions types for generating metagene profile plot. e.g. Genomic features: genebody, TSS or TES; Custom regions: peaks or peakslist. REQUIRED.

-p [peaks file]: please specify the peak file. REQUIRED only for 'custom regions: peaks' mode.

-f [peaks list]: please specify the peak files list. e.g. peaks_cluster1.bed,peaks_cluster2.bed,peaks_cluster3.bed. REQUIRED only for 'custom regions: peakslist' mode. Note: only one KAS-seq bigWig file is needed.

-c [colors]: please specify the color list for (sp)KAS-seq data in metagene profile plot. Note: the number of colors in the profile plot needs to be consistent with the number of KAS-seq bigWig files. REQUIRED. Note: the list of valid color names https://matplotlib.org/examples/color/named_colors.html.

-l [labels.txt]: please input the text file containing the labels of (sp)KAS-seq data or peaks files (-f need to be specified) that used for generating metagene profile. Default: basename of (sp)KAS-seq bigWig files or peaks files.
Example:
WT_rep1          Cluster1
WT.rep2          Cluster2
KO.rep1          Cluster3
KO.rep2   or                        ---labels.txt

-k [KAS-seq.txt]: please input the text file containing normalized (sp)KAS-seq bigWig files, which can be generated with 'KAS-pipe2 normalize' and 'KAS-pipe2 bedGraphToBigWig' shell scripts. The order and number of (sp)KAS-seq bigWig files should be the consistent with the labels file when provided. REQUIRED.
Example:
KAS-seq_WT_rep1.nor.bigWig
KAS-seq_WT_rep2.nor.bigWig
KAS-seq_KO_rep1.nor.bigWig
KAS-seq_KO_rep2.nor.bigWig          ---KAS-seq.txt

-h/--help: print this help and exit.
Note: The 'KAS-pipe2 profile' shell script is applied to generate metagene profile for (sp)KAS-seq data on genomic features( genebody, TSS or TES) or provided custom regions. 'KAS-pipe2 metageneprofile' shell script mainly invoke deeptools 'computeMatrix' and 'plotProfile', please refer to https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html for more information.
```

#### heatmap
#### Usage: KAS-pipe2 heatmap [ -h/--help ] [ -t threads ] [ -e length ] [ -s assembly id ] [ -q ] [ -u samples using ] [ -m maximum value ] [ -o prefix ] [ -r regions ] [ -p peaks ] [ -l labels ] [ -c colors ] [ -k KAS-seq ]
```
Example: 
Genomic features(genebody, TSS or TES):
nohup KAS-pipe2 heatmap -t 10 -s hg19 -o KAS-seq_heatmap -r genebody -q -c Reds -l labels.txt -k KAS-seq.txt &
Custom regions(peaks. e.g. enhancers.bed):
nohup KAS-pipe2 heatmap -t 10 -o KAS-seq_heatmap -r peaks -q -p KAS-seq_peaks.bed -c Reds,Reds,Blues,Blues -l labels.txt -k KAS-seq.txt &

-t [threads]: please specify the number of threads used for generating (sp)KAS-seq heatmap plot. Default: 1.

-s [assemblyid]: please specify the genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED only for 'genomic features' mode.

-q: please specify to sort the regions. Default: off.

-u [sample using]: please specify the samples list to sort the regions. e.g. -u 1,2,3

-e [length]: please specify the distance upstream of the start site of the regions defined in the region file. If the regions are genebody, this would be the distance upstream of the transcription start site. DEFAULT: 3000.

-m [maximum value]: please specify the maximum value of the heatmap intensities. e.g. -m 15,20,60. Note: the number of values need to be consistent with KAS-seq samples.

-o [KAS-seq_heatmap]: please specify the prefix (basename) of 'KAS-pipe2 heatmap' output files. REQUIRED.

-r [regions]: please specify the regions types to generate heatmap plot. e.g. Genomic features: genebody, TSS or TES; Custom regions: peaks. REQUIRED.

-p [peakslist]: please input the peak list file. REQUIRED only for 'custom regions' mode.

-c [colors]: please specify the colors for (sp)KAS-seq data in heatmap plot. Note: you can specify only one colors (e.g. Reds) or specify a color list (e.g. Reds,Greens,Blues,Purples), which needs to be consistent with the number of KAS-seq bigWig files. REQUIRED. Note: please refer to http://matplotlib.org/users/colormaps.html to get the list of valid colors names.

-l [labels.txt]: please input the text file containing the labels of (sp)KAS-seq data to generate metagene profile. Default: basename of (sp)KAS-seq bigWig files.
Example:
WT_rep1
WT.rep2
KO.rep1
KO.rep2                        ---labels.txt

-k [KAS-seq.txt]: please input the text file containing normalized (sp)KAS-seq bigWig files, which can be generated with 'KAS-pipe2 normalize' and 'KAS-pipe2 bedGraphToBigWig' shell scripts. The order and number of (sp)KAS-seq bigWig files should be the consistent with the labels file when specified. REQUIRED.
Example:
KAS-seq_WT_rep1.nor.bigWig
KAS-seq_WT_rep2.nor.bigWig
KAS-seq_KO_rep1.nor.bigWig
KAS-seq_KO_rep2.nor.bigWig     ---KAS-seq.txt

-h/--help: print this help and exit.
Note: The 'KAS-pipe2 heatmap' shell script is applied to generate heatmap plots for (sp)KAS-seq data on genomic features( genebody, TSS or TES) or provided custom regions. 'KAS-pipe2 heatmap' shell script mainly invoke deeptools 'computeMatrix' and 'plotHeatmap', please refer to https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html for more information.
```

### Differential analysis sub-commands:
---------------------------------------------------------------------------

#### KASexpre
#### Usage: KAS-pipe2 KASexpre [ -h/--help ] [ -t threads ] [ -o prefix ] [ -s assembly id ] [ -r regions ] [ -p peaks ] [ -l labels ] [ -k KAS-seq ] 
```
Example: nohup KAS-pipe2 KASexpre -o KAS-seq_expression -t 10 -s mm10 -r TSS -l labels.txt -k KAS-seq_data.txt &

-t [threads]: please specify the number of threads used for calculating KAS expression. Default: 1.

-o [prefix]: please input the prefix (basename) of 'KAS-pipe2 KASexpre' output files. Default: basename of KAS-seq data.

-s [assembly id]: please specify the genome assembly id of KAS-seq data. -s [assembly id]. e.g. Human: hg18, hg19, hg38; Mouse: mm9, mm10, mm39; C.elegans: ce10, ce11; D.melanogaster: dm3, dm6; Rat: rn6, rn7; Zebra fish: danRer10, danRer11. REQUIRED.

-r [regions]: please specify the types of genomic regions. e.g. promoter, genebody, gene or peak. REQUIRED.

-p [peaks]: please specify the custom regions file used for KAS index calculation. OPTIONAL.

-l [labels]: please input the text file containing the labels of (sp)KAS-seq data that show in output files. REQUIRED.
Example:
WT.rep1
WT.rep2
KO.rep1
KO.rep2                       ---labels.txt

-k [KAS-seq]: please input the text file containing the indexed bam files of (sp)KAS-seq data that used to calculate KAS expression. REQUIRED.
Example:
KAS-seq_WT.rep1.bam
KAS-seq_WT.rep2.bam
KAS-seq_KO.rep1.bam
KAS-seq_KO.rep2.bam           ---KAS-seq_data.txt

-h/--help: print this help and exit.
Note: The 'KAS-pipe2 KASexpre' shell script is applied to calculate normalized KAS-seq expression levels on promoter, genebody, genes or custom regions.
```


### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/Ruitulyu/KAS-pipe2/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
