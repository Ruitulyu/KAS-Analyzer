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
- [Configure](#Configure)
```	
   download	    Downlaod the index of reference genome for aligners (bowtie, bowtie2, bwa, star).
   build	    Build the index of reference genome for aligners (bowtie, bowtie2, bwa, star).
   install	    Install and check the 'KAS-pipe2' conda environment; check the installation of tools.
   uninstall       Uninstall the 'KAS-pipe2' conda environment.
   activate        Activate the 'KAS-pipe2' conda environment.
   deactivate      Deactivate the 'KAS-pipe2' conda environment.
```	
- #### [Fastqc](#Fastqc)
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
- #### [Mapping](#Mapping)
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

- #### [Plot](#Plot)
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

### Configure

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
