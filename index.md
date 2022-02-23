## KAS-pipe2: a user-friendly toolkit for exploring KAS-seq data

### Author: Ruitu Lyu, Chemistry department, the University of Chicago 
---------------------------------------------------------------------------
<script type="text/javascript" src="//rf.revolvermaps.com/0/0/6.js?i=5f2zpl0mkwd&amp;m=7&amp;c=e63100&amp;cr1=ffffff&amp;f=arial&amp;l=0&amp;bv=90&amp;lx=-420&amp;ly=420&amp;hi=20&amp;he=7&amp;hc=a8ddff&amp;rs=80" async="async"></script>
![Image](https://raw.githubusercontent.com/Ruitulyu/KAS-pipe2/main/image/KAS-pipe2.jpg)

### Introduction

KAS-pipe2 is a collection of command line tools specifically developped for exploring KAS-seq or strand-specific (sp)KAS-seq data, including the basic processing tools for quality control, reference genome index, raw reads mapping, and heatmaps and summary plots. KAS-pipe2 also includes many novel features and completely new frame design compared to KAS-pipe. e.g. time-courese(TC) KAS-seq differential analysis, R-loop identification (only for spKAS-seq), ss-enhancers, motif, index calculation, termination length and so on.

KAS-pipe2 is still on active development and the source code is hosted on [GitHub](https://github.com/Ruitulyu/KAS-pipe2).

### Installation

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
	
### Overview

KAS-pipe2 is a collection of command line tools for KAS-seq or strand specific KAS-seq(spKAS-seq) data analysis.
```
Version:   v2.0
About:     KAS-pipe2 is mainly developed by Ruitu Lyu, postdoc fellow in Prof. Chuan He's group at the University of Chicago.
Docs:      https://ruitulyu.github.io/KAS-pipe2/
Code:      https://github.com/Ruitulyu/KAS-pipe2
Mail:      https://github.com/Ruitulyu/KAS-pipe2/discussions
```
#### Usage:
KAS-pipe2 sub-command [options]

#### The KAS-pipe2 sub-commands include:

##### Configure
```	
   download        Downlaod the index of reference genome for aligners (bowtie, bowtie2, bwa, star).
   build           Build the index of reference genome for aligners (bowtie, bowtie2, bwa, star).
   install         Install and check the 'KAS-pipe2' conda environment; check the installation of tools.
   uninstall       Uninstall the 'KAS-pipe2' conda environment.
   activate        Activate the 'KAS-pipe2' conda environment.
   deactivate      Deactivate the 'KAS-pipe2' conda environment.
```
	
##### Fastqc
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
##### Mapping
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

##### Plot
```
   profile         Generate metagene profile for KAS-seq data (normalized bigWig files are needed).
   heatmap         Generate heatmap for KAS-seq data (normalized bigWig files are needed).
```

##### Differential KAS-seq analysis
```
   KASexpre        Calculate normalized KAS-seq expression levels on promoter, genebody, genes or custom regions. 
   diff            Perform differential KAS-seq analysis on promoter, genebody, gene, bin or custom regions.        
   TC              Perform 'case-only' or 'case-control' differential time course(TC) analysis for (sp)KAS-seq data.
   PCA             Perform and plot PCA analysis for (sp)KAS-seq data.
```
##### R-loops
```
   R-loop          Identify R-loops regions with spKAS-seq data.
```   
  
##### Single-stranded enhancers identification
```
   ss_enhancer     Identify the single stranded (ss) enhancers.
   motif           Identify enriched TF binding motifs on ss_enhancers.
```   
  
##### Termination length
```
   termilength     Calculate the transcription termination length of protein coding genes.
```   
 
##### KAS-seq index
```
   index           Calculate the pausing or termination index.
```   

##### General help
```
   --help          Print this help menu.
   --version       Print the version of KAS-pipe2 you are using.
   --contact       Feature requests, bugs, mailing lists, etc.
```   

### Tutorial

#### Configure sub commands:

##### <details><summary>Download</summary>
<p>
##### Usage and option summary:

###### Usage: KAS-pipe2 download [ -l ] [ -h ] [ -a aligner ] [ -g assembly id ] [ -d directory to save index of aligner ]
```   
Example: KAS-pipe2 download -a bowtie2 -g hg19 -d /Software/reference_genome/ 

-l list all of the available aligner index for reference genomes.

-a [aligner]: aligner name you want to use, e.g. bowtie, bowtie2 or bwa. REQUIRED.

-g [assembly id]: reference genome assembly id, e.g. Human: hg18, hg19, hg38; Mouse: mm9,mm10,mm39; Fruitfly: dm3, dm6; Rat: rn6, rn7; C.elegans: ce10, ce11; Zebra fish: danRer10, danRer11. REQUIRED.

-d [directory to save index of aligner]: directory to save downloaded reference genome index. REQUIRED.

-h\-help: print this help and exit.
```   
</p>
</details>

You can use the [editor on GitHub](https://github.com/Ruitulyu/KAS-pipe2/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.


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
