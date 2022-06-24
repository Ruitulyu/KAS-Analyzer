# KAS-pipe2
## An integrated and flexible toolkit for exploring (sp)KAS-seq data ##
<img src="https://github.com/Ruitulyu/KAS-pipe2/blob/main/image/KAS-pipe2.png"  height="180" align="right" />

----------------------------------------
KAS-pipe2 is a collection of command line tools specifically developped for exploring KAS-seq or strand-specific (sp)KAS-seq data, including the basic processing tools for quality control, reference genome index, raw reads mapping, and heatmaps and summary plots. KAS-pipe2 also includes many novel command line tools and completely new frame design compared to the first version of KAS-pipe. e.g. time-courese(TC) KAS-seq differential analysis, R-loop identification (only for spKAS-seq), single-stranded enhancers identification, motif calling,  pausing and termination index calculation, termination length and so on. 

KAS-seq is a kethoxal-assisted single-stranded DNA sequencing (KAS-seq) approach, based on the fast and specific reaction between N3-kethoxal and guanines in ssDNA. KAS-seq allows rapid (within 5 min), sensitive and genome-wide capture and mapping of ssDNA produced by transcriptionally active RNA polymerases or other processes in situ using as few as 1,000 cells. KAS-seq can also enable definition of a group of enhancers that are single-stranded. Overall, KAS-seq facilitates fast and accurate analysis of transcription dynamics and enhancer activities in both low-input and high-throughput manner. Strand-specific KAS-seq (spKAS-seq) was further developped to identify in vivo genome-wide R-loops with strand specificity using as few as 50,000 cells by mapping the exposed single-stranded DNA (ssDNA) caused by RNA-DNA duplex R-loop formation.

### KAS-seq:
<img src="https://github.com/Ruitulyu/KAS-pipe2/blob/main/image/regular_KAS-seq_procedure.png">


### Strand-specific KAS-seq (spKAS-seq):
<img src="https://github.com/Ruitulyu/KAS-pipe2/blob/main/image/spKAS-seq_procedure.jpg" width="470" height="320">

----------------------------------------

## Citation:

1. Tong, Wu, Ruitu Lyu, Qiancheng You, and Chuan He. [Kethoxal-assisted single-stranded DNA sequencing captures global transcription dynamics and enhancer activity in situ.](https://www.nature.com/articles/s41592-020-0797-9) Nature methods 17, no. 5 (2020): 515-523.

2. Ruitu, Lyu, Tong Wu, Allen C. Zhu, Diana C. West-Szymanski, Xiaocheng Weng, Mengjie Chen, and Chuan He. [KAS-seq: genome-wide sequencing of single-stranded DNA by N3-kethoxal–assisted labeling.](https://www.nature.com/articles/s41596-021-00647-6) Nature protocols (2022): 1-19.

----------------------------------------

### Documentation:

Our [documentation](https://ruitulyu.github.io/KAS-pipe2/) contains more details on the [individual command line tool and usages](https://ruitulyu.github.io/KAS-pipe2/)

>Please see also the [FAQ](https://github.com/Ruitulyu/KAS-pipe2/wiki), which we update regularly. 

>If you have any questions about the KAS-pipe2, please post your questions to [github discussions](https://github.com/Ruitulyu/KAS-pipe2/discussions).

>For more specific **troubleshooting, feedback, and suggestions**, please post to [github issue](https://github.com/Ruitulyu/KAS-pipe2/issues) or send emails to the authors: lvruitu@gmail.com.

-----------------------------------------

### Installation

**Install by cloning KAS-pipe2 git repository on github:**

You can install KAS-pipe2 using command line (linux) by cloning git repository on github:
```
git clone https://github.com/Ruitulyu/KAS-pipe2.git
cd KAS-pipe2
bash ./setup.sh
source ~/.bashrc
	
# If anaconda or miniconda was not installed on your system. #OPTIONAL.
KAS-pipe2 install -conda
	
# Install conda 'KAS-pipe2' environment, which may last 1-2 hours. 
KAS-pipe2 install -KAS-pipe2
	
# Activate conda 'KAS-pipe2' environment.
conda activate KAS-pipe2
```	

------------------------------------	
### Quick start

**Download test data**

Users can download example KAS-seq data in HEK293T cells from Gene Expression Omnibus (GEO):

```
Note: install sra-tools using "conda install -c bioconda sra-tools" if fastq-dump is not available in your system.

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10349532/SRR10349532 ./ &
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10349533/SRR10349533 ./ &
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10349534/SRR10349534 ./ &
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10349535/SRR10349535 ./ &

mv SRR10349532 HEK293T_KAS-Input.rep1.sra
mv SRR10349533 HEK293T_KAS-Input.rep2.sra
mv SRR10349534 HEK293T_KAS-seq.rep1.sra
mv SRR10349535 HEK293T_KAS-seq.rep2.sra
        
fastq-dump HEK293T_KAS-Input.rep1.sra & 
fastq-dump HEK293T_KAS-Input.rep2.sra &
fastq-dump HEK293T_KAS-seq.rep1.sra &
fastq-dump HEK293T_KAS-seq.rep2.sra &	
```	

**Trimming of adapter and poor quality sequence**

```
Note: the "nohup" command is optional, which means "no hang up" and writes the screen output to "nohup.out" file.
nohup KAS-pipe2 trim -a illumina -t 10 -1 HEK293T_KAS-Input.rep1.fastq.gz &
nohup KAS-pipe2 trim -a illumina -t 10 -1 HEK293T_KAS-Input.rep2.fastq.gz &
nohup KAS-pipe2 trim -a illumina -t 10 -1 HEK293T_KAS-seq.rep1.fastq.gz &
nohup KAS-pipe2 trim -a illumina -t 10 -1 HEK293T_KAS-seq.rep2.fastq.gz &
```
	
**Read alignment of KAS-seq and Input control data**

```        
nohup KAS-pipe2 KAS-seq -t 10 -i /absolute path/hg19_Bowtie2Index/hg19 -o HEK293T_KAS-Input.rep1 -s hg19 -1 HEK293T_KAS-Input.rep1_trimmed.fq.gz &
nohup KAS-pipe2 KAS-seq -t 10 -i /absolute path/hg19_Bowtie2Index/hg19 -o HEK293T_KAS-Input.rep2 -s hg19 -1 HEK293T_KAS-Input.rep2_trimmed.fq.gz &
nohup KAS-pipe2 KAS-seq -t 10 -i /absolute path/hg19_Bowtie2Index/hg19 -o HEK293T_KAS-seq.rep1 -s hg19 -1 HEK293T_KAS-seq.rep1_trimmed.fq.gz &
nohup KAS-pipe2 KAS-seq -t 10 -i /absolute path/hg19_Bowtie2Index/hg19 -o HEK293T_KAS-seq.rep2 -s hg19 -1 HEK293T_KAS-seq.rep2_trimmed.fq.gz &
```
	
Here, the user can generate a read alignment summary report: 

```
# Run from "Summary" directory.
KAS-pipe2 statistics -o HEK293T_KAS-seq_statistics -s /absolute path/Summary/ &
```
 
Example summary report: 

|Samples         | Raw reads      | Clean reads   | Mapped reads | Deduplicated reads | Mapping ratios | Duplication ratios |
|     :---:      |     :---:      |     :---:     |    :---:     |        :---:       |      :---:     |        :---:       |
|  KAS-seq.rep1  | 39,249,037     | 39,194,455    | 38,608,066   | 31,874,927         | 98.50%         |  17.44%            |
|  KAS-seq.rep2  | 37,235,447     | 37,195,072    | 36,623,472   | 30,659,014         | 98.46%         |  16.29%            |
|   Input.rep1   | 45,182,939     | 45,162,826    | 44,306,192   | 41,261,883         | 98.10%         |  6.87%             |
|   Input.rep2   | 39,911,067     | 39,886,186    | 39,065,103   | 36,158,950         | 97.94%         |  7.44%             |


**KAS-seq peaks calling**	
         
Call merged KAS-seq peaks with two KAS-seq replicates:	 

```
# Run from "Bed_files" directory.	 
nohup KAS-pipe2 peakscalling -t HEK293T_KAS-seq.rep1.ext150.bed,HEK293T_KAS-seq.rep2.ext150.bed -c HEK293T_KAS-Input.rep1.ext150.bed,HEK293T_KAS-Input.rep2.ext150.bed -o HEK293T_KAS-seq -g hg19 &
```
	
Call KAS-seq peaks with KAS-seq data individually:

```
# Run from "Bed_files" directory. 
nohup KAS-pipe2 peakscalling -t HEK293T_KAS-seq.rep1.ext150.bed -c HEK293T_KAS-Input.rep1.ext150.bed -o HEK293T_KAS-seq.rep1 -g hg19 &
nohup KAS-pipe2 peakscalling -t HEK293T_KAS-seq.rep2.ext150.bed -c HEK293T_KAS-Input.rep2.ext150.bed -o HEK293T_KAS-seq.rep2 -g hg19 &
```
	
**Quality control**

Fingerprint plot of KAS-seq data:

```
# Run from "Bam_files" directory.
nohup KAS-pipe2 fingerprint -t 10 -s hg19 -o HEK293T_KAS-seq_fingerprint -l labels.txt -k KAS-seq_data.txt &

Example "labels.txt" and "KAS-seq_data.txt" :

KAS-seq.rep1                            HEK293T_KAS-seq.rep1_rmdup.bam             
KAS-seq.rep2                            HEK293T_KAS-seq.rep2_rmdup.bam
Input.rep1                              HEK293T_Input.rep1_rmdup.bam
Input.rep2       ---labels.txt          HEK293T_Input.rep2_rmdup.bam        ---KAS-seq_data.txt
```

Example fingerprint plot:

<img src="https://github.com/Ruitulyu/KAS-pipe2/blob/main/image/KAS-seq_fingerprint_plot.png"  height="240" align="middle" > 


Fraction of reads in peaks (FRiP) scores:

```
nohup KAS-pipe2 FRiP -o HEK293T_KAS-seq_FRiP -p peaks_files.txt -l labels.txt -k KAS-seq.txt &

Example "peaks_files.txt", "labels.txt", and "KAS-seq.txt" :

HEK293T_KAS-seq.rep1_peaks.broadPeak                                 KAS-seq.rep1                            HEK293T_KAS-seq.rep1.ext150.bed
HEK293T_KAS-seq.rep2_peaks.broadPeak     ---peaks_files.txt          KAS-seq.rep2     ---labels.txt          HEK293T_KAS-seq.rep2.ext150.bed     ---KAS-seq.txt
```
	
Example fraction of reads in peaks (FRiP) scores:

<img src="https://github.com/Ruitulyu/KAS-pipe2/blob/main/image/KAS-seq_FRiP.png"  height="240" align="middle" > 


Calculate library complexity metrics for KAS-seq data, including PCR Bottlenecking Coefficient and Non-Redundant Fraction (NRF):

```
# Run from "Bam_files" directory.        
nohup KAS-pipe2 complexity -o HEK293T_KAS-seq_complexity -l labels.txt -k KAS-seq.txt &

Example "labels.txt" and "KAS-seq.txt" :

KAS-seq.rep1                              HEK293T_KAS-seq.rep1_sorted.bam
KAS-seq.rep2      ---labels.txt           HEK293T_KAS-seq.rep2_sorted.bam      ---KAS-seq.txt
```
Example library complexity metrics of KAS-seq data in HEK293T cells:
|Samples         | PBC     | Bottlenecking level | NRF     | Complexity |
|     :---:      | :---:   |     :---:           | :---:   | :---:      |
|  KAS-seq.rep1  | 0.88    | None                | 0.84    | Ideal      |
|  KAS-seq.rep2  | 0.89    | None                | 0.85    | Ideal      |
|   Input.rep1   | 0.96    | None                | 0.95    | Ideal      |
|   Input.rep2   | 0.95    | None                | 0.95    | Ideal      |

Calculate the correlation coefficient and pvalue, generate scatterplot for replicates of KAS-seq data:

``` 
# Run from "Bam_files" directory.
nohup KAS-pipe2 correlation -m pearson -t 10 -s hg19 -r bin -p heatmap -o KAS-seq -l labels.txt -k KAS-seq.txt &

Example "labels.txt" and "KAS-seq.txt" :

KAS-seq.rep1                             HEK293T_KAS-seq.rep1.bam
KAS-seq.rep2     ---labels.txt           HEK293T_KAS-seq.rep2.bam     ---KAS-seq.txt
```

Example scatterplot between two replicates of KAS-seq data:

<img src="https://github.com/Ruitulyu/KAS-pipe2/blob/main/image/KAS-seq_pearson_scatterplot.png"  height="260" align="middle" > 


Generate KAS-seq read density files that can be viewed in the UCSC genome browser:

```
# Run from "BedGraph_files" directory.
nohup KAS-pipe2 UCSC -k KAS-seq.txt -n UCSC_track_names.txt &

Example "KAS-seq.txt" and "UCSC_track_names.txt" :

HEK293T_KAS-seq.rep1.ext150.bg                             KAS-seq.rep1
HEK293T_KAS-seq.rep2.ext150.bg     ---KAS-seq.txt          KAS-seq.rep2     ---UCSC_track_names.txt
```

Example snapshot of KAS-seq data custom tracks from UCSC Genome Browser:

<img src="https://github.com/Ruitulyu/KAS-pipe2/blob/main/image/KAS-seq_snapshot.png"  height="240" align="middle" > 


**Normalization of KAS-seq data**

Using the number of uniquely mapped reads:

```
# Run from "BedGraph_files" directory.
nohup KAS-pipe2 normalize -m ratios -k HEK293T_KAS-seq_data.txt -r ratios.txt -b -s hg19 &

Example "HEK293T_KAS-seq_data.txt" and "ratios.txt" :

HEK293T_KAS-seq.rep1.ext150.bg                                           1.2
HEK293T_KAS-seq.rep1.ext150.bg      ---HEK293T_KAS-seq_data.txt          1.4       ---ratios.txt
```

Using Reads Per Kilobase per Million mapped reads (RPKM):

```
# Run from "Bam_files" directory.
nohup KAS-pipe2 normalize -m RPKM -k HEK293T_KAS-seq_data.txt -b -s hg19 &

Example "HEK293T_KAS-seq_data.txt":
Note: the bam files need to be deduplicated.

HEK293T_KAS-seq.rep1.rmdup.bam
HEK293T_KAS-seq.rep2.rmdup.bam      ---HEK293T_KAS-seq_data.txt
```
		
**Generate summary plots of KAS-seq signals**

Generate metagene profile of KAS-seq read density on gene coding regions (TSS, genebody or TES):

```
# Run from "BedGraph_files" directory.
nohup KAS-pipe2 profile -t 10 -s hg19 -o HEK293T_KAS-seq_profile -r genebody -c red,blue,green,purple -l labels.txt -k KAS-seq.normalized.bigWig.txt &
```

Generate heatmap of KAS-seq read density on gene coding regions (TSS, genebody or TES):

```
# Run from "BedGraph_files" directory.
nohup KAS-pipe2 heatmap -t 10 -s hg19 -o HEK293T_KAS-seq_heatmap -r genebody -q -c Reds -l labels.txt -k KAS-seq.normalized.bigWig.txt &
```

Example heatmap and metagene profile of KAS-seq signals on gene coding regions:

 <img src="https://github.com/Ruitulyu/KAS-pipe2/blob/main/image/KAS-seq_metagene_heatmap.png"  height="550" align="middle" > <img src="https://github.com/Ruitulyu/KAS-pipe2/blob/main/image/KAS-seq_metagene_profile.png"  height="240" align="middle" >           


**Calculate transcription-related metrics**

Calculate pausing index (PI) of RNA Pol II:

```
# Run from "Bam_files" directory.
nohup KAS-pipe2 index -o HEK293T_pausing_index -t 10 -s hg19 -i pausing -l labels.txt -k KAS-seq.txt &
```
      
Calculate elongation index (EI) of RNA Pol II:

```
# Run from "BedGraph_files" directory.
nohup KAS-pipe2 KASindex -o HEK293T_elongation_index -t 10 -s hg19 -r gene -l labels.txt -k KAS-seq.txt &
```
	
Calculate termination index (TI) of RNA Pol II:

```
# Run from "Bam_files" directory.
nohup KAS-pipe2 index -o HEK293T_termination_index -t 10 -s hg19 -i termination -l labels.txt -k KAS-seq.txt &
```
	
Calculate termination length (TL) of RNA Pol II:

```
# Run from "Bam_files" directory.
nohup KAS-pipe2 termilength -o HEK293T_termination_length -t 10 -g hg19 -p peaks.txt -l labels.txt -k KAS-seq.txt &
```
	
**Identify single-stranded transcribing enhancers (SST enhancers)**

```
# Run from "BedGraph_files" directory.
nohup KAS-pipe2 SST_enhancer -o HEK293T_SST_enhancers -t 10 -s hg19 -e H3K27ac_enhancers.bed -p HEK293T_KAS-seq_peaks.broadPeak -k HEK293T_KAS-seq.rep1.ext150.nor.bigWig,HEK293T_KAS-seq.rep1.ext150.nor.bigWig &
```

**For more details on the above tools, and other features implemented in KAS-pipe2, such as "identification of genome-wide R-loops" (only for spKAS-seq) and "differential RNA Pols activity analysis" (two-conditions or time-course KAS-seq data), please check the [tutorial](https://ruitulyu.github.io/KAS-pipe2/).**
------------------------------------	

This tool suite is developed by the [Dr. Ruitu Lyu](https://scholar.google.com/citations?user=7nt2ezgAAAAJ&hl=en) at [Prof. Chuan He's lab](https://he-group.uchicago.edu/) of [the University of Chicago](https://www.uchicago.edu/).

[Documentation](https://ruitulyu.github.io/KAS-pipe2/) | [FAQ](https://github.com/Ruitulyu/KAS-pipe2/wiki)
