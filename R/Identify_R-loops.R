### Ref links
# http://tuxette.nathalievilla.org/?p=1696
# https://swcarpentry.github.io/r-novice-inflammation/05-cmdline/

## -----------------------------------------------------------------------------
# get the input passed from the shell script
args <- commandArgs(TRUE)

## -----------------------------------------------------------------------------
# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least two argument must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input:  ","option1:", args[1]," option2:", args[2]))
}

## -----------------------------------------------------------------------------
# use shell input
Rloops.matrix <- round(as.matrix(read.table(args[1], header = T)))
coldata <- read.table(args[2], header = T)

## -----------------------------------------------------------------------------
# install and load DESeq2 package
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager",repos = "http://cran.us.r-project.org")

# BiocManager::install("DESeq2")
library("DESeq2")

## -----------------------------------------------------------------------------
#generate DESeqDataSet
coldata$condition <- factor(coldata$condition)
dds <- DESeqDataSetFromMatrix(countData = Rloops.matrix,
                              colData = coldata,
                              design = ~ condition)

#specify the reference level
dds$condition <- relevel(dds$condition, ref = "minus")

#Differential analysis is based on normalized KAS-seq data, so we need to set the size factor of every sample to 1.
sizeFactors(dds)<- rep(1, times=ncol(Rloops.matrix))

#Differential KAS-seq analysis
dds <- DESeq(dds)
res <- results(dds, name="condition_plus_vs_minus")
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
Rloops=resOrdered[abs(resOrdered$log2FoldChange)>log2(1.5) & resOrdered$padj <0.05 ,]

## -----------------------------------------------------------------------------
#output the results of differential KAS-seq analysis
write.csv(resOrdered,"spKAS-seq_plus_vs_minus_DESeq2_R-loops_output.csv")
write.csv(Rloops,"spKAS-seq_plus_vs_minus_DESeq2_R-loops_Fold1.5_padj0.05_output.csv")

## ----sessionInfo--------------------------------------------------------------
sessionInfo()
