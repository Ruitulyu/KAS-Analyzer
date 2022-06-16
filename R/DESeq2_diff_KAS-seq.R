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
  print(paste0("Arg input:  ","option1:", args[1]," option2:", args[2]," option3:", args[3]))
}

## -----------------------------------------------------------------------------
# use shell input
KAS.matrix <- round(as.matrix(read.table(args[1], header = T)))
coldata <- read.table(args[2], header = T)

## -----------------------------------------------------------------------------
# install and load DESeq2 package
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager",repos = "http://cran.us.r-project.org")
all_packages <- c("BiocManager") 

for (package in all_packages){
  if (!require(package, character.only = TRUE)){
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}


all_packages <- c("DESeq2")

for (package in all_packages){
  if (!require(package, character.only = TRUE)){
    BiocManager::install(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}


# BiocManager::install("DESeq2")
# library(DESeq2)

## -----------------------------------------------------------------------------
#generate DESeqDataSet
coldata$condition <- factor(coldata$condition)
dds <- DESeqDataSetFromMatrix(countData = KAS.matrix,
                              colData = coldata,
                              design = ~ condition)

#Differential analysis is based on normalized KAS-seq data, so we need to set the size factor of every sample to 1.
sizeFactors(dds)<- rep(1, times=ncol(KAS.matrix))

#Differential KAS-seq analysis
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
foldchange <- as.double(print(paste0("",args[3])))
DE.KAS=resOrdered[abs(resOrdered$log2FoldChange)>log2(foldchange) & resOrdered$padj < 0.05,]

## -----------------------------------------------------------------------------
#output the results of differential KAS-seq analysis
write.csv(resOrdered,"KAS-seq_DESeq2_output.csv")
write.csv(DE.KAS,"DE.KAS-seq_DESeq2_Fold_padj0.05_output.csv")

## ----sessionInfo--------------------------------------------------------------
sessionInfo()


