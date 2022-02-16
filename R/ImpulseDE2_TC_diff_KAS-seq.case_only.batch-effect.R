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
KASseq <- as.matrix(read.table(args[1], header = T))
annotation <- as.data.frame(read.table(args[2], header = T))

## -----------------------------------------------------------------------------
# install and load ImpulseDE2 package
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager",repos = "http://cran.us.r-project.org")

# BiocManager::install("ImpulseDE2")
library(ImpulseDE2)

# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

## -----------------------------------------------------------------------------
# perform case-only time course(TC) KAS-seq differential analysis using ImpulseDE2 package with batch effect normalization.
objectImpulseDE2 <- runImpulseDE2(
     matCountData    = KASseq,
     dfAnnotation    = annotation,
     boolCaseCtrl    = FALSE,
     boolIdentifyTransients = TRUE,
     vecConfounders  = c("Batch"),
     scaNProc        = 1 )

## -----------------------------------------------------------------------------
#output the results of time course differential KAS-seq analysis
write.csv(objectImpulseDE2$dfImpulseDE2Results,"KAS-seq_case_only_TC_ImpulseDE2_nor-batch_output.csv")

## -----------------------------------------------------------------------------
# Visualization of global TC KAS-seq patterns via a heatmap

# Visualize global TC KAS-seq patterns and save it as png plot.
png(file="KAS-seq_case_only_TC_ImpulseDE2_nor-batch_heatmap.png", bg="transparent")
plotHeatmap(
  objectImpulseDE2       = objectImpulseDE2,
  strCondition           = "case",
  boolIdentifyTransients = TRUE,
  scaQThres              = 0.01)
dev.off()

# Visualize global TC KAS-seq patterns and save it as svg plot.
png(file="KAS-seq_case_only_TC_ImpulseDE2_nor-batch_heatmap.svg", bg="transparent")
plotHeatmap(
  objectImpulseDE2       = objectImpulseDE2,
  strCondition           = "case",
  boolIdentifyTransients = TRUE,
  scaQThres              = 0.01)
dev.off()

# https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/ImpulseDE2/inst/doc/ImpulseDE2_Tutorial.html#transiently-regulated-genes

## ----sessionInfo--------------------------------------------------------------
sessionInfo()
