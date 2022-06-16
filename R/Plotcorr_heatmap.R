
### Ref links
# http://tuxette.nathalievilla.org/?p=1696
# https://swcarpentry.github.io/r-novice-inflammation/05-cmdline/

## -----------------------------------------------------------------------------
# get the input passed from the shell script.
args <- commandArgs(TRUE)

## -----------------------------------------------------------------------------
# test if there is at least one argument: if not, return an error.
if (length(args) == 0) {
  stop("At least two argument must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input:  ", args[1]))
}

## -----------------------------------------------------------------------------
# use shell input.
KAS.matrix <- as.data.frame(read.table(args[1], header = T))

## -----------------------------------------------------------------------------
# install and load corrplot package.
all_packages <- c("corrplot","RColorBrewer")

for (package in all_packages){
  if (!require(package, character.only = TRUE)){
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# install.packages("corrplot",repos = "http://cran.us.r-project.org")
# install.packages("RColorBrewer",repos = "http://cran.us.r-project.org")
# library("corrplot")
# library("RColorBrewer")

## -----------------------------------------------------------------------------
# please refer to http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram.
# calculate the correlation of KAS-seq matrix.
KAS.matrix_corr <- round(cor(KAS.matrix),2)

# calculate the pvalue of correlation of KAS-seq matrix.
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.matrix <- cor.mtest(KAS.matrix)

write.csv(KAS.matrix_corr,"KAS-seq_correlation.csv")
write.csv(p.matrix,"KAS-seq_correlation_pvalue.csv")

## -----------------------------------------------------------------------------
# Plot the correlation heatmap plot and save it to png format.
# If you want to save to other format, please refer:https://d.cosx.org/d/16930-16930/4
png(file="KAS-seq_corr_heatmap.png", bg="transparent")
# The visualization method : "circle", "color", "number", etc.
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = KAS.matrix_corr, col = col, symm = TRUE)
dev.off()

# Plot the correlation heatmap plot and save it to svg format.
# If you want to save to other format, please refer:https://d.cosx.org/d/16930-16930/4
svg(file="KAS-seq_corr_heatmap.svg", bg="transparent")
# The visualization method : "circle", "color", "number", etc.
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = KAS.matrix_corr, col = col, symm = TRUE)
dev.off()

# Plot the correlation circle plot and save it to png format.
# If you want to save to other format, please refer:https://d.cosx.org/d/16930-16930/4
png(file="KAS-seq_corr_circle.png", bg="transparent")
# The visualization method : "circle", "color", "number", etc.
col<- colorRampPalette(c("blue", "white", "red"))(20)
corrplot(KAS.matrix_corr, type="full", order="hclust", 
         col=brewer.pal(n=8, name="PuOr"), tl.col="black", tl.srt=45)
dev.off()

# Plot the correlation circle plot and save it to svg format.
svg(file="KAS-seq_corr_circle.svg", bg="transparent")
# The visualization method : "circle", "color", "number", etc.
corrplot(KAS.matrix_corr, type="full", order="hclust", 
         col=brewer.pal(n=8, name="PuOr"), tl.col="black", tl.srt=45)
dev.off()

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

