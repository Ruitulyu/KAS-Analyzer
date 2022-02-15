### Ref links
# http://tuxette.nathalievilla.org/?p=1696
# https://swcarpentry.github.io/r-novice-inflammation/05-cmdline/

## -----------------------------------------------------------------------------
# get the input passed from the shell script.
args <- commandArgs(TRUE)
  
# test if there is at least one argument: if not, return an error.
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input:  ", args[1]))
}

## -----------------------------------------------------------------------------
#install and load package.
# install.packages("factoextra")
library(factoextra)
# install.packages("FactoMineR")
library(FactoMineR)

## -----------------------------------------------------------------------------
# use shell input
KAS <- read.table(args[1])
KAS<-t(KAS)
KAS.res.pca <- PCA(KAS, graph = FALSE)


# Visualize eigenvalues/variances and save it as png plot.
png(file="KAS-seq_percentage_of_variances.png", bg="transparent")
fviz_screeplot(KAS.res.pca, addlabels = TRUE)
dev.off()

# Visualize eigenvalues/variances and save it as svg plot.
svg(file="KAS-seq_percentage_of_variances.svg", bg="transparent")
fviz_screeplot(KAS.res.pca, addlabels = TRUE)
dev.off()


# Graph of individuals
# 1. Use repel = TRUE to avoid overplotting
# 2. Control automatically the color of individuals using the cos2
    # cos2 = the quality of the individuals on the factor map
    # Use points only
# 3. Use gradient color
# plot the individuals in the PCA plot and save it as png plot.
png(file="KAS-seq_PCA_plot.png", bg="transparent")
fviz_pca_ind(KAS.res.pca, col.ind = "cos2", pointsize = 3,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
             )
dev.off()

# plot the individuals in the PCA plot and save it as svg plot.
svg(file="KAS-seq_PCA_plot.svg", bg="transparent")
fviz_pca_ind(KAS.res.pca, col.ind = "cos2", pointsize = 3,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
             )
dev.off()

# http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-data-analyses-and-elegant-visualization

## ----sessionInfo--------------------------------------------------------------
sessionInfo()
