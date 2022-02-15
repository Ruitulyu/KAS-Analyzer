### Ref links
# http://tuxette.nathalievilla.org/?p=1696
# https://swcarpentry.github.io/r-novice-inflammation/05-cmdline/

## -----------------------------------------------------------------------------
# get the input passed from the shell script.
args <- commandArgs(TRUE)

# test if there is at least one argument: if not, return an error.
if (length(args) == 0) {
  stop("At least two argument must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input:  ", args[1], args[2]))
}

## -----------------------------------------------------------------------------
#install and load ggpubr package.
#install.packages("ggpubr",repos = "http://cran.us.r-project.org")
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")

library("ggpubr")

## -----------------------------------------------------------------------------
# use shell input.
KAS.matrix <- as.data.frame(read.table(args[1], header = T))

## -----------------------------------------------------------------------------
# Plot the correlation scatterplot and save it to png format.
png(file="KAS-seq_corr_scatterplot.png", bg="transparent")
ggscatter(KAS.matrix, x = colnames(KAS.matrix)[1], y = colnames(KAS.matrix)[2],
          color = "lightgray", add = "reg.line", conf.int = TRUE,
          add.params = list(color = "steelblue",fill = "steelblue")) +
          stat_cor(method = args[2])

dev.off()

# Plot the correlation scatterplot and save it to svg format.
svg(file="KAS-seq_corr_scatterplot.svg", bg="transparent")
ggscatter(KAS.matrix, x = colnames(KAS.matrix)[1], y = colnames(KAS.matrix)[2],
          color = "lightgray", add = "reg.line", conf.int = TRUE,
          add.params = list(color = "steelblue",fill = "steelblue")) +
          stat_cor(method = args[2])
dev.off()

#http://www.sthda.com/english/wiki/wiki.php?title=ggplot2-scatterplot-easy-scatter-plot-using-ggplot2-and-r-statistical-software 

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

