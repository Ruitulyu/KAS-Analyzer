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
#install and load easyggplot2 package.
all_packages <- c("devtools")

for (package in all_packages){
  if (!require(package, character.only = TRUE)){
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# install.packages("devtools")
# library(devtools)
Sys.setenv(TAR = "/bin/tar")
Sys.getenv("TAR")

all_packages <- c("kassambara/easyGgplot2")

for (package in all_packages){
  if (!require("easyGgplot2", character.only = TRUE)){
    devtools::install_github(package, dependencies = TRUE)
    library("easyGgplot2", character.only = TRUE)
  }
}

# devtools::install_github("kassambara/easyGgplot2")
# library(easyGgplot2)

## -----------------------------------------------------------------------------
# use shell input
fragmentsize <- as.data.frame(read.table(args[1], header = TRUE))

## -----------------------------------------------------------------------------
# Plot the fragment size of KAS-seq data and save it to png format.
png(file="KAS-seq_fragment_size_density_plot.png", bg="transparent")
ggplot2.density(data=fragmentsize, xName='Length', groupName='Samples',
    legendPosition="right",
    backgroundColor="white",
    xtitle="Length (bp)", ytitle="Density",
    mainTitle="Fragment size density curve",
    removePanelGrid=TRUE,removePanelBorder=TRUE,
    axisLine=c(0.5, "solid", "black"), xlim=c(0,500),
    fillGroupDensity=TRUE, alpha=0.5, addMeanLine=TRUE )
dev.off()

# Plot the fragment size of KAS-seq data and save it to svg format.
svg(file="KAS-seq_fragment_size_density_plot.svg", bg="transparent")
ggplot2.density(data=fragmentsize, xName='Length', groupName='Samples',
    legendPosition="right",
    backgroundColor="white",
    xtitle="Length (bp)", ytitle="Density",
    mainTitle="Fragment size density curve",
    removePanelGrid=TRUE,removePanelBorder=TRUE,
    axisLine=c(0.5, "solid", "black"), xlim=c(0,500),
    fillGroupDensity=TRUE, alpha=0.5, addMeanLine=TRUE )
dev.off()

# http://www.sthda.com/english/wiki/wiki.php?title=ggplot2-density-easy-density-plot-using-ggplot2-and-r-statistical-software

## ----sessionInfo--------------------------------------------------------------
sessionInfo()
