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
Sys.setenv(TAR = "/bin/tar")
Sys.getenv("TAR")
packages <- c("devtools")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
install.packages(setdiff(packages, rownames(installed.packages())))}

# install.packages("devtools")
library(devtools)

packages <- c("kassambara/easyGgplot2")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
devtools::install_github(setdiff(packages, rownames(installed.packages())))}

# devtools::install_github("kassambara/easyGgplot2")
library(easyGgplot2)

## -----------------------------------------------------------------------------
# use shell input

FRiP <- as.data.frame(read.table(args[1], header = TRUE))

## -----------------------------------------------------------------------------
# Plot the fraction of reads in peaks (FRiP) scores and save it to png format.
png(file="KAS-seq_fraction_of_reads_in_peaks.png", bg="transparent")
ggplot2.barplot(data=FRiP, xName='labels', yName="Percentage",
      groupName='Types', groupColors=c('#999999','#E69F00'),
      #background and line colors
      backgroundColor="white", color="black",
      xtitle="spKAS-seq samples", ytitle="Percentage (%)",
      mainTitle="Fraction of reads in peaks (FRiP)",
      removePanelGrid=TRUE,removePanelBorder=TRUE,
      axisLine=c(0.5, "solid", "black")
      )
dev.off()

# Plot the fraction of reads in peaks (FRiP) scores and save it to svg format.
svg(file="KAS-seq_fraction_of_reads_in_peaks.svg", bg="transparent")
ggplot2.barplot(data=FRiP, xName='labels', yName="Percentage",
      groupName='Types', groupColors=c('#999999','#E69F00'),
      #background and line colors
      backgroundColor="white", color="black",
      xtitle="spKAS-seq samples", ytitle="Percentage (%)",
      mainTitle="Fraction of reads in peaks (FRiP)",
      removePanelGrid=TRUE,removePanelBorder=TRUE,
      axisLine=c(0.5, "solid", "black")
      )
dev.off()

# http://www.sthda.com/english/wiki/ggplot2-barplot-easy-bar-graphs-in-r-software-using-ggplot2
## ----sessionInfo--------------------------------------------------------------
sessionInfo()