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
packages <- c("devtools") 

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
install.packages(setdiff(packages, rownames(installed.packages())))}

# install.packages("devtools")
library(devtools)
Sys.setenv(TAR = "/bin/tar")

devtools::install_github("kassambara/easyGgplot2")
library(easyGgplot2)

## -----------------------------------------------------------------------------
# use shell input

saturation <- as.data.frame(read.table(args[1], header = TRUE))

## -----------------------------------------------------------------------------
# Plot the KAS-seq saturation analysis and save it to png format.
png(file="KAS-seq_saturation_plot.png", bg="transparent")
ggplot2.lineplot(data=saturation, xName="Percentage", yName='Percent', 
                 size=2,addPoint=FALSE, # pointSize=4, pointColor = '#999999',
		 backgroundColor="white", color='#999999',
                 xtitle="Uniquely mapped reads (%)", ytitle="Percentage of peaks coverage (%)",
                 mainTitle="KAS-seq saturation analysis",            
                 removePanelGrid=TRUE, removePanelBorder=TRUE,
                 axisLine=c(0.5, "solid", "black"))
dev.off()

# Plot the KAS-seq saturation analysis and save it to svg format.
svg(file="KAS-seq_saturation_plot.svg", bg="transparent")
ggplot2.lineplot(data=saturation, xName="Percentage", yName='Percent',
                 size=2, addPoint=FALSE, # pointSize=4, pointColor = '#999999',
                 backgroundColor="white", color='#999999',
                 xtitle="Uniquely mapped reads (%)", ytitle="Percentage of peaks coverage (%)",
                 mainTitle="KAS-seq saturation analysis",      
                 removePanelGrid=TRUE, removePanelBorder=TRUE,
                 axisLine=c(0.5, "solid", "black"))
dev.off()

# http://www.sthda.com/english/wiki/ggplot2-barplot-easy-bar-graphs-in-r-software-using-ggplot2

## ----sessionInfo--------------------------------------------------------------
sessionInfo()
