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
# install and load ggplot2 package.
all_packages <- c("ggplot2") 

for (package in all_packages){
  if (!require(package, character.only = TRUE)){
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# install.packages("ggplot2")
# library(ggplot2)

## -----------------------------------------------------------------------------
# use shell input

Distribution <- as.data.frame(read.table(args[1], header = TRUE))
Distribution$Features <- factor(Distribution$Features, levels=c("Promoter", "Exon", "Intron","Terminal3kb","Intergenic"), ordered=TRUE)

blank_theme <- theme_minimal()+
     theme(
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         panel.border = element_blank(),
         panel.grid=element_blank(),
         axis.ticks = element_blank(),
         plot.title=element_text(size=14, face="bold"))

## -----------------------------------------------------------------------------
# Plot the genomic distribution pie chart and save it to png format.
png(file="KAS-seq_peaks_genomic_distribution_pie_chart.png", bg="transparent")
ggplot(Distribution, aes(x="", y=Percentage, fill=Features)) + geom_bar(width = 2, stat = "identity") + coord_polar("y", start=0) + blank_theme + theme(axis.text.x=element_blank()) + scale_fill_manual(values=c("#86b0cc", "#5ba58f", "#bc5090", "#ed7781", "#ffa600"))
dev.off()

# Plot the genomic distribution pie chart and save it to png format.
svg(file="KAS-seq_peaks_genomic_distribution_pie_chart.svg", bg="transparent")
ggplot(Distribution, aes(x="", y=Percentage, fill=Features)) + geom_bar(width = 2, stat = "identity") + coord_polar("y", start=0) + blank_theme + theme(axis.text.x=element_blank()) + scale_fill_manual(values=c("#86b0cc", "#5ba58f", "#bc5090", "#ed7781", "#ffa600"))
dev.off()

# http://www.sthda.com/english/wiki/wiki.php?title=ggplot2-pie-chart-quick-start-guide-r-software-and-data-visualization

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

