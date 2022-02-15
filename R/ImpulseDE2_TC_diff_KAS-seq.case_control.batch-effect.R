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
KAS.matrix <- round(as.matrix(read.table(args[1], header = T)))
coldata <- read.table(args[2], header = T)

## -----------------------------------------------------------------------------
# install and load ImpulseDE2 package
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager",repos = "http://cran.us.r-project.org")

# BiocManager::install("ImpulseDE2")
library(ImpulseDE2)

## -----------------------------------------------------------------------------
