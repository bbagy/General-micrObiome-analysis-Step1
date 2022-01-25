#!/usr/bin/env Rscript
# 190525
# working good
# GoBlastCount.R -i [input] -o [output]
# output control
# merge_metaphlan_tables.py *.txt > merged_blast.txt"
# install.packages("optparse")
version <- "1.0"

library("optparse", warn.conflicts = FALSE) #  suppress the message
library("dplyr",warn.conflicts = FALSE)
library("data.table",warn.conflicts = FALSE)

# command-line options.
option_list <- list(
  
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Location of input blast out file (required)." , metavar="path"),
  
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Location of output path (required)." , metavar="path"),
  
  make_option(c("--version"), action = "store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)


opt_parser <- OptionParser(
  option_list=option_list,
  usage = "%prog [options] -i blastfile -o count.txt",

description = paste(
  "\nVersion v1.0",
  "\nThis is a counting tool for the blast format 6 tabular output file to build the count table.",
   sep=" ")
)

opt <- parse_args(opt_parser)

# Print out version if --version flag set.
if (opt$version) {
  cat("Wrapper version:", version, "\n")
  options_tmp <- options(show.error.messages=FALSE) 
  on.exit(options(options_tmp)) 
  stop()
}

# Check if path to input blast file set, if not stop job.
if(is.null(opt$input)) {
    stop("path to input blast file needs to be set.")
}


# Check if path to out file set, if not stop job.
if(is.null(opt$output)) {
  stop("path to output needs to be set.")
}

# run counting blast file
bt <- read.table(opt$input, header=F, as.is=T, sep="\t", comment.char="")
sample.names <- sapply(strsplit(opt$input, "_top.txt"), `[`, 1)
bt <- data.frame(bt$V2)
bc <- bt %>% group_by(bt.V2) %>% tally()
setnames(bc, old=c("bt.V2","n"), new=c("gene", sample.names))
write.table(bc, file = opt$output, quote=F, row.names=F, col.names=T, sep="\t")
print(sample.names)


# system("merge_metaphlan_tables.py *.txt > merged_blast.txt")



