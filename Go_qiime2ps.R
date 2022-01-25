#!/usr/bin/env Rscript
# 20200909
# it works after Go_qiime2.sh
# Go_qiime2ps.R -p [project name] -b [tab location] -x [tax location] -t [tree location] -o [tree location]
# output control
# install.packages("optparse")
version <- "1.0"

library("optparse", warn.conflicts = FALSE) #  suppress the message
library("dplyr",warn.conflicts = FALSE)
library("data.table",warn.conflicts = FALSE)
library(dada2)
library(phyloseq)



# command-line options.
option_list <- list(
  
  make_option(c("-p", "--project"), type="character", default=NULL,
              help="add project name (required)." , metavar="name"),
  
  make_option(c("-b", "--tab"), type="character", default=NULL,
              help="Location of input table file path (required)." , metavar="fastq"),
  
  make_option(c("-x", "--tax"), type="character", default=NULL,
              help="Location of tax file path (required)." , metavar="fastq"),
  
  make_option(c("-t", "--tree"), type="character", default=NULL,
              help="Location of tree file path (required)." , metavar="fastq"),
  
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="Location of tree file path (required)." , metavar="fastq"),
  
  make_option(c("--version"), action = "store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)


opt_parser <- OptionParser(
  option_list=option_list,
  usage = "%prog [options] -p [project name] -b [tab location] -x [tax location] -t [tree location] -o [tree location]",

description = paste(
  "\nVersion v1.0",
  "\nThis is a tool for creating phyloseq object from qiime2 output files.",
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
if(is.null(opt$project)) {
  stop("add project name.")
}

if(is.null(opt$tab)) {
    stop("path to forward FASTQs needs to be set.")
}

if(is.null(opt$tax)) {
  stop("tax needs to be set.")
}

if(is.null(opt$tree)) {
  stop("treeneeds to be set.")
}


if(is.null(opt$out)) {
  stop("tree location needs to be set.")
}


#variation
project <- opt$project
tab <- opt$tab
tax <- opt$tax
tree <- opt$tree
location <- opt$out


# table
otu <- read.table(tab, header=T, as.is=T, sep="\t", row.names=1, comment.char="", quote="")


# for tree file
tree <- read_tree(tree, errorIfNULL = T)


# for taxa table
tax <- read.table(tax, sep="\t", row.names = 1, header=T)

tax2 <- tidyr::separate(data.frame(text = tax$Taxon), text, into = c("Kingdom", "Phylum", "Class","Order","Family","Genus","Species"), sep = ";", fill = "right", extra = "drop")
row.names(tax2) = row.names(tax)


taxmat= as.matrix(tax2)

taxmat = gsub(" ", "", taxmat, fixed = TRUE)
taxmat = gsub("NA", "", taxmat, fixed = TRUE)

taxmat_null <- taxmat
#row.names(taxmat_null) <- NULL
taxmat_null[is.na(taxmat_null)] <- ""

### Merge all info for phyloseq object
tax <- as.matrix(taxmat_null)
otu <- as.matrix(otu)

OTU <- otu_table(otu, taxa_are_rows = TRUE)
TAX <- tax_table(tax)
ps <- phyloseq(otu_table(OTU, taxa_are_rows=FALSE), tax_table(TAX), tree);ps

saveRDS(ps, sprintf("%s/ps.%s.%s.rds", location, project, format(Sys.Date(), "%y%m%d")))

