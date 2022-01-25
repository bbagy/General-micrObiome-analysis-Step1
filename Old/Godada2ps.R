#!/usr/bin/env Rscript
# 191025
# Godada2ps.R -p [project name] -s [seqtab.nochim] -x [tax] -tr [tree] -m [map]
# output control
# install.packages("optparse")
version <- "1.0"

library("optparse", warn.conflicts = FALSE) #  suppress the message
library("dplyr",warn.conflicts = FALSE)
library("data.table",warn.conflicts = FALSE)

# command-line options.
option_list <- list(
  
  make_option(c("-p", "--project"), type="character", default=NULL,
              help="add project name (required)." , metavar="name"),
  
  make_option(c("-s", "--seqtab"), type="character", default=NULL,
              help="Location of input seqtab file path (required)." , metavar="path"),
  
  make_option(c("-x", "--taxa"), type="character", default=NULL,
              help="Location of taxa file path (required)." , metavar="path"),
  
  make_option(c("-t", "--tree"), type="character", default=NULL,
              help="Location of input tree file path (required)." , metavar="path"),
  
  make_option(c("-a", "--alpha"), type="character", default=NULL,
              help="Location of output path (required)." , metavar="value"),
  
  make_option(c("--version"), action = "store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)


opt_parser <- OptionParser(
  option_list=option_list,
  usage = "%prog [options] -p [project name] -s [seqtab.nochim] -x [tax] -t [tree] -m [map] -o [output]",

description = paste(
  "\nVersion v1.0",
  "\nThis is a tool for creating phyloseq object from dada2 output files.",
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

if(is.null(opt$seqtab)) {
    stop("path to input seqtab file needs to be set.")
}

if(is.null(opt$taxa)) {
  stop("path to input taxa file needs to be set.")
}

if(is.null(opt$tree)) {
  stop("path to input tree file needs to be set.")
}

if(is.null(opt$map)) {
  stop("path to input map file (.csv) needs to be set.")
}

if(is.null(opt$alpha)) {
  stop("value for filter.")
}

# run 

project <- opt$project
# reads RDS
# ps
seqtab.nochim <- readRDS(opt$seqtab)
tax <- readRDS(opt$tax)
treefile <- read_tree(opt$tree, errorIfNULL = T)
sampledata <- read.csv(opt$map,row.names=1,check.names=FALSE);head(sampledata)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(tax));ps

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(tax), phy_tree(tree));ps

saveRDS(ps, sprintf("2_rds/ps.%s.%s.rds",project, format(Sys.Date(), "%y%m%d")))

phylo_relabun <- transform_sample_counts(ps, function(x) x / sum(x))
phylo_filter = filter_taxa(phylo_relabun, function(x) mean(x) < sprintf("%s", opt$alpha),TRUE) #.00005
rmtaxa = taxa_names(phylo_filter)
alltaxa = taxa_names(phylo_relabun)
myTaxa = alltaxa[!alltaxa %in% rmtaxa]
ps_relabun_filtered <- prune_taxa(myTaxa,phylo_relabun); ps_relabun_filtered
ps_filtered <- prune_taxa(myTaxa,ps2);ps_filtered


saveRDS(ps_filtered, sprintf("2_rds/ps_filtered.%s.%s.rds",project, format(Sys.Date(), "%y%m%d")))
