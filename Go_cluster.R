#!/usr/bin/env Rscript
# 20220114
# working good
# Go_cluster.R -p [project] -i [input] -o [output] -d [db] -c [cluster] -t [nproc]
# output control
# merge_metaphlan_tables.py *.txt > merged_blast.txt"
# install.packages("optparse")
version <- "1.0"

bioconductors <- c("dada2","phyloseq","DECIPHER")

for (bioconductor in bioconductors){
  if(!bioconductor %in% installed.packages()){
    library(BiocManager)
    BiocManager::install(bioconductor)
  }else{library(bioconductor, character.only = TRUE)}
}

packages <- c("optparse","dplyr", "data.table","tibble") #"venneuler",

for (package in packages){
  if(!package %in% installed.packages()){
    install.packages(package)
  }else{library(package, character.only = TRUE)}
}




install.packages("remotes")
remotes::install_github("mikemc/speedyseq")
library(speedyseq)
# Packages that are required but not loaded:
# library(DECIPHER)
# library(Biostrings)


library("optparse", warn.conflicts = FALSE) #  suppress the message
library("dplyr",warn.conflicts = FALSE)
library("data.table",warn.conflicts = FALSE)

# command-line options.
option_list <- list(
  
  make_option(c("-p", "--project"), type="character", default=NULL,
              help="project name (required)." , metavar="path"),
  
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Location of input ps object file (required)." , metavar="path"),
  
  make_option(c("-d", "--db"), type="character", default=NULL,
              help="Location of database (required)." , metavar="path"),
  
  make_option(c("-c", "--clustering"), type="character", default=NULL,
              help="% of desired clustering (required)." , metavar="path"),
  
  make_option(c("-t", "--nproc"), type="character", default=NULL,
              help="number of processor (required)." , metavar="path"),
  
  make_option(c("--version"), action = "store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)


opt_parser <- OptionParser(
  option_list=option_list,
  usage = "%prog [options] -p [project] -i [input] -d [db] -c [cluster] -t [nproc]",
  
  description = paste(
    "\nVersion v1.0",
    "\nThis is a clustering tool for DADA2 AVSs table.",
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


if(is.null(opt$project)) {
  stop("needs project name.")
}


# Check if path to input blast file set, if not stop job.
if(is.null(opt$input)) {
  stop("path to input blast file needs to be set.")
}

if(is.null(opt$db)) {
  stop("path to input db needs to be set.")
}

# Check if path to out file set, if not stop job.

if(is.null(opt$clustering)) {
  stop("need clustering %")
}

if(is.null(opt$nproc)) {
  stop("need number of nproc")
}

Go_cluster <- function(psIN, project,db, percent){
}

project <- opt$project;project
psIN <- opt$input
psIN <- readRDS(psIN);psIN
db <- opt$db
percent <- as.numeric(opt$clustering);percent
nproc <- as.numeric(opt$nproc);nproc



# out dir
out <- file.path(sprintf("%s_%s_clustered%s",project, format(Sys.Date(), "%y%m%d"), percent)) 
if(!file_test("-d", out)) dir.create(out)


# ----- Input ------#
project.name <-sprintf("%s_%s",project,percent);project.name


ps <- psIN
x <- 1-percent/100

seqtab <- otu_table(psIN)




asv_sequences <- colnames(seqtab);head(asv_sequences)
sample_names <- rownames(seqtab);(sample_names)
dna <- Biostrings::DNAStringSet(asv_sequences)


## Find clusters of ASVs to form the new OTUs
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
clusters <- DECIPHER::IdClusters(
  d, 
  method = "complete",
  cutoff = x , # use `cutoff = 0.03` for a 97% OTU 
  processors = nproc)

## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
cluster <- clusters %>%
  add_column(sequence = asv_sequences)

merged_seqtab <- seqtab %>% 
  t %>%
  rowsum(clusters$cluster) %>%
  t


# rebuilt ASVs table 
clustered <- distinct(cluster, cluster, .keep_all = TRUE)

merged_seqtab.t <- data.frame(t(merged_seqtab))
merged_seqtab.t$seqs <- factor(clustered$sequence[match(as.factor(rownames(merged_seqtab.t)), as.factor(clustered$cluster))])

rownames(merged_seqtab.t) <- merged_seqtab.t$seqs
merged_seqtab.t$seqs <- NULL

seqtab <- as.matrix(t(merged_seqtab.t))

#----- save seqs.fna for tree  -----#
seqs <- getSequences(seqtab)
headers <- paste(">", seqs, sep="")
fasta <- c(rbind(headers, seqs))
write(fasta, file=sprintf("%s/ASVs%s.%s.%s.seqs.fna",  out, percent,project, format(Sys.Date(), "%y%m%d"),sep="/"))


#----- re classification  -----#
tax <- assignTaxonomy(seqtab, db, 
                      taxLevels = c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species"), 
                      minBoot = 80, verbose = TRUE, multithread = TRUE)


#----- merge  -----# 
ps_clustered <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(tax));ps_clustered

#----- setting fix species names  -----# 

tax <- data.frame(tax_table(ps_clustered))
tax$Species.1 <- paste(tax$Genus,tax$Species)

tax$Species <- NULL
colnames(tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_table(ps_clustered) <- as.matrix(tax)




cat("#--  Before clustered  --#\n")
print(psIN)
cat("\n")
cat(sprintf("#--  After clustered by %s   --#\n",percent))
print(ps_clustered)

saveRDS(ps_clustered, sprintf("%s/ps%s_filtered.%s.%s.rds", out, percent, project.name,format(Sys.Date(), "%y%m%d")))

otu <- as.data.frame(t(otu_table(ps_clustered)));dim(otu)
tax <- tax_table(ps_clustered);dim(tax)
otuTable <- cbind(otu,tax)

write.csv(otuTable, quote = FALSE,col.names = NA,row.names = T, file=sprintf("%s/ASVs%s_clustered.%s.%s.csv", out, percent, project.name,format(Sys.Date(), "%y%m%d"), sep="/"))


