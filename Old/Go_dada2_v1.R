#!/usr/bin/env Rscript
# 191025
# Go_dada2.R -p [project name] -f [file location] -t [trimLeft] -d [file location]
# output control
# install.packages("optparse")
version <- "1.0"

library("optparse", warn.conflicts = FALSE) #  suppress the message
library("dplyr",warn.conflicts = FALSE)
library("data.table",warn.conflicts = FALSE)
library(dada2)
library(phyloseq)
multiplot <- function(..., plotlist=NULL, file, cols=1, rows=1) {
  require(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  
  i = 1
  while (i < numPlots) {
    numToPlot <- min(numPlots-i+1, cols*rows)
    layout <- matrix(seq(i, i+cols*rows-1), ncol = cols, nrow = rows, byrow=T)
    if (numToPlot==1) {
      print(plots[[i]])
    } else {

      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

      for (j in i:(i+numToPlot-1)) {

        matchidx <- as.data.frame(which(layout == j, arr.ind = TRUE))
        print(plots[[j]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
    i <- i+numToPlot
  }
}




# command-line options.
option_list <- list(
  
  make_option(c("-p", "--project"), type="character", default=NULL,
              help="add project name (required)." , metavar="name"),
  
  make_option(c("-f", "--fastq"), type="character", default=NULL,
              help="Location of input fastq file path (required)." , metavar="fastq"),
  
  make_option(c("-d", "--db"), type="character", default=NULL,
              help="Location of input fastq file path (required)." , metavar="fastq"),
  
  make_option(c("-t", "--trimLeft"), type="integer", default=NULL,
              help="Location of taxa file path (required)." , metavar="path"),
  
  make_option(c("--version"), action = "store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)


opt_parser <- OptionParser(
  option_list=option_list,
  usage = "%prog [options] -p [project name] -f [file location] -t [trimLeft] -d [file location]",

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

if(is.null(opt$fastq)) {
    stop("path to forward FASTQs needs to be set.")
}


if(is.null(opt$trimLeft)) {
  stop("path to input taxa file needs to be set.")
}


#variation
project <- opt$project
path <- opt$fastq
trimleft.value <- opt$trimLeft
db <- opt$db


# output


dada2 <- file.path(sprintf("%s_dada2",project))
if(!file_test("-d", dada2)) dir.create(dada2)
out <- file.path(sprintf("%s_dada2/%s",project,"1_out")) 
if(!file_test("-d", out)) dir.create(out)
rds <- file.path(sprintf("%s_dada2/%s",project,"2_rds")) 
if(!file_test("-d", rds)) dir.create(rds)


# run 


#--------------------------#
#      qulity control      #
#--------------------------#
#-- input --#

fnFs <- list.files(path, pattern="_L001_R1_001.fastq.gz");head(fnFs)
fnRs <- list.files(path, pattern="_L001_R2_001.fastq.gz");head(fnRs)
sample.names <- sapply(strsplit(fnFs, "_L001_R1_001.fastq.gz"), `[`, 1)

#-- plot quality profiles before filter--#
pdf(sprintf("%s/%s.qualityProfiles.%s.pdf", dada2,project, format(Sys.Date(), "%y%m%d")))
for (sn in sample.names) {
  plist <- list()
  fn <- list.files(path, pattern=sprintf("%s_L001_R1_001.fastq.gz", sn), full.names =T)
  plist[[length(plist)+1]] <- plotQualityProfile(fn)
  fn <- list.files(path, pattern=sprintf("%s_L001_R2_001.fastq.gz", sn), full.names =T)
  plist[[length(plist)+1]] <- plotQualityProfile(fn)
  print(multiplot(plotlist=plist, cols=1, rows=2))
}  
dev.off()



#-------------------------#
#   sequences filter      #
#-------------------------#
filt_path <- file.path(sprintf("%s/3_DADA2_filtered_%s", dada2,format(Sys.Date(), "%y%m%d")))
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R2_filt.fastq.gz"))

fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq.gz", full.names = TRUE))


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,250), trimLeft=c(trimleft.value, trimleft.value), 
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, 
                     multithread=T, verbose=TRUE, matchIDs=TRUE)
head(out)

#- plot quality profiles after filter -#

pdf(sprintf("%s/%s.qualityProfiles.filt.%s.pdf", dada2,project, format(Sys.Date(), "%y%m%d")))
for (sn in sample.names) {
  plist <- list()
  fn <- sprintf("%s/3_DADA2_filtered_%s/%s_R1_filt.fastq.gz", dada2,format(Sys.Date(), "%y%m%d"), sn)
  plist[[length(plist)+1]] <- plotQualityProfile(fn)
  fn <- sprintf("%s/3_DADA2_filtered_%s/%s_R2_filt.fastq.gz", dada2,format(Sys.Date(), "%y%m%d"), sn)
  plist[[length(plist)+1]] <- plotQualityProfile(fn)
  print(multiplot(plotlist=plist, cols=1, rows=2))
}
dev.off()


#-------------------------------------#
# Dereplication for unique sequences  #   
#_____________________________________#

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE) 


names(derepFs) <- sample.names
names(derepRs) <- sample.names


# option 1 fan 방법 인듯 하다. 좀더 오래된 방법. Apaia 분석 떄문에 찾아 보았다. batch 이전 sorce코드로 사용할떄 썼던 방법 (20210316 작성)
set.seed(100)
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out

dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out

# option 2 보통 사용하던거 같은데...언제 부터 사용된것인지 확인은 안되지만, batch 코드 초기 버전인듯 하다. (20210316 작성)
#set.seed(100)
#errF <- learnErrors(filtFs, multithread=TRUE)
#set.seed(100)
#errR <- learnErrors(filtRs, multithread=TRUE)






dadaFs <- dada(derepFs, err=errF, multithread=TRUE) #pool = TRUE, ‘pooled’
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


pdf(sprintf("%s/%s.splotErrors.errF1.%s.pdf",dada2,project ,format(Sys.Date(), "%m%d%y")))
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf(sprintf("%s/%s.splotErrors.errF2.%s.pdf",dada2,project ,format(Sys.Date(), "%m%d%y")))
plotErrors(errR, nominalQ=TRUE)
dev.off()


# Mergers
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#head(mergers[[1]])


seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE, multithread=TRUE)


#---------------------------------#
# sequences 처리에 따른의 갯수변화#
#---------------------------------#
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

out <- file.path(sprintf("%s_dada2/%s",project,"1_out")) 
if(!file_test("-d", out)) dir.create(out)
write.csv(track, quote = FALSE,col.names = NA, file=sprintf("%s/%s.%s.track.csv",out, project, format(Sys.Date(), "%y%m%d"), sep="/"))

# save seqtab.nochim
rds <- file.path(sprintf("%s_dada2/%s",project,"2_rds")) 
if(!file_test("-d", rds)) dir.create(rds)
saveRDS(seqtab.nochim, sprintf("%s/seqtab.nochim.%s.%s.rds", rds, project, format(Sys.Date(), "%y%m%d")))


#---------------------------------#
#     save seqs.fna for tree      #
#---------------------------------#
seqs <- getSequences(seqtab.nochim)
headers <- paste(">", seqs, sep="")
fasta <- c(rbind(headers, seqs))
write(fasta, file=sprintf("%s/%s.%s.seqs.fna",  out,project, format(Sys.Date(), "%y%m%d"),sep="/"))

print("fna file is saved for tree.")


tax <- assignTaxonomy(seqtab.nochim, db, 
                      taxLevels = c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species"), 
                      minBoot = 80, verbose = TRUE, multithread = TRUE)
rds <- file.path(sprintf("%s_dada2/%s",project,"2_rds")) 
if(!file_test("-d", rds)) dir.create(rds)
saveRDS(tax, sprintf("%s/tax.%s.%s.rds", rds,project, format(Sys.Date(), "%y%m%d")))

#---------------------------------#
#             for ps              #
#---------------------------------#
taxa.print <- tax # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(tax))
rds <- file.path(sprintf("%s_dada2/%s",project,"2_rds")) 
if(!file_test("-d", rds)) dir.create(rds)
saveRDS(ps, sprintf("%s/ps.%s.%s.rds", rds, project, format(Sys.Date(), "%y%m%d")))

#--------------------------------
# Extracting the standard goods from DADA2
#--------------------------------
seqs <- colnames(seqtab.nochim)
headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final OTU seqs:
#fasta <- c(rbind(headers, seqs))
#out <- file.path(sprintf("%s_dada2/%s",project,"1_out")) 
#if(!file_test("-d", out)) dir.create(out)
#write(fasta, file=sprintf("%s/%s.%s.fasta", out, project, format(Sys.Date(), "%m%d%y"),sep="/"))

# count table:
otu <- t(seqtab.nochim)
row.names(otu) <- sub(">", "", headers)
out <- file.path(sprintf("%s_dada2/%s",project,"1_out")) 
if(!file_test("-d", out)) dir.create(out)
write.csv(otu, 
          quote = FALSE,
          #row.names = FALSE, 
          col.names = NA,
          file=sprintf("%s/%s.%s.asv.csv",  out,project, format(Sys.Date(), "%y%m%d"),sep="/"))

# tax table:
taxatab <- tax
row.names(taxatab) <- sub(">", "", headers)
out <- file.path(sprintf("%s_dada2/%s",project,"1_out")) 
if(!file_test("-d", out)) dir.create(out)
write.csv(taxatab, 
          quote = FALSE,
          #row.names = FALSE, 
          col.names = NA,
          file=sprintf("%s/%s.%s.tax.csv", out,project,format(Sys.Date(), "%y%m%d"), sep="/"))


otuTable <- cbind(otu,tax)
out <- file.path(sprintf("%s_dada2/%s",project,"1_out")) 
if(!file_test("-d", out)) dir.create(out)
write.csv(otuTable, 
          quote = FALSE,
          #row.names = FALSE, 
          col.names = NA,
          file=sprintf("%s/%s.%s.asvTable.csv",  out,project,format(Sys.Date(), "%y%m%d"), sep="/"))


