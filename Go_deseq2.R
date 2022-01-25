#!/usr/bin/env Rscript
# Wrote by Heekuk Park
# 210310
# Go_deseq2.R -p [project name] -i [psIN] -t [trimLeft] -d [file location] # reverse strand trimLeft 5 fixed
# output control
# install.packages("optparse")
# 20200929
# qiime2 와 같은 과정으로 제작
# foreign install problem
# install.packages("https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-76.tar.gz")


version <- "1.1"

bioconductors <- c("dada2","phyloseq","DESeq2","dplyr")

for (bioconductor in bioconductors){
  if(!bioconductor %in% installed.packages()){
    library(BiocManager)
    BiocManager::install(bioconductor)
  }else{library(bioconductor, character.only = TRUE)}
}

packages <- c("optparse","dplyr", "data.table") #"venneuler",

for (package in packages){
    if(!package %in% installed.packages()){
        install.packages(package)
    }else{library(package, character.only = TRUE)}
}



# command-line options.
option_list <- list(
  
  make_option(c("-p", "--project"), type="character", default=NULL,
              help="add project name (required)." , metavar="project"),
  
  make_option(c("-i", "--psIN"), type="character", default=NULL,
              help="Location of input psIN file path (required)." , metavar="psIN"),
  
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Location of input metadata file path (required)." , metavar="metadata"),
  
  make_option(c("-t", "--type"), type="character", default=NULL,
              help="Location of data type  (required)." , metavar="type"),
  
  make_option(c("-d", "--data_type"), type="character", default=NULL,
              help="Location of datatype file path (required)." , metavar="datatype"),
  
  make_option(c("-o", "--order"), type="character", default=NULL,
              help="Location of datatype order path (required)." , metavar="order"),
  
  make_option(c("-a", "--adjust"), type="character", default=NULL,
              help="Location of datatype order path (required)." , metavar="adjust"),
  
  make_option(c("-s", "--description"), type="character", default=NULL,
              help="Location of datatype order path (required)." , metavar="description"),
  
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="Location of datatype order path (required)." , metavar="name"),
  
  make_option(c("-c", "--alpha"), type="character", default=NULL,
              help="Location of datatype order path (required)." , metavar="alpha"),
  
  make_option(c("--version"), action = "store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)


opt_parser <- OptionParser(
  option_list=option_list,
  usage = "%prog [options] -p [project name] -i [psIN] -m [metadata] -t [type] -d [datatype] -o [order] -a [ajdust] -s [description] -n [name] -c [alpha]",

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
if(is.null(opt$psIN)) {
  stop("add psIN .")
}

if(is.null(opt$data_type)) {
    stop("add data_type ")
}


#if(is.null(opt$order)) { stop("add order")}

if(is.null(opt$adjust)) {
  stop("add adjst")
}


if(is.null(opt$description)) {
  stop("add des")
}

if(is.null(opt$name)) {
  stop("add nane")
}


if(is.null(opt$alpha)) {
  stop("d alpha")
}

#variation
project <- opt$project
psIN <- opt$psIN
metadata <- opt$metadata
type <- opt$type
data_type <- opt$data_type
# order <- opt$order
adjust <- opt$adjust
des <- opt$des
name <- opt$name
alpha <- opt$alpha



#==========#
#=  START =#
#==========#

# out dir
out <- file.path(sprintf("%s_%s",project, format(Sys.Date(), "%y%m%d"))) 
if(!file_test("-d", out)) dir.create(out)
out_path <- file.path(sprintf("%s_%s/table",project, format(Sys.Date(), "%y%m%d"))) 
if(!file_test("-d", out_path)) dir.create(out_path)
out_deseq2 <- file.path(sprintf("%s_%s/table/deseq2",project, format(Sys.Date(), "%y%m%d"))) 
if(!file_test("-d", out_deseq2)) dir.create(out_deseq2)

metadata <- read.csv(sprintf("%s",metadata),header=T,as.is=T,row.names=1,check.names=F)


# map 정리
psIN <- readRDS(psIN)

mapping <- data.frame(sample_data(psIN))
sel <- intersect(rownames(metadata), colnames(mapping)); head(sel, "3")
metadata.sel <- metadata[sel,, drop=F];head(metadata.sel)
mapping.sel <- mapping[rownames(mapping), sel, drop=F];head(mapping.sel)

dim(mapping.sel)

print("pass1")

# 이전 버전 for strafied
# if(type == "function"){
#   # remove colume sum 0 and psIN 재구성(20201027)
#   a <- data.frame(otu_table(psIN))
#   b <- a[, -which(numcolwise(sum)(a) < 1)]
#   OTU.sta <- otu_table(b, taxa_are_rows = TRUE);head(OTU.sta)
#   colnames(OTU.sta) <- gsub("X", "", colnames(OTU.sta))
#   otu_table(psIN) <-  OTU.sta
# }else if(type == "taxanomy"){
#   psIN <- psIN
# }

# 최근 버전 for unstrafied (20210112 확인)
if(type == "function"){
  # remove colume sum 0 and psIN 재구성(20201027)
  a <- data.frame(otu_table(psIN))*10000
  a.ceiling <- ceiling(a[-c(99),])
  b <- a.ceiling[, -which(numcolwise(sum)(a.ceiling) < 1)]
  if (length(b) == 0){
    OTU.sta <- otu_table(a, taxa_are_rows = TRUE);head(OTU.sta)
    colnames(OTU.sta) <- gsub("X", "", colnames(OTU.sta))
    otu_table(psIN) <-  OTU.sta
  }else if(length(b) > 1){
    OTU.sta <- otu_table(b, taxa_are_rows = TRUE);head(OTU.sta)
    colnames(OTU.sta) <- gsub("X", "", colnames(OTU.sta))
    otu_table(psIN) <-  OTU.sta
  }
}else if(type == "taxanomy"){
  psIN <- psIN
}else if(type == "bacmet"){
  psIN <- psIN
}



# start
res <- {}
for (mvar in rownames(subset(metadata.sel, Go_deseq2 =="yes"))) {
  if (length(unique(mapping.sel[, mvar])) == 1) {
    next
  }
  
  #na remove
  mapping.sel <- data.frame(sample_data(psIN))
  mapping.sel[mapping.sel==""] <- "NA"
  mapping.sel.na <- mapping.sel[!is.na(mapping.sel[,mvar]), ]
  na.count <- length(mapping.sel.na)
  psIN.na <- prune_samples(rownames(mapping.sel[!is.na(mapping.sel[,mvar]), ]), psIN)
  mapping.sel.na.rem <- data.frame(sample_data(psIN.na ))
  if (length(unique(mapping.sel.na.rem[,mvar])) == 1 )
    next
  
  if (length(des) == 1) {
    print(sprintf("##-- %s-%s (total without NA: %s/%s) --##",
                  des,mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
    
  } else{
    print(sprintf("##-- %s (total without NA: %s/%s) --##",
                  mvar, dim(mapping.sel.na.rem)[1], dim(mapping.sel)[1]))
  }
  
  if (length(mapping.sel.na.rem[,mvar]) < 4){
    next
    print(sprintf("%s is removed because length(%s) less than 4", mvar, length(mapping.sel.na.rem[,mvar])))
  }
  
  # integer control
  if (class(mapping.sel.na.rem[,mvar]) == "character"){
    mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
    sample_data(psIN.na) <- mapping.sel.na.rem
  }
  if (class(mapping.sel.na.rem[,mvar]) == "integer" | class(mapping.sel.na.rem[,mvar]) == "numeric"){
    mapping.sel.na.rem[,mvar] <- factor(mapping.sel.na.rem[,mvar])
    sample_data(psIN.na) <- mapping.sel.na.rem
  }
  
  #-- DESeq2 for phyloseq --#
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  
  if (length(adjust) >= 1) {
    form <-as.formula(sprintf("~ %s + %s", mvar, paste(setdiff(adjust, "SampleType"), collapse="+")))
    print(form)
    dds = phyloseq_to_deseq2(psIN.na, form)
    
  }    else {
    dds = phyloseq_to_deseq2(psIN.na, as.formula(sprintf("~ %s", mvar)))
    print(sprintf("~ %s", mvar))
  }
  
  
  geoMeans = apply(counts(dds), 1, gm_mean)
  dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  dds = estimateDispersions(dds)
  vst = getVarianceStabilizedData(dds)
  dds = DESeq(dds, fitType="local")
  resultsNames(dds)
  
  
  # go back to the default order
  
  # mapping.sel[,mvar] <- factor(mapping.sel[,mvar], levels = orders)
  
  mapping.sel[,mvar] <- factor(mapping.sel[,mvar])
  cbn <- combn(x = levels(mapping.sel[,mvar]), m = 2)
  
  my_comparisons <- {}
  for(i in 1:ncol(cbn)){
    x <- cbn[,i]
    my_comparisons[[i]] <- x
  };my_comparisons
  
  # calsculation
  
  for(i in 1:length(my_comparisons)){
    print(my_comparisons[i])
    
    combination <- unlist(my_comparisons[i]);combination
    
    basline <-combination[1]
    smvar <- combination[2]
    print("pass2")
    tmp <- results(dds, contrast = c(mvar, smvar, basline))
    tmp$taxa <- unlist(lapply(rownames(tmp), function(x) {
      tmp <- unlist(strsplit(x, ";"))
      tmp[length(tmp)]
    }))
    
    tmp$dir <- ifelse(tmp$padj < alpha, ifelse(sign(tmp$log2FoldChange)==1, "up", "down"), "NS")
    tmp$mvar <- mvar
    tmp$basline<-basline
    tmp$smvar <- smvar
    if (length(des) == 1) {
      tmp$des <- des
    }
    
    #-- give taxa name --#
    res <- cbind(as(tmp, "data.frame"), as(tax_table(psIN)[rownames(tmp), ], "matrix"))
    print("pass3")
    #res$TaxaName <- paste("p__",res$Phylum,";c__",res$Class,";o__",res$Order,";f__",res$Family,";g__",res$Genus,";s__",res$Species)
    # removing no name and NA
    #for(taxa in c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
    
    
    if(type == "taxonomy"){
      for(taxa in c("Phylum","Class","Order","Family","Genus","Species")){
        res[,taxa] == "NA"
        res[,taxa]<- as.character(res[,taxa])
        res[,taxa][is.na(res[,taxa])] <- "__"
        for(i in 1:length(res[,taxa])){
          if (res[,taxa][i] == "s__" || res[,taxa][i] == "g__" || res[,taxa][i] == "f__" || res[,taxa][i] == "o__" || res[,taxa][i] == "c__"|| res[,taxa][i] == "p__"|| res[,taxa][i] == "__"){
            res[,taxa][i] <- ""
          }
        }
      }
      print("pass4")
      res$TaxaName <- paste(res$Phylum,"",res$Class,"",res$Order,"",res$Family,"",res$Genus,"",res$Species)
      #res$ShortName <- paste(res$Phylum,res$Family," ",res$Genus," ",res$Species)
      
      if (data_type == "dada2" | data_type == "DADA2") {
        res$ShortName <- paste(res$Genus,"",res$Species)
      }
      else if (data_type == "Nephele" | data_type == "nephele") {
        res$ShortName <- paste(res$Genus,"",res$Species)
      }
      else if (data_type == "other" | data_type == "Other") {
        res$ShortName <- paste(res$Species)
      }
      
      
      # use last taxa name
      for(taxa in c("Family", "Order", "Class","Phylum")){
        for(i in 1:length(res[,taxa])){
          if (res$ShortName[i] != "  "){
            next
          }      else if (res$ShortName[i] == "  " & res[,taxa][i] != ""){
            res$ShortName[i] <- paste(res[,taxa][i])
          }
        }
      }
    } else if(type == "function"){
      for(taxa in c("KO", "KO.des","Path","Path.des")){
        res[,taxa] == "NA"
        res[,taxa]<- as.character(res[,taxa])
        res[,taxa][is.na(res[,taxa])] <- "__"
        for(i in 1:length(res[,taxa])){
          if (res[,taxa][i] == "s__" || res[,taxa][i] == "g__" || res[,taxa][i] == "f__" || res[,taxa][i] == "o__" || res[,taxa][i] == "c__"|| res[,taxa][i] == "p__"|| res[,taxa][i] == "__"){
            res[,taxa][i] <- ""
          }
        }
      }
      print("pass4")
      res$KOName <- paste(res$Path,"",res$KO)
      res$ShortName <- paste(res$Path.des,"",res$KO.des)
      
      # use last taxa name
      for(taxa in c("KO", "KO.des","Path","Path.des")){
        for(i in 1:length(res[,taxa])){
          if (res$ShortName[i] != "  "){
            next
          }      else if (res$ShortName[i] == "  " & res[,taxa][i] != ""){
            res$ShortName[i] <- paste(res[,taxa][i])
          }
        }
      }
    }else if(type == "bacmet"){
      for(taxa in c("Gene",	"Organism",	"Compound",	"NCBI_annotation")){
        res[,taxa] == "NA"
        res[,taxa]<- as.character(res[,taxa])
        res[,taxa][is.na(res[,taxa])] <- "__"
        for(i in 1:length(res[,taxa])){
          if (res[,taxa][i] == "s__" || res[,taxa][i] == "g__" || res[,taxa][i] == "f__" || res[,taxa][i] == "o__" || res[,taxa][i] == "c__"|| res[,taxa][i] == "p__"|| res[,taxa][i] == "__"){
            res[,taxa][i] <- ""
          }
        }
      }
      print("pass4")
      res$TaxaName <- paste(res$Compound,"",res$Gene,"",res$Organism)
      res$ShortName <- paste(res$Compound,"",res$Gene,"",res$Organism)
    }
    
    #--- give simple name to res---#
    headers <- vector(dim(res)[2], mode="character")
    for (i in 1:dim(res)[1]) {
      headers[i] <- paste("ASV", i, sep="_")
    }
    res$taxa <- headers
    print("pass5")
    #-- create table --#
    res <- as.data.frame(res)
    res$padj <- p.adjust(res$pvalue, method="fdr")
    res$dir <- ifelse(res$padj < alpha, ifelse(sign(res$log2FoldChange)==1, "up", "down"), "NS")
    
    if (length(des) == 1) {
      if (length(name) == 1) {
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/(%s.vs.%s).%s.%s.%s.%s.Deseq2.csv",out_deseq2,basline,smvar, mvar, des, name, project,sep="/"))
      } else {
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/(%s.vs.%s).%s.%s.%s.Deseq2.csv",out_deseq2,basline,smvar,mvar,des,project,sep="/"))
      }
    } else {
      if (length(name) == 1) {
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/(%s.vs.%s).%s.%s.%s.Deseq2.csv",out_deseq2, basline,smvar,mvar,name,project,sep="/"))
      } else{
        write.csv(res, quote = FALSE,col.names = NA,file=sprintf("%s/(%s.vs.%s).%s.%s.Deseq2.csv",out_deseq2,basline, smvar,mvar, project, sep="/"))
      }
    }
  }
}
