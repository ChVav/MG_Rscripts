#!/usr/bin/env Rscript

# Author: Charlotte Vavourakis
# This script creates from depth coverage files used for Metabat(2) separate abundance files and lists for Maxbin2
# Pattern assumed for naming depth coverage files: depth_<contigset>.txt
# For each contig set a folder will be created, containing one abundance file per mapped read set

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  print("No input directory specified, reading all .txt files in working directory")
  dir.in = "./"
} else if (length(args)==1){
  print("Output will be written to working directory.")
  dir.in=args[1]
}

#----------------------

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))

#dir.in <- "./"
depthList <- list.files(dir.in, pattern=".txt")
contigsetList <- gsub("depth_","",depthList)
contigsetList <- gsub(".txt","",contigsetList)

for (i in 1:length(contigsetList)){
  dir.create(paste0(dir.in,contigsetList[i]))
  dat <- read.table(paste0(file.path(dir.in),depthList[i]), header=TRUE)
  dat <- dat %>% select(!ends_with(".var"))
  abund_list <- c()
  for (j in 1:(length(colnames(dat))-3)){
    readset <- colnames(dat[3+j])
    abund_file <- dat %>% select(all_of(c("contigName",readset)))
    abund_filename <- paste0("abund_",gsub("_sorted.bam","",readset))
    abund_list <- c(abund_list,abund_filename)
    output.file <- file(paste0(dir.in,contigsetList[i],"/",abund_filename),"wb")
    write.table(abund_file,
                row.names=FALSE,
                col.names=FALSE,
                file = output.file,
                quote=FALSE,
                append=TRUE,
                sep="\t")
    close(output.file)
  }
  writeLines(abund_list, paste0(dir.in,contigsetList[i],"/abund_list"))
}
