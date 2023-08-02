#!/usr/bin/env Rscript

# Author: Charlotte Vavourakis
# assumes individual reads were mapped with minimap2 to a combined set of contigs from different samples
# assumes BBmap pileup.sh script was then run to extract coverage information from the .bam files (one for each read set)
# This script parses such individual depth files in specified directory into depth files for each separate contig set as input for Metabat(2) 

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("No input directory specified, reading all bam.depth files in working directory")
	dir.in = "."
} else if (length(args)==1){
  print("Output will be written to working directory.")
	dir.in=args[1]
}

#----------------------

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))

#dir.in <- "."
alignedList <- list.files(dir.in, pattern="depth")
headers <- scan(alignedList[1], what = "", nlines = 1, sep="\t")[-1]
dat <- as.data.frame(matrix(nrow=0,ncol=13))
colnames(dat) <- c("contigset", "contigname", headers, "readset")
dat[] <- sapply(dat, as.integer)
dat$contigset <- as.character(dat$contigset)
dat$contigname <- as.character(dat$contigname)
dat$Avg_fold <- as.numeric(dat$Avg_fold)
dat$Ref_GC <- as.numeric(dat$Ref_GC)
dat$Covered_percent <- as.numeric(dat$Covered_percent)
dat$Read_GC <- as.numeric(dat$Read_GC)
dat$Std_Dev <- as.numeric(dat$Std_Dev)
dat$readset <- as.character(dat$readset)

for (i in 1:length(alignedList)){
  dat2 = read.table(alignedList[i]) %>%
    separate(V1,
             c("contigset", "contigname"),
             sep = "C",
             remove = TRUE)
  dat2$readset <- gsub(".bam.depth","",alignedList[i])
  colnames(dat2) = c("contigset", "contigname", headers, "readset")
  dat <- bind_rows(dat,dat2)
}
colnames(dat) <- c("contigset", "contigname", headers, "readset")

contigsetList <- dat %>% pull(contigset) %>% unique()

for (i in 1:length(contigsetList)){
  dat2 = dat %>%
    filter(contigset==contigsetList[i]) %>%
    select(contigname,Length,readset,Avg_fold,Std_Dev) %>%
    pivot_wider(names_from = readset,
                names_glue = "{readset}_{.value}",
                values_from = c(Avg_fold,Std_Dev)) %>%
    mutate(totalAvgDepth = rowSums(select(., ends_with('Avg_fold')))) %>%
    dplyr::rename(contigName=contigname,contigLen=Length)
  colnames(dat2) = gsub("Avg_fold","sorted.bam",colnames(dat2))
  colnames(dat2) = gsub("Std_Dev","sorted.bam-var",colnames(dat2))
  dat2 = dat2[,order(colnames(dat2))] %>%
    relocate(totalAvgDepth) %>%
    relocate(contigLen) %>%
    relocate(contigName)
  output.file <- file(paste0(dir.in,"/depth_",contigsetList[i],".txt"),"wb")
  write.table(dat2,
              row.names=FALSE,
              col.names=TRUE,
              file = output.file,
              quote=FALSE, 
              append=TRUE,
              sep="\t")
  close(output.file)
}
