
#!/usr/bin/env Rscript

# Author: Charlotte Vavourakis
# This script assigns silvadb taxonomy to a set of 16S rRNA gene sequences
# download first databases to working directory: 
# wget --content-disposition "https://zenodo.org/records/3986799/files/silva_nr99_v138_train_set.fa.gz?download=1"
# wget --content-disposition "https://zenodo.org/records/3986799/files/silva_species_assignment_v138.fa.gz?download=1"
# to do:amend adding also gtdb taxonomy?

# arg, lib load ##----

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("No input fasta specified")
} else {
  fastafile = args[1]
  cat("Assuming databases are in working directory.\nOutput will be written to working directory as well.")
}

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)})

# silva taxonomy ##----

#fastafile <- "./7-summaries-phyloflash/derep_assembled_16S.fasta"
fasta <- readDNAStringSet(fastafile)
seqs <- c()
for (i in 1:length(fasta)){
  seqs[i] <- as.character(toString(fasta[i]))
}
headers <- names(fasta)

## assign species does not work for non-bases ##----

# identify problematic sequences and separate them out
check_non_base <- function(string) {
  contains_non_base <- grepl("[^ACGTU]", string)
  return(contains_non_base)
}
check <- unlist(lapply(seqs, check_non_base))
sum(check) 
names(check) <- headers
true_names <- names(check[check])

filtered_seqs <- fasta[!names(fasta) %in% true_names]
seqs_not_assigned <- fasta[names(fasta) %in% true_names]
writeXStringSet(seqs_not_assigned, "seqs_no_species_assigned.fasta")

## Down to genus ##----
seqs1 <- c()
for (i in 1:length(filtered_seqs)){
  seqs1[i] <- as.character(toString(filtered_seqs[i]))
}
headers1 <- names(filtered_seqs)

seqs2 <- c()
for (i in 1:length(seqs_not_assigned)){
  seqs2[i] <- as.character(toString(seqs_not_assigned[i]))
}
headers2 <- names(seqs_not_assigned)

set.seed(100) # Initialize random number generator for reproducibility
taxa1 <- assignTaxonomy(seqs1, "./silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
taxa2 <- assignTaxonomy(seqs2, "./silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
rownames(taxa2) <- headers2
saveRDS(taxa2, file = paste0(fastafile, "_assignedtaxsilva_onlygenus.Rds"))

## Add species ##---
taxa1 <- addSpecies(taxa1, "./silva_species_assignment_v138.fa.gz")
rownames(taxa1) <- headers1

saveRDS(taxa1, file = paste0(fastafile, "_assignedtaxsilva_complete.Rds"))

# GTDB taxonomy ##----

