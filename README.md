Various Rscripts to deal with data-wrangling during metagenomic analysis.

* parsePileup.R : parse coverage depth files generated with BBMAP pileup.sh, after mapping individual read sets onto a combined contig set with minimap2 for usage with Metabat2
* jgidepthToAbundance.R : create abundance files and lists from Metabat-like coverage files
* assign_taxonomy_16S.R : assign taxonomy using DADA2 functions and silvadb to a set of (full) 16S rRNA sequences
