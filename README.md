General_pop_gen_functions
=========================

Here are some resources for you in this Git reprository:

1) A collection of functions useful for population genetics analyses
      - Note: Some I have created myself and others I have modified from the library(SOLOMON)
2) A script to provide you with some basic population genetics statistics describing your data
      - Statistics include:
          - Number of samples (N)
          - Number of samples missing genotypes (N.missing)
          - Percent individuals genotyped (P.GT)
          - Number of alleles (A)
          - Allelic richness (Ar)
          - Number of alleles unique to a population (Au)
          - Observed heterozygosity (Ho)
          - Expected heterozygosity (He)
          - Fis
      - These statistics are summarized per locus per population, and averaged across all loci per population
      - Some basic graphics are provided
3) A script that performs a simple PCA on your genotypes
      - Provides a 3D interactive graphic describing your data
4) A script to calculate pairwise Fst between sample sites
      - Provides Fst estimates
      - Provides a heatmap of Fst values between populations
