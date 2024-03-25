Data and scripts in this directory reproduce the plots and results of my thesis, "The use of 3D chaos game representation for DNA sequence analysis with applications for building phylogenetic trees"

The purpose of this analysis was to:

1) demonstrate that a mapping of DNA sequences to a 3D chaos game representation and capturing information encoded in this representation can produce phylogenetic trees that are competitive with the results using multiple sequence alignment
2) demonstrate techniques inspired by shape comparison methods that does not require any stretching, padding, or truncating of the DNA sequences in order to be compared


This directory contains:

- example.Rmd = the script to produce the phylogenetic tree results and various explanatory plots in the paper
- simulations.Rmd = the script to produce the simulation results comparing substitution and deletion mutation effects on overall distance between the a sequence and its progenitor sequence
- data/ = directory containing the DNA sequences used to generate the plots and results
- R/ = directory containing functions implementing the 3D CGR