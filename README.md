# chromium
R package for the analysis and visualization of Hi-C data.
--------------------------------

* We import preprocessed Hi-C data and convert them into the InteractionSet format (from Bioconductor). 
* We provide functions to easily subset the InteractionSet by only the genomic region required, discarding unneccessary 
genomic annotation to reduce object size.
* The Hi-C data can be binned at any resolution, starting from either binned or restriction fragment resolution data.
* We normalize the interactions using iterative proportional fitting, a matrix balancing algorithm.
* The interactions can be visualized in a triangular visualization alongside ChIPSeq tracks, the genomic axis, etc.

For more info, please visit the package vignette.
