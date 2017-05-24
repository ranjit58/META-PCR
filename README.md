# unique-snp-id-finder
Find unique snp regions for discriminating samples using a sliding window-based approach for VCF data.

vcf_tools.py
- contains general tools for VCF files

region_search.py
- command line parameters:
	- in a terminal type python region_search.py
- contains various functions for manipulating the VCF data after it has been read in
- contains a main function for running the pipeline on a vcf file
- uses the command line interface tool to read in arguments when the program is run from the terminal

This pipeline has been designed with RAM limitations in mind, including: 
- reading in files line by line (to serve as the "moving window") (files are commonly GB in size so reading into memory is not feasible)
- performing computationally expensive operations as infrequently as possible
- writing output that is easily accessible and parsable for direct analysis

Output File
- the output file is structured such that is has some header and file information
- the analysis information is sorted from largest to smallest based on the magnitude of the filtered data
- the analysis information is structured in the following way:
  | filtered score | filtered total samples | unfiltered score | genomic bp bounds of the window |
  - the filtered data uses a simple % threshold (if a specified % of the data for that window is missing information, then discard it)
  - the unfiltered data does not threshold (meaning that '.' for all samples would be valid as unique)
  - the bounds are the actual basepair values of the genome for which that data returns the filtered score and unfiltered score
	- the filtered total samples are the total number of samples with a percentage of basepairs below the threshold that have not been checked for uniqueness
	- the filtered score is the total number of samples with a percentage of basepairs below the threshold that have been checked for uniqueness
	- the unfiltered score is the total number of samples that are unique (without thresholding)


