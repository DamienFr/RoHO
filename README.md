# SARS-CoV-2 Ratio of Homoplasic Offspring (RoHO)

![schematics](https://github.com/DamienFr/RoHO/blob/master/schematics.png)

This repository contains scripts and command line used to compute the RoHO index of a set of high quality filtered homoplasies. It has been used in the publication  

**No evidence for increased transmissibility from recurrent mutations in SARS-CoV-2**  
Lucy van Dorp*, Damien Richard*, Cedric CS. Tan, Liam P. Shaw, Mislav Acman, François Balloux   
\*contributed equally   
https://doi.org/10.1101/2020.05.21.108506 

Homoplasies are mutations that emerged repeatedly and independently. They are good candidates for sites under natural selection. In our study, we use a phylogenetic index to assess whether particular homoplasic mutations increase transmissibility of SARS-CoV-2. To do so, we quantify the relative number of descendants in sister clades with and without a specific allele. This github repository contains the code used to do so.

The RoHO index is the ratio of the number of descendents of sister clades with and without a specific mutation over all independent emergences of a homoplasic allele in a phylogeny.

The script 'Main_homoplasy.sh" is to be run in a unix command line step by step to ensure proper functioning.
You are likely to have to make some modifications to it in order to modify paths etc.


** This script is very likely to crash as it was not developped for general use but for the study of one particular dataset **

## Inputs
- A fasta reference genome (not provided, Wuhan-Wu-1 download here :https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id=NC_045512.2)
- A list of high quality homoplasy positions (provided)
- Homoplasyfinder Consistency Index Report (provided, output file of Homoplasyfinder)
- Annotated tree file (not provided due to copyright issues, output file of Homoplasyfinder)
- Alignement file (not provided due to copyright issues, fasta file of aligned sequences from Gisaid)

## Outputs
- Outputs several figures, tables and statistics 

## Dependencies
Perl core  
ape R package  
phangorn R package  
reshape2 R package  
ggplot2 R package  


## This is a two steps pipeline :

### First step in bash / perl:
Main_homoplasy.sh

From a filtered homoplasy list, this script :
* Converts a raw homoplasyfinder annotated tree to a filtered homoplasyfinder annotated tree
* Produces a matrix of alleles for each isolates at variable positions from an alignement fasta file

The matrix has the following form :

		Isolate_1	Isolate_2	Isolate_3	Isolate_4
	1912	"ref"	"ref"	"not_ref"	"ref"
	11083	"ref"	"not_ref"	"ref"	"undef"
	23043	"ref"	"ref"	"not_ref"	"ref"



### Second step in R:
Homoplasy_ratios.R

In the filtered phylogeny (tree file), 
For each filtered node of the phylogeny annotated by HomoplasyFinder as corresponding to an ancestor that acquired a homoplasy, this script counts the number of offsprings having "ref", "not_ref" or "undef" alleles based on the input matrix.
To be considered, an internal node must meet the following criterion :
* No children nodes themselves annotated as carrying the homoplasy.   
* Having at least three descendant tips.    
* An homoplasic position is only considered for RoHo score computation if at least n=5 or n=10 nodes satisfy the two first criteria.

# Acknowledgements

This study has been possible only thanks to the very large number of scientists in originating and submitting labs who have readily made available SARS-CoV-2 assemblies to the research community. They are listed in originating-laboratories.txt and submitting-laboratories.txt files.

A test dataset with no copyright issue will be provided in order to get the scripts running
