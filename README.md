# Ratio of Homoplasic Offspring (RoHO)

![schematics](https://github.com/DamienFr/RoHO/blob/master/schematics.png)

This repository contains scripts and command line used to compute the RoHO index of a set of high quality filtered homoplasies. It has been used in the publication  

No evidence for increased transmissibility from recurrent mutations in SARS-CoV-2  
Lucy van Dorp*, Damien Richard*, Cedric CS. Tan, Liam P. Shaw, Mislav Acman, François Balloux   
\*contributed equally   
..link..

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
- output several figures, tables and statistics 

This is a two steps pipeline :

## first step in bash / perl:
00.Main_homoplasy.sh

Production of a matrix like

		EPI_ISL_43074	EPI_ISL_43074	EPI_ISL_43074	EPI_ISL_43074
	1912	"ref"	"ref"	"not_ref"	"ref"
	11083	"ref"	"not_ref"	"ref"	"undef"
	23043	"ref"	"ref"	"not_ref"	"ref"

I've included the matrix in each dataset folder (sometimes zipped because too heavy for github) so you can directly run step 2 if you want to.

## second step in R:
00.homoplasy_ratios.R

I go through the phylogeny and for each node annotated as giving rise to an homoplasy i simply count the number of offsprings having "ref", "not_ref" or "undef" alleles based on the input matrix

## Possible caveats:

See publication

