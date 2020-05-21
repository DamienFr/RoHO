# Ratio of Homoplasic Offspring (RoHO)

This is a two steps pipeline :

## first step in bash / perl:
00.homoplasyv2.sh

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

I always consider that the "not_ref" allele is the "homoplasic one". This should be ok as the ref allele is near the root ...

I do not take into account nested reported homoplasic nodes

