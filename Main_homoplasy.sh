

folder=dataset_17.1 # 24k iqtree NextStrain masking

high_qual_homoplasies=./${folder}/homoplasies_to_keep.txt
total_consistencyIndexReport=(./${folder}/consistencyIndexReport_*.txt)
annotated_tree=(./${folder}/annotatedNewickTree_*.tree)
alignement_file=(./${folder}/*.aln)
reference_fasta_genome=refseq_sars.fasta

echo $high_qual_homoplasies
echo $total_consistencyIndexReport
echo $annotated_tree
echo $alignement_file
echo $reference_fasta_genome

############################################################################################################
###############  SCRIPT 1 To use only if you want to use a subset of the homoplasies only ##################
###### This script produce a tree file with spurious homoplasies removed from the annotations ##############
############################################################################################################
env high_qual_homoplasies=$high_qual_homoplasies total_consistencyIndexReport=$total_consistencyIndexReport perl -F',' -ane 'BEGIN{
open (GA, "<", $ENV{high_qual_homoplasies}) or die "Cannot open input file $!";
while ($line = <GA>) {
$line =~ s/\r?\n//;
$h_to_keep{$line} = "keep";
}
close GA;

open (FIL, ">", $ENV{total_consistencyIndexReport} . "_filtered")  or die "Cannot open output file $!";

open (GB, "<", $ENV{total_consistencyIndexReport}) or die "Cannot open input file $!";
while ($line = <GB>) {
$a++;
if($a > 1){
@fields = split /\t/,$line;
if(!exists($h_to_keep{$fields[0]})){ $h_to_delete{$fields[0]} = "ok" }else{ print FIL $line }
}else{ print FIL $line }
};
close GB;
close FIL;

};
foreach $key (keys %h_to_delete){
$_ =~ s/(\)|-)$key(:|-)/$1$2/g; s/--/-/g ;  s/\)-:/\):/g ; s/\)-/\)/g ;
} print ' $annotated_tree > ${annotated_tree}_filtered

############################################################################################################
###########  If you don't want to filter out any homoplasies, just do these two cp #########################
############################################################################################################
cp ${total_consistencyIndexReport} ${total_consistencyIndexReport}_filtered
cp ${annotated_tree} ${annotated_tree}_filtered

############################################################################################################
###########  SCRIPT 2 produces a matrix with individuals in columns, homoplasic sites in lines #############
########### and the matrix is filled with the allele types : "ref" not_ref or undef ########################
############################################################################################################

# output example :
	EPI_ISL_43074	EPI_ISL_43074	EPI_ISL_43074	EPI_ISL_43074
1912	"ref"	"ref"	"not_ref"	"ref"
11083	"ref"	"not_ref"	"ref"	"undef"
23043	"ref"	"ref"	"not_ref"	"ref"

###############################################################################################################
########### This script produces a matrix giving the allele of each homoplasy site for each isolate. ##########
############################ Alleles can be "ref" or "not_ref" or "undef" #####################################
###############################################################################################################
############# Detail : Alignement file is read, genomes are stored, ref genome is opened and also stored ######
#################### then i compare isolate and ref allele for each selected homoplasies ######################
#################### "Ns" originating from ambiguous sites in isolates are noted "undef" ######################
###############################################################################################################
env total_consistencyIndexReport=$total_consistencyIndexReport alignement_file=$alignement_file reference_fasta_genome=$reference_fasta_genome perl -e '
open (GC, "<", $ENV{alignement_file}) or die "Cannot open input file $!";
while ($line = <GC>) {
$line =~ s/\r?\n//;
if($line =~ /^>/){ $line =~ s/>//; $a = $line}
else{ s/\r?\n//g ; $h{$a} .= $line };
}
close GC;

foreach $key (sort { $h{$b} <=> $h{$a} or $a cmp $b } keys %h){
print "\t" . $key;
} 
print "\n";

open (GA, "<", $ENV{reference_fasta_genome}) or die "Cannot open input file $!";
while ($line = <GA>) {
$line =~ s/\r?\n//;
if($line !~ /^>/){ $ref .= $line };
}
close GA;

open (GB, "<", $ENV{total_consistencyIndexReport} . "_filtered") or die "Cannot open input file $!";
while ($line = <GB>) {
$line_count++;
if($line_count > 1){
@fields = split /\t/,$line;

$allele_ref = substr $ref, ($fields[0] -1),1;

print $fields[0];

foreach $key (sort { $h{$b} <=> $h{$a} or $a cmp $b } keys %h){
$temp = substr $h{$key}, ($fields[0] -1),1; ;
if(lc($temp) eq "n"){ print "\tundef" ; $nnn++ }else{
if(lc($temp) ne lc($allele_ref)){ print "\tnot_ref"}else{ print "\tref" }}
} 
print "\n";
# print $fields[0] . "\t" . $allele_ref . "\n";
}

}; print STDERR "number of cumulative n at homo positions : $nnn\n";
close GB; ' > ./${folder}/input_matrix


######################################################################################################
############## This R script computes the actual Ratio of Homoplasic Offspring (RoHO) ################
######################################################################################################
######################################################################################################
#### Detail: For each node of the phylogeny, this script counts the number of offsprings #############
########## having "ref", "not_ref" or "undef" alleles based on the input matrix. #####################
########## To be considered, an internal node must meet the following criterion : ####################
############# * No children nodes themselves annotated as carrying the homoplasy. ####################
############# * Having at least two descendant tips of each allele ###################################
############# * An homoplasic position is only considered for RoHo score computation if ##############
################ at least n=3 or n=5 nodes satisfy the two first criteria. ###########################
######################################################################################################

Rscript --vanilla Homoplasy_ratios.R ${folder}

