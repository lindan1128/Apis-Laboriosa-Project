## Quality control and genome assembly
## Gene prediction
## Gene family construction
	# Count the number of all types of gene families
	python count_all_gene_family --gene_family_directory /home/lin/gene_families/ --file_suffix cds.muscle --gene_family_number 10460
## Phylogenetic analysis
### Species tree construction
	# Concatenate all single-copy gene families
	cat *.conserved_gene_family.muscle >> all_conserved_gene_family.fasta
	# Change fasta format into phy format
	perl fasta2philip.pl all_conserved_gene_family.fasta > 	all_conserved_gene_family.phy
	# Run phyml
	phyml -i all_conserved_gene_family.phy -b 1000 -d nt -m GTR -c 4 -a e -o tlr
### Genetic distance estimate 
###  Divergence time estimation
## Positive selection analysis
## Gene duplication/loss analysis