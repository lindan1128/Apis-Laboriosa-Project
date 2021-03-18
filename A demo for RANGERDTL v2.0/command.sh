# Run optRoot
./OptRoot.linux -i 1.gene_family.unroot.newick -D 2 -T 1000 -L 1 | grep ";" >> 1.gene_family.opt.newick
# Cat species tree and gene tree
cat speices_tree_rooted.newick 1.gene_family.opt.newick >> 1.gene_family.newick
# Run Ranger-DTL 
mkdir OutputFiles1
for((i=1;i<=50;i++)); do
	./ranger-dtl-U.linux --seed $i -i 1.gene_family.newick -D 2 -T 1000 -L 1 -o OutputFiles1/rangerOutput$i
done
# Run AggregateRanger
./AggregateRanger.linux OutputFiles1/rangerOutput >> AggregateOutput1.txt