## Quality control 
## Genome assembly
## RepeatMasker estimation
We use search engines — RMBlast, library — Repbase and other default parameters to screen Apis laboriosa and Apis dorsata genomes for interspersed repeats and low complexity DNA sequences. 

	RepeatMasker -parallel 20 -species honeybee -gff -dir laboriosa_repeatmasker Apis_laboriosa.fasta
	RepeatMasker -parallel 20 -species honeybee -gff -dir dorsata_repeatmasker Apis_dorsata.fasta
	
## Gene prediction
## Gene family construction
	# Count the number of all types of gene families
	python count_all_gene_family --gene_family_directory /home/lin/gene_families/ --file_suffix fasta --gene_family_number 10460

## Multiple alignment
    # Run muscle
	ls *.pep.fasta|awk '{print "muscle -in "$0" -out "$1".pep.muscle"}' | sh
	ls /home/lin/gene_families | grep "pep.fasta$" | awk -F'.' '{print "perl pepMfa_to_cdsMfa.pl "$1".pep.muscle "$1".fasta > "$1".gene_family.fasta"}' | sh
	
## Phylogenetic analysis
### Genetic distance estimation
	# Change fasta format into phy format for all single-copy gene file
	ls *.conserved_gene_family_for_genetic_distance.fasta | awk '{print "perl fasta2philip.pl "$0" > "$0".phy "}' | sh
	# Calculate genetic distance
	calculate_genetic_distance.r --gene_family_directory /home/lin/gene_families/ 
### Species tree construction
	# Concatenate all single-copy gene families into a fasta file
	cat *.conserved_gene_family.fasta >> all_conserved_gene_family.fasta
	# Change fasta format into phy format
	perl fasta2philip.pl all_conserved_gene_family.fasta > 	all_conserved_gene_family.phy
	# Run phyml
	phyml -i all_conserved_gene_family.phy -b 1000 -d nt -m GTR -c 4 -a e -o tlr
###  Divergence time estimation
## Positive selection analysis
### PAML free-ratio model

"The branch models allow the ω (dN/dS, the synonymous rate / the non-synonymous rate) to vary among branches in the phylogeny and are useful for detecting positive selection acting on particular lineages."

#### Step 1: Define phylogenetic tree
As we assume no clock and rates are entirely free to vary from branch to branch (clock = 0 at the ctl file), the unrooted tree is used under this model. The tree is shown as follows:
		
	(Temnothorax, Dinoponera, ((Florea, ((Laboriosa, Dorsata), (Cerana, Mellifera))), (Terrestris, Impatiens)));

#### Step 2: Prepare the control(ctl) file
We use the free-ratios model (model = 1) to compute an independent ω ratio for each branch. The ctl file is set as follows:
     
      seqfile = *.conserved_gene_family.fasta
     treefile = *.conserved_gene_family.nwk
      outfile = *.conserved_gene_family.free_ratio.mlc

        noisy = 3   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        model = 1   * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
      NSsites = 0   * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                    * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                    * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                    * 13:3normal>0
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1   * initial or fixed omega, for codons or codon-transltd AAs
    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 4   * # of categories in the dG or AdG models of rates
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 	RateAncestor = 1   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
 	 fix_blength = 0    * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
 	 
 Run the free_ratio.ctl file:
 
 	ls *.conserved_gene_family.free_ratio_branch.ctl | awk '{print "codeml "$0""}' | sh

### PAML branch-site model

"The branch-site models allow ω to vary both among sites in the protein and across branches on the tree and aim to detect positive selection affecting a few sites along particular lineages."

#### Step 1: Define phylogenetic tree
We need to define the lineages of interest (called foreground branches) in the phylogentic tree. Different ω can apply to the foreground branches, while  the background branches share the same distribution of ω value among sites.

We set Apis laboriosa and Apis dorsata as the foreground branch, respectively. The trees are shown as follows:

	# Set Apis laboriosa lineage as foreground branch
	(Temnothorax, Dinoponera, ((Florea, ((Laboriosa #1, Dorsata), (Cerana, Mellifera))), (Terrestris, Impatiens)));
	# Set Apis dorsata lineage as foreground branch
	(Temnothorax, Dinoponera, ((Florea, ((Laboriosa, Dorsata #1), (Cerana, Mellifera))), (Terrestris, Impatiens)));
	
#### Step 2: Prepare the control(ctl) files
To calculate the likelihood value, we need to compare an alternative model (model = 2 NSsites = 2), where the foreground branch could have a proportion of sites under positive selection, with the corresponding null model (fix_omega = 1 and omega = 1), in which the foreground branch may have different proportions of sites under neutral selection and purifying selection to the background branches (see the following table).

The ctl file of the alternative mdoel is set as follows:

*model = 2 (different dN/dS for branches)

*NSsites = 2 (3 categories for sites: purifying, neutral and positive selection).
	      		
	  seqfile = *.conserved_gene_family.fasta
     treefile = *.conserved_gene_family.nwk
      outfile = *.conserved_gene_family.alternative.mlc

        noisy = 3   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        model = 2  * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
      NSsites = 2   * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                    * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                    * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                    * 13:3normal>0
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1   * initial or fixed omega, for codons or codon-transltd AAs
    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0   * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 4   * # of categories in the dG or AdG models of rates
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 		RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
  		fix_blength = 0    * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
  		
 Run the alternative.ctl file:
 
 	ls  *.conserved_gene_family.branch_site_alternative.ctl |awk '{print "codeml "$0""}' | sh

The ctl file of the null mdoel is set as follows:

*fix_omega = 1 (dN/dS is fixed to 1)

*omega = 1

	  seqfile = *.conserved_gene_family.fasta
     treefile = *.conserved_gene_family.nwk
      outfile = *.conserved_gene_family.null.mlc

        noisy = 3   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        model = 2  * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
      NSsites = 2   * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                    * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                    * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                    * 13:3normal>0
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
    fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1   * initial or fixed omega, for codons or codon-transltd AAs
    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0   * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 4   * # of categories in the dG or AdG models of rates
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 		RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
  		fix_blength = 0    * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
  		
  Run the null.ctl file:

 	ls  *.conserved_gene_family.branch_site_null.ctl |awk '{print "codeml "$0""}' | sh
  		
#### Step 3: Calculate the p value

We retieve the likelihood values lnL1 and lnL0 from the alternative.mlc and the null.mlc files, respectively. Then, we construct the ΔLRT = 2×(lnL1 - lnL0). Finally, we run chi2 1  ΔLRT to get a p value. The final p value is this p value divided by 2.

#### Step 4: Retrieve sites under positive selection if statistically significant (p value < 0.05)

PAML uses the Bayes empirical Bayes (BEB) to calculate the posterior probabilities for site classes and identify sites under positive selection if the likelihood ratio test is significant. The program prints out an * if the posterior probability is >95%, and ** if the probability is > 99%.

### BUSTED (Branch-Site Unrestricted Statistical Test for Episodic Diversification) 
"BUSTED provides a gene-wide (not site-specific) test for positive selection by asking whether a gene has experienced positive selection at at least one site on at least one branch."

Unlike PAML, BUSTED assumes that positive selection may also be present in the background branches, allowing a more adequate description of complex biological data. 

#### Step 1: Define phylogenetic tree

We set Apis laboriosa and Apis dorsata as the foreground branch, respectively. The trees are shown as follows:

	# Set Apis laboriosa lineage as foreground branch
	(Temnothorax, Dinoponera, ((Florea, ((Laboriosa{Foreground}, Dorsata), (Cerana, Mellifera))), (Terrestris, Impatiens)));
	# Set Apis dorsata lineage as foreground branch
	(Temnothorax, Dinoponera, ((Florea, ((Laboriosa, Dorsata{Foreground}), (Cerana, Mellifera))), (Terrestris, Impatiens)));

#### Step 2: Run hyphy BUSTED
  
  	ls  *.conserved_gene_family.muscle.fasta |awk -F '.' '{print "hyphy busted --alignment "$0" --tree *.conserved_gene_family.nwk --branches Foreground >> "$1".conserved_gene_family.hyphy.log"}' | sh
	
BUSTED outputs the likelihood under the unconstrained alternative model (allows positive selection, neutral selection and purifying selection for both foreground and backgrpund branches), and under the constrained null model (disallows positive selection among foreground branches). Also, BUSTED outputs the p value for likelihood ratio test. If p value < 0.05,  we can reject the null hypothesis that there is no at least one site that has, at least some of the time, experienced positive selection on the foreground branches.

In addition to this output, BUSTED also calculates "Evidence Ratios" (ERs) for each site. The ER reports the likelihood ratio (in a log-scale format) that the unconstrained alternative model is a better fit to the data compared to the constrained null model. Unlike the BEB posterior probabilities provided by PAML, the ER for each site just provides descriptive information about whether a given site could have evolved under positive selection. 

## Gene duplication/loss analysis
### RANGERDTL v2.0

"RANGER-DTL (short for Rapid ANalysis of Gene family Evolution using Reconciliation-DTL) is a software package for inferring gene family evolution by speciation, gene duplication, horizontal gene transfer, and gene loss. The software takes as input a gene tree (rooted or unrooted) and a rooted species tree and reconciles the two by postulating speciation, duplication, transfer,
and loss events."

#### Step 1: Construct unrooted gene trees
	
	# Change fasta format into phy format
	ls *.gene_family_for_dl.fasta | awk '{print "perl fasta2philip.pl "$0" >  "$1".gene_family_for_dl.phy"}' | sh
	# Run phyml
	ls *.gene_family_for_dl.phy | awk '{print " phyml -i "$0" -b 100 -d nt -m GTR -c 4 -a e"}' | sh

#### Step 2: Compute the optimal roots for unrooted gene trees to minimize the DTL reconciliation cost

	ls  *.gene_family_for_dl.newick |awk -F '.' '{print "./OptRoot.linux -i "$0" -D 2 -T 1000 -L 1 | grep ";" >> "$1".gene_family_for_dl.opt.newick"}' | sh

#### Step 3: Compute an optimal DTL reconciliation for the rooted gene tree and rooted species tree and summarize the results

We use the following Bash script to compute 100 optimal reconciliations (sampled uniformly at random), and then aggregate the samples into a single reconciliation file (AggregateOutput.txt). 

We use default cost for duplication and loss, and set the HGT cost as 1,000 to exclude the HGT events.

	for j in *.gene_family_for_dl.opt.newick; do
		cat species.newick >> $j.gene_family.newick
		cat $j.gene_family.opt.newick >> $j.gene_family.newick
		mkdir OutputFiles$j
		for((i=1;i<=100;i++)); do
			./ranger-dtl-U.linux --seed $i -i $j.gene_family.newick -D 2 -T 1000 -L 1 -o OutputFiles$j/rangerOutput$i
		done
		./AggregateRanger.linux OutputFiles$j/rangerOutput >> AggregateOutput$j.txt
	done
	
### TREERECS

"Treerecs is a very easy-use phylogenetic software based on duplication-loss reconciliation."

	ls *.gene_family_for_dl.newick | awk -F '.' '{print "treerecs -g "$0" -s speices_tree_rooted.nwk -O newick:svg -O recphyloxml -t 98 -o "$1".treerecs_out -f"}' | sh
