# Description

The iRIGS (integrative risk gene selector) is a framework to identify risk genes from GWAS loci.


# installation

## First install the dependent package in R

	library("devtools")
	install_github("crotoc/iRIGS")

## Second copy the binary iRIGS to your PATH

The binary iRIGS in this repository is the command to run iRIGS

# Quick start

## Download resource file

To run iRIGS, at least we need:

--cand_loci: Specify a candidate loci list
--network: A precomputed network files saved in RData. 
--gene_symbol: A file including the start and end of all genes.
--generic_evi_file: formated features. Will use in burnin round.
--extra_evi_file: formated features. Will use in after-burnin round.

All the files used in our study can be downloaded here:

	https://vanderbilt365-my.sharepoint.com/:f:/g/personal/rui_chen_1_vanderbilt_edu/EmTLMwVI5DxLtN0dz8ITJJoBeBAwk7Yw6f6vyTE36N85Mg?e=RTWEOa

## network only mode 
	
	Gibbs --cand_loci scz.145.sub.1e-6.5e5.bed.fmt \
		  --network go_propogation_probality_rp_0.3.RData \
		  --gene_symbol 53934_transcribed_genes_GencodeV12_with_official_symbol \
		  --threads 20 \
		  --res_prefix output \
		  --res_path gibbs_result \
		  --distance TRUE \
		  --burnin_round 1000\
		  --after_burnin_round 1000 \
		  --max_gene 20 \
		  --evaluate_region FALSE \
		  --window_size 1e6



## feature integrating mode

	Gibbs --cand_loci scz.145.sub.1e-6.5e5.bed.fmt \
		  --network go_propogation_probality_rp_0.3.RData \
		  --gene_symbol 53934_transcribed_genes_GencodeV12_with_official_symbol \
	      --generic_evi_file generic.evi \
		  --extra_evi_file generic.evi \
		  --threads 20 \
		  --res_prefix output \
		  --res_path gibbs_result \
		  --distance TRUE \
		  --burnin_round 1000 \
		  --after_burnin_round 1000 \
		  --max_gene 20 \
		  --evaluate_region FALSE \
		  --window_size 1e6

 

