# download all necessary data from Figshare

all: \
	phylogeny/bac120_r95.tree \
	ancom/ancom_FC_16S_time_conditions.txt.gz \
	ancom/ancom_FC_metaG_time_conditions.txt.gz \
	ancom/ancom_FC_metaP_time_conditions.txt.gz \
	ancom/ancom_FC_metaT_time_conditions.txt.gz \
	annotations/out.annot.2.0.1.emapper.annotations.gz \
	annotations/eggnogv4.funccats.txt.gz \
	colours/colour_scheme.tsv.gz \
	differential_expression/DE_genes_and_proteins_unfiltered.tsv.gz \
	differential_expression/DE_genes_LRT_unfiltered.tsv.gz \
	gene_count_tables/metagenomics_gene_counts_normalised_ko.tsv.gz \
	gene_count_tables/metagenomics_gene_counts_raw.tsv.gz \
	gene_count_tables/metaproteomics_gene_counts_normalised_ko.tsv.gz \
	gene_count_tables/metaproteomics_gene_counts_raw.tsv.gz \
	gene_count_tables/metatranscriptomics_gene_counts_normalised_ko.tsv.gz \
	gene_count_tables/metatranscriptomics_gene_counts_normalised_per_species.tsv.gz \
	gene_count_tables/metatranscriptomics_gene_counts_raw.tsv.gz \
	growth_curve/values.txt.gz \
	kegg_pathways/kegg_pathway_compound.tsv.gz \
	kegg_pathways/kegg_pathway_coverage_per_sample_metaB_specific_filteredann.tsv.gz \
	kegg_pathways/kegg_pathway_coverage_per_sample_metaB.tsv.gz \
	kegg_pathways/kegg_pathway_coverage_per_sample_metaG.tsv.gz \
	kegg_pathways/kegg_pathway_coverage_per_sample_metaP.tsv.gz \
	kegg_pathways/kegg_pathway_coverage_per_sample_metaT.tsv.gz \
	kegg_pathways/kegg_pathwayCLASS_coverage_per_reaction_per_sample_metaG.tsv.gz \
	kegg_pathways/kegg_pathwayCLASS_coverage_per_reaction_per_sample_metaP.tsv.gz \
	kegg_pathways/kegg_pathwayCLASS_coverage_per_reaction_per_sample_metaT.tsv.gz \
	kegg_pathways/kegg_pathwayCLASS_coverage_per_sample_metaB_specific.tsv.gz \
	kegg_pathways/kegg_pathways_filtered_by_community_members.csv.gz \
	kegg_pathways/keggCLASS_num_reactions.tsv.gz \
	kegg_pathways/keggCLASS_nummet.tsv.gz \
	kegg_pathways/keggPTWenrichment_fischer_metaB.csv.gz \
	kegg_pathways/keggPTWenrichment_fischer_metaP.csv.gz \
	kegg_pathways/keggPTWenrichment_fischer_metaT.csv.gz \
	metabolomics/metabolomics_joint_ann_filtered.tsv.gz \
	metabolomics/metabolomics_jointdata_quantnorm_ann_filtered.tsv.gz \
	metabolomics/metabolomics_negative_annotations.tsv.gz \
	metabolomics/metabolomics_negative_ions_abundance.tsv.gz \
	metabolomics/metabolomics_negative_ions_abundance.without_bacteria.tsv.gz \
	metabolomics/metabolomics_positive_annotations.tsv.gz \
	metabolomics/metabolomics_positive_ions_abundance.tsv.gz \
	metabolomics/metabolomics_positive_ions_abundance.without_bacteria.tsv.gz \
	metadata/32_genomes.prokka.eggnog_annotations.gff.gz \
	metadata/all_sample_metadata.tsv.gz \
	metadata/species_annotations.tsv.gz \
	pathway_enrichment_analysis/species_coverage_order_updated_extendedcolor.txt.gz \
	pathway_enrichment_analysis/species_tree_order_updated_extendedcolor.txt.gz \
	ranova/FC.metaP.tsv.gz \
	ranova/ranova.metaB.runAB.post_hoc_ttest.tsv.gz \
	ranova/ranova.metaP.runAB.post_hoc_ttest.tsv.gz \
	ranova/ranova.separate.metaB.log10.runAB.tsv.gz \
	ranova/ranova.separate.metaP.log10.runAB.tsv.gz \
	relative_abundance_tables/amplicon_sequencing_relabun.tsv.gz \
	relative_abundance_tables/metagenomics_metaphlan3.0.tsv.gz \
	relative_abundance_tables/metagenomics_motus2.5_relabun.tsv.gz \
	relative_abundance_tables/metagenomics_relabun.tsv.gz \
	relative_abundance_tables/metaproteomics_relabun.tsv.gz \
	relative_abundance_tables/metatranscriptomics_metaphlan3.0.tsv.gz \
	relative_abundance_tables/metatranscriptomics_motus2.5_relabun.tsv.gz \
	relative_abundance_tables/metatranscriptomics_relabun.tsv.gz \
	relative_abundance_tables/rel_abund_colour_scheme.tsv.gz \
	seq_stats/seq_stat.tsv.gz \
	rm_data_zip

%.gz: | data.zip
	mkdir -p $(dir $@)
	unzip data.zip $(notdir $@) -d $(dir $@)

phylogeny/bac120_r95.tree:
	cd phylogeny && make

data.zip: 
	wget https://figshare.com/ndownloader/articles/21667763/versions/2 -O data.zip

rm_data_zip:
	rm data.zip

clean:
	-rm */*gz phylogeny/bac120_r95.tree
	-rmdir *

.PHONY: clean rm_data_zip
