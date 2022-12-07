pacman::p_load(tidyverse)
pacman::p_load(patchwork)

# Read in data
## Transcript counts and annotation

# Get all transcripts that are differential abundant in niclosamide
## get p-value from over all timepoints
metaT_differential_abundant_pval <- read_tsv("../data/differential_expression/DE_genes_LRT_unfiltered.tsv.gz") %>% 
  filter(padj <= 0.05) 

## Get log2foldchange from per timepoint analysis
## Only use the first few timepoints
metaT_differential_abundant_all <- read_tsv("../data/differential_expression/DE_genes_and_proteins_unfiltered.tsv.gz") %>% 
  filter(omic == "metatranscriptomics")
metaT_differential_abundant_log2fold <- metaT_differential_abundant_all %>% 
  filter(timepoint %in% c("timepoint15", "timepoint30", "timepoint60", "timepoint180")) %>%   
  group_by(gene, drug) %>% 
  summarise(max_log2fc = suppressWarnings(max(abs(log2FoldChange), na.rm = T))) %>% 
  filter(max_log2fc >= 2) %>% 
  select(gene, drug, max_log2fc)

metaT_differential_abundant <- inner_join(metaT_differential_abundant_pval, metaT_differential_abundant_log2fold) 

# Counts
metaT_count <- read_tsv("../data/gene_count_tables/metatranscriptomics_gene_counts_raw.tsv.gz") %>%   # Filter on differential abundant genes only
  filter(gene_name %in% unique(metaT_differential_abundant$gene)) %>% 
  pivot_longer(c(-species_name, -gene_name), names_to = "sample", values_to = "count") %>% 
  separate(sample, into = c("run", "drug", "timepoint", "repeat")) %>% 
  filter(!is.na(`repeat`)) %>% 
  mutate(timepoint = factor(timepoint, levels = c("T0", "T15", "T30", "T1h", "T3h", "T48h", "T96h")))

# Curated list of pathways
kegg_pathways_curated <- read_tsv("../data/kegg_pathways/kegg_pathways_filtered_by_community_members.csv.gz") %>% 
  mutate(ID = str_replace(ID, "path:", "")) %>% 
  # Filter out the Gloval and overview maps
  filter(! (str_detect(ID, "map011") | str_detect(ID, "map012"))) %>% 
  rename(KEGG_Pathway = ID)

# Read in eggnog annotations for metaT to map genes to KEGG pathways
eggnog <-
  read_tsv(
    "../data/annotations/out.annot.2.0.1.emapper.annotations.gz",
    comment = "#",
    col_names = c(
      "query_name",
      "seed_eggNOG_ortholog",
      "seed_ortholog_evalue",
      "seed_ortholog_score",
      "best_tax_level",
      "Preferred_name",
      "GOs",
      "EC",
      "KEGG_ko",
      "KEGG_Pathway",
      "KEGG_Module",
      "KEGG_Reaction",
      "KEGG_rclass",
      "BRITE",
      "KEGG_TC",
      "CAZy",
      "BiGG_Reaction",
      "taxonomic_scope",
      "eggNOG_OGs",
      "best_eggNOG_OG",
      "COG_Functional_cat",
      "eggNOG_free_text_desc"
    )
  ) %>% 
  mutate(feature = str_extract(query_name, "(?<=\\|)[A-z]+_[0-9]+(?=|)")) 

eggnog_KEGG_annotation <- eggnog %>% 
  select(feature, KEGG_Pathway) %>% 
  filter(!is.na(feature), !is.na(KEGG_Pathway)) %>% 
  separate_rows(KEGG_Pathway, sep = ",") %>% 
  filter(str_detect(KEGG_Pathway, "map")) %>% 
  left_join(kegg_pathways_curated) %>% 
  filter(!is.na(Name)) 

# Niclosamide => Nitrogen metabolism
transcripts_of_interest <- eggnog %>% 
  mutate(OG = str_extract(eggNOG_OGs, "COG[0-9]{4}@2")) %>% 
  filter(OG %in% c("COG0334@2", "COG1151@2", "COG0549@2")) %>% 
  filter(!is.na(feature)) %>% 
  select(feature, OG)

pdf("Supp_Figure_9.pdf", width=8, height=4)

metaT_differential_abundant_all %>% 
  filter(drug == "niclosamide") %>% 
  left_join(transcripts_of_interest %>% rename(gene = feature)) %>% 
  filter(!timepoint %in% c("timepoint2880", "timepoint5760")) %>% 
  mutate(OG = case_when(
    OG == "COG0334@2" ~ "NAD-specific glutamate dehydrogenase (COG0334)",
    OG == "COG1151@2" ~ "Hydroxylamine reductase (COG1151)",
    OG == "COG0549@2" ~ "Carbamate kinase (COG0549)"
  ))  %>% 
  mutate(timepoint = factor(timepoint, levels = c("timepoint15", "timepoint30", "timepoint60", "timepoint180"))) %>% 
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = OG)) +
  geom_point(data = function(x) subset(x, is.na(OG)), alpha = 0.6) +
  geom_point(data = function(x) subset(x, !is.na(OG)), size = 2) +
  facet_wrap(~timepoint, scales = "free") +
  scale_colour_brewer(palette = "Dark2", na.value = "grey") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  ggtitle("COGs involved in Nitrogen metabolism")

dev.off()
