pacman::p_load(tidyverse)
pacman::p_load(clusterProfiler)
pacman::p_load(patchwork)
pacman::p_load(ggrepel)

organism_names_all_levels <-
  readxl::read_xlsx("../data/phylogeny/Lineage_of_32_sp.xlsx",
                          skip = 2,
                          col_names = c(
                            "ID",
                            "kingdom_gtdb",
                            "phylum_gtdb",
                            "class_gtdb",
                            "order_gtdb",
                            "family_gtdb",
                            "genus_gtdb",
                            "species_gtdb",
                            "nearest_id_in_GTDB",
                            "reference_genomes_in_use",
                            "kingdom_ncbi",
                            "phylum_ncbi",
                            "class_ncbi",
                            "order_ncbi",
                            "family_ncbi",
                            "genus_ncbi",
                            "species_ncbi"
                          )) %>%
  rename(accession_number = ID,
         species_name = species_gtdb)

organism_names <- organism_names_all_levels %>%
  select(accession_number, species_name)

# Load eggnog data
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
  mutate(feature = str_extract(query_name, "(?<=\\|)[A-z]+_[0-9]+(?=|)")) %>%
  mutate(COG = str_extract(eggNOG_OGs, "COG[0-9]{4}@2")) %>%
  filter(!is.na(feature), !is.na(COG))

gene_metadata <- read_tsv("../data/metadata/32_genomes.prokka.eggnog_annotations.gff.gz",
                            comment = "#",
                            col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")) %>%
  select(attribute) %>%
  mutate(accession_number = str_extract(attribute, "(?<=ID=)[^.]+"),
         product = str_extract(attribute, "(?<=product=).+$"),
         gene_name = str_extract(attribute, "(?<=locus_tag=)[^;]+(?=;)" )) %>%
  select(-attribute)  %>%
  left_join(organism_names) %>%
  left_join(eggnog %>% rename(gene_name = feature))

metaT_differential_abundant_all <- read_tsv("../data/differential_expression/DE_genes_and_proteins_unfiltered.tsv.gz") %>% 
  filter(omic == "metatranscriptomics")

eggnoglongnames <- read_tsv("../data/annotations/eggnogv4.funccats.txt.gz") %>% 
  rename(COG_Functional_cat = letter)

## Build a df of COG to category mapping
eggnog_categories_dereplicated <- eggnog %>%
  select(COG, COG_Functional_cat) %>%
  mutate(COG_Functional_cat = str_replace_all(COG_Functional_cat, "(?<=.)(?=.)", ",")) %>%
  separate_rows(COG_Functional_cat, sep = ",") %>%
  filter(!is.na(COG_Functional_cat)) %>%
  select(COG_Functional_cat, COG) %>% 
  distinct()

## Build a df of COG to gene mapping
eggnog_COGs_dereplicated <- eggnog %>% 
  select(feature, COG) %>% 
  filter(!is.na(COG)) %>% 
  rename(geneID = feature) %>% 
  select(COG, geneID)

# Define functions
perform_COG_enrichments <- function(dataset){
  # Filter for upregulated genes
  upreg_genes <- dataset %>% 
    filter(log2FoldChange > 2) %>% 
    left_join(eggnog %>% rename(gene = feature)) %>% 
    filter(!is.na(COG)) %>% 
    pull(gene)

  # Filter T15 for downregulated genes
  downreg_genes <- dataset %>% 
    filter(log2FoldChange < 2) %>% 
    left_join(eggnog %>% rename(gene = feature)) %>% 
    filter(!is.na(COG)) %>% 
    pull(gene)

  # Perform analyses
  upreg_enriched <- enricher(gene = upreg_genes, 
                   TERM2GENE = eggnog_COGs_dereplicated,
                   minGSSize = 10) %>% 
    as.data.frame() %>% 
    mutate(direction = "upregulated")

  downreg_enriched <- enricher(gene = downreg_genes, 
                 TERM2GENE = eggnog_COGs_dereplicated,
                 minGSSize = 10) %>% 
    as.data.frame() %>% 
    mutate(direction = "downregulated")
  
  # Merge results
  upreg_enriched %>% 
    bind_rows(downreg_enriched) %>% 
    arrange(-Count)
}

bacteroidota <- organism_names_all_levels %>% 
  filter(phylum_gtdb == "Bacteroidota") %>% 
  pull(accession_number)

# Perform analysis
COG_enrichment_analysis_species <- metaT_differential_abundant_all %>% 
  # Keep DEG only
  filter(abs(log2FoldChange) > 2,
         padj < 0.001,
         drug == "chlorpromazine",
          timepoint %in% c("timepoint15", "timepoint30", "timepoint60", "timepoint180")) %>%
  mutate(bacteroidota = if_else(accession_number %in% bacteroidota, "Bacteroidota", "Non-Bacteroidota")) %>% 
  # Make a df by drug and timepoint
  group_by(drug, timepoint, bacteroidota) %>% 
  nest() %>% 
  mutate(enriched_COGs = map(data, perform_COG_enrichments)) 

plot_data_species <- COG_enrichment_analysis_species %>% 
  ungroup() %>% 
  select(-data) %>% 
  unnest(enriched_COGs) %>% 
  left_join(eggnog_categories_dereplicated %>% rename(Description = COG)) %>% 
  left_join(eggnoglongnames) %>% 
  mutate(Count = if_else(direction == "downregulated", -Count, Count)) %>% 
  mutate(timepoint = factor(timepoint, levels = c("timepoint15", "timepoint30", "timepoint60", "timepoint180"))) %>% 
  group_by(timepoint, bacteroidota, description, direction) %>% 
  summarise(count = sum(Count)) %>% 
  ungroup() %>% 
  complete(timepoint, bacteroidota, description, direction, fill = list(count = 0))

order_species <- plot_data_species %>% 
  filter(timepoint == "timepoint15", bacteroidota == "Bacteroidota") %>% 
  group_by(description) %>% 
  summarise(summed_count = sum(count)) %>% 
  arrange(summed_count) %>% 
  pull(description) 

## Add all other levels to order if they do not exist in T15, to avoid conversion to NA
order_species <- c(unique(plot_data_species$description)[!unique(plot_data_species$description) %in% order_species], order_species)

pdf("Supp_Figure_13.pdf", height=8, width=12)

plot_data_species  %>% 
  mutate(description = factor(description, levels = order_species)) %>% 
  mutate(colouring = case_when(
    bacteroidota == "Bacteroidota" & direction == "upregulated" ~ "Bacteroidota: upregulated",
    bacteroidota == "Bacteroidota" & direction != "upregulated" ~ "Bacteroidota: downgulated",
    bacteroidota != "Bacteroidota" & direction == "upregulated" ~ "Non-Bacteroidota: upgulated",
    bacteroidota != "Bacteroidota" & direction != "upregulated" ~ "Non-Bacteroidota: downgulated")) %>% 
  ggplot(aes(x = count, y = description, fill = colouring, group = bacteroidota)) +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_col(position = "dodge", width = 0.5) + 
  facet_grid(~timepoint, scales = "free_x", space = "free") +
  ggtitle("Enriched COGs in differentially expressed genes") +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4", "lightgrey", "darkgrey")) +
  theme_bw() + 
  theme(strip.text.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()
        ) 

dev.off()
