pacman::p_load(tidyverse)
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

normalised_count <- read_tsv("../data/gene_count_tables/metatranscriptomics_gene_counts_normalised_per_species.tsv.gz")

hs_plot <- normalised_count %>%
  filter(gene_name %in% c("AFJCKDGM_03987", "AFJCKDGM_04224", "AFJCKDGM_04223", "AFJCKDGM_04766", "AFJCKDGM_00012")) %>%
  pivot_longer(c(-species_name, -gene_name), names_to = "sample", values_to = "count_norm") %>%
  separate(sample, into = c("run", "drug", "timepoint", "repeat")) %>%
  filter(!str_detect(drug, "Transfer")) %>%
  left_join(gene_metadata %>% select(gene_name, Preferred_name, species_name)) %>%
  mutate(gene_name = str_c(Preferred_name, "\n", species_name, "\n", gene_name)) %>%
  mutate(timepoint = factor(timepoint, levels = c("T0", "T15", "T30", "T1h", "T3h", "T48h", "T96h"))) %>%
  ggplot(aes(x = timepoint, y = count_norm, colour = drug, group = drug)) +
  geom_point(position=position_dodge(width=0.3)) +
  stat_summary(fun=mean, geom="line") +
  facet_grid(gene_name~run, scales = "free_y") +
  theme(strip.text.y = element_text(angle = 0)) +
  ggtitle("Metatranscriptomics: heat shock response genes in Escherichia coli")

os_plot <- normalised_count %>%
  filter(gene_name %in% c("AFJCKDGM_04512", "AFJCKDGM_04492", "AFJCKDGM_04631", "AFJCKDGM_04630", "AFJCKDGM_04455", "AFJCKDGM_01789", "AFJCKDGM_01780", "AFJCKDGM_01867", "AFJCKDGM_00577", "AFJCKDGM_00578")) %>%
  pivot_longer(c(-species_name, -gene_name), names_to = "sample", values_to = "count_norm") %>%
  separate(sample, into = c("run", "drug", "timepoint", "repeat")) %>%
  filter(!str_detect(drug, "Transfer")) %>%
  left_join(gene_metadata %>% select(gene_name, Preferred_name, species_name)) %>%
  mutate(gene_name = str_c(Preferred_name, "\n", species_name, "\n", gene_name)) %>%
  mutate(timepoint = factor(timepoint, levels = c("T0", "T15", "T30", "T1h", "T3h", "T48h", "T96h"))) %>%
  ggplot(aes(x = timepoint, y = count_norm, colour = drug, group = drug)) +
  geom_point(position=position_dodge(width=0.3)) +
  stat_summary(fun=mean, geom="line") +
  facet_grid(gene_name~run, scales = "free_y") +
  theme(strip.text.y = element_text(angle = 0)) +
  ggtitle("Metatranscriptomics: oxidative stress response genes in Escherichia coli")

pdf("Supp_Figure_12.pdf", height=12, width=12)

hs_plot + os_plot + plot_layout(ncol=1) + plot_layout(guides = "collect", heights = c(0.4, 0.6))

dev.off()
