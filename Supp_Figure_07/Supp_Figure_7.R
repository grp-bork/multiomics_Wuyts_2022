pacman::p_load(tidyverse)
pacman::p_load(rlang)
pacman::p_load(patchwork)
pacman::p_load(ggpubr)
pacman::p_load(RColorBrewer)
pacman::p_load(vegan)
pacman::p_load(gghighlight)

# Link genome accession numbers to species name
taxonomy <- 
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
                          )) 
genomes <- taxonomy %>% 
  rename(accession_number = ID,
         species_name = species_gtdb,
         phylum = phylum_gtdb,
         genus = genus_gtdb) %>% 
  select(accession_number, species_name, phylum, genus)

# Metadata table
metadata <- read_tsv("../data/metadata/all_sample_metadata.tsv.gz")

colour_scheme <- read_tsv("../data/colours/colour_scheme.tsv.gz")

count_data <- tibble(
  omic = c("amplicon_sequencing", "metagenomics", "metatranscriptomics", "metaproteomics"),
  relabun = list(
    read_tsv("../data/relative_abundance_tables/amplicon_sequencing_relabun.tsv.gz") %>% 
      select(-accession_number) %>% 
      pivot_longer(-species_name, names_to = "sample", values_to = "relabun"),
    read_tsv("../data/relative_abundance_tables/metagenomics_relabun.tsv.gz") %>% 
      select(-accession_number) %>% 
      pivot_longer(-species_name, names_to = "sample", values_to = "relabun"),
    read_tsv("../data/relative_abundance_tables/metatranscriptomics_relabun.tsv.gz") %>% 
      select(-accession_number) %>% 
      pivot_longer(-species_name, names_to = "sample", values_to = "relabun"),
    read_tsv("../data/relative_abundance_tables/metaproteomics_relabun.tsv.gz")  %>%     select(-accession_number) %>% 
      pivot_longer(-species_name, names_to = "sample", values_to = "relabun")
  )
)

count_data <- count_data %>% 
  mutate(n_samples = map_int(relabun, ~ .x$sample %>% unique %>% length),
         n_taxa = map_int(relabun, ~ .x$species_name %>% unique %>% length)
  )

samples_to_keep <- count_data %>%
  select(omic, relabun) %>% 
  unnest(cols = c(relabun)) %>% 
  select(omic, sample) %>% 
  distinct() %>% 
  group_by(sample) %>% 
  summarise(sample_abundance = n()) %>% 
  ungroup() %>% 
  filter(sample_abundance == max(sample_abundance)) %>% 
  pull(sample)

count_data <- count_data %>% 
  mutate(relabun_shared_samples = map(relabun,
                                      ~ .x %>% filter(sample %in% samples_to_keep)))

count_data <- count_data %>% 
  mutate(n_shared_samples = map_int(relabun_shared_samples, ~ .x$sample %>% unique %>% length)
  )

count_data <- count_data %>% 
  mutate(relabun_shared_samples = map(relabun_shared_samples,
                    ~ .x %>% 
                      separate(sample, into = c("run", "drug", "timepoint", "rep")) %>% 
                      group_by(species_name, run, drug, timepoint) %>% 
                      summarise(relabun = mean(relabun)) %>% 
                      ungroup() %>% 
                      mutate(sample = str_c(run, drug, timepoint, sep = "_")) %>% 
                      select(species_name, sample, relabun)))

mean_abundance <- count_data %>% 
  filter(omic == "metagenomics") %>%
  select(relabun) %>% 
  unnest(relabun) %>% 
  separate(sample, into = c("run", "drug", "time", "tech_rep"), remove = T) %>% 
  filter(!is.na(tech_rep)) %>% 
  group_by(run, drug, time, species_name) %>% 
  summarise(mean_abundance = mean(relabun))

# Get distance matrix
distmat <- mean_abundance %>% 
  ungroup() %>% 
  mutate(sample = str_c(run, drug, time, sep = "_")) %>% 
  select(-run, -time, - drug) %>% 
  pivot_wider(names_from = sample, values_from = mean_abundance) %>%
  column_to_rownames("species_name") %>% 
  # Calculate Bray-curtis distance
  as.matrix() %>% 
  t() %>% 
  vegdist(method = "bray", upper = T, diag = T)

pdf("Supp_Figure_7.pdf", width=16, height=16)

distmat %>%
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "sample1")  %>% 
  pivot_longer(-sample1, names_to = "sample2", values_to = "dist")%>% 
  separate(sample1, into = c("run1", "drug1", "time1"), remove = F) %>%
  separate(sample2, into = c("run2", "drug2", "time2"), remove = F) %>%
  # Calculate similarity
  ggplot(aes(x = sample1, y = sample2, fill = dist, label = round(dist, 2))) +
  geom_tile(colour = "lightgrey", size = 0.5) +
  # Add extra col for drugs
  geom_tile(aes(x = 35), fill = "white") +
  geom_point(aes(x = 35, colour = drug2), size = 4) +
  # Add extra row for drugs
  geom_tile(aes(y = 35), fill = "white") +
  geom_point(aes(y = 35, colour = drug1), size = 4) +
  # Add values in each box
  geom_text(size = 3, alpha = 0.8) + 
  
  
  expand_limits(fill = 1) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "",
    y = "",
    fill = "Bray-Curtis dissimilarity index"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    
    # Gridlines
    panel.grid = element_blank(),
    
    # Backgrounds
    panel.background = element_blank(),
    plot.background = element_blank(),
    
    # Facet titles
    strip.text.x = element_text(face = "bold", size = 14),
    strip.background = element_blank()
    
    )

dev.off()
