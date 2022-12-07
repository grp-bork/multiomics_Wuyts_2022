pacman::p_load(tidyverse)
pacman::p_load(ggridges)

# Load data

count_data <- 
  tibble(
    omic = c("metagenomics", "metatranscriptomics", "metaproteomics"),
    counts = list(
      read_tsv("../data/gene_count_tables/metagenomics_gene_counts_raw.tsv.gz"),
      read_tsv("../data/gene_count_tables/metatranscriptomics_gene_counts_raw.tsv.gz"),
      read_tsv("../data/gene_count_tables/metaproteomics_gene_counts_raw.tsv.gz")
    )
  )

# Get the coverage

## Sanity check: how many genes per genome?

count_data <- count_data %>% 
  mutate(total_gene_count_species = map(counts,
                                        ~.x %>% 
                                          group_by(species_name) %>% 
                                          mutate(total_gene_count_species = n_distinct(gene_name)) %>% 
                                          select(species_name, total_gene_count_species) %>% 
                                          distinct()))

# Get total gene count per species
genes_per_genome <- count_data %>% 
  filter(omic == "metagenomics") %>% 
  select(total_gene_count_species) %>% 
  unnest(total_gene_count_species) 

count_data <- count_data %>% 
  mutate(genome_coverage = map(counts,
                               ~.x %>% 
                                 mutate_if(is.numeric, ~ if_else(.x > 0, 1, 0)) %>%
                                 mutate_if(is.numeric , replace_na, replace = 0) %>% 
                                 group_by(species_name) %>% 
                                 gather("sample", "gene_present", -species_name, -gene_name) %>% 
                                 left_join(genes_per_genome) %>% 
                                 group_by(species_name, sample) %>% 
                                 summarise(genes_with_count = sum(gene_present),
                                           genes_with_count_perc = sum(gene_present / total_gene_count_species)))
         
         
  )

# Create an order on highest metaG value
order <- count_data %>% 
  filter(omic == "metagenomics") %>% 
  select(genome_coverage) %>% 
  unnest(genome_coverage) %>%
  group_by(species_name) %>% 
  summarise(mean = mean(genes_with_count_perc)) %>% 
  arrange(mean) %>% 
  pull(species_name)

order_alternate_colour <- tibble(omic = c(rep("metagenomics", length(order)),
                                          rep("metatranscriptomics", length(order)),
                                          rep("metaproteomics", length(order))),
                                 species_name = rep(order, 3),
                                 linenumber = rep(1:length(order), 3)
                                 ) %>% 
  mutate(colour_variable = case_when(
    omic == "metagenomics" ~ "colour4",
    omic == "metatranscriptomics" ~ "colour5",
    omic == "metaproteomics" ~ "colour6",
  ))

svg("Figure_2_C1.svg")

count_data %>% 
  select(omic, genome_coverage) %>% 
  unnest(genome_coverage) %>%
  left_join(order_alternate_colour) %>% 
  mutate(species_name = factor(species_name, levels = order)) %>% 
  ggplot(aes(x = genes_with_count_perc, y=species_name, fill = colour_variable)) +
  geom_density_ridges(colour = "white", alpha = .8, from = 0, to = 1, scale = 2) +
  theme_ridges() +
  scale_fill_manual(breaks = c("colour4", "colour5", "colour6"),
                    values = c("#2b8cbe", "#e34a33", "#2ca25f", "#045a8d", "#b30000", "#006d2c"),
                    labels = c("colour4" = "metagenomics", "colour5" = "metatranscriptomics", "colour6" = "metaproteomics"),
                    name = "omic") +
  ylab("") +
  xlab("Fraction of species-specific proteins covered")

dev.off()
