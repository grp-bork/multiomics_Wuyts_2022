pacman::p_load(tidyverse)
pacman::p_load(GGally)

# Load data
# Read in taxonomy information
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

# Read in datasets with species estimators
datasets <- tibble(
  omic = c("ampliconsequencing_currentpipeline", 
           "metagenomics_currentpipeline", "metagenomics_metaphlan3.0", "metagenomics_motus2.5",
           "metatranscriptomics_currentpipeline", "metatranscriptomics_metaphlan3.0", "metatranscriptomics_motus2.5",
           "metaproteomics_currentpipeline"),
  
  counts = list(
    read_tsv("../data/relative_abundance_tables/amplicon_sequencing_relabun.tsv.gz") %>% 
      select(-accession_number),
    read_tsv("../data/relative_abundance_tables/metagenomics_relabun.tsv.gz") %>% 
      select(-accession_number),
    read_tsv("../data/relative_abundance_tables/metagenomics_metaphlan3.0.tsv.gz"),
    read_tsv("../data/relative_abundance_tables/metagenomics_motus2.5_relabun.tsv.gz"),
    read_tsv("../data/relative_abundance_tables/metatranscriptomics_relabun.tsv.gz") %>% 
      select(-accession_number),
    read_tsv("../data/relative_abundance_tables/metatranscriptomics_metaphlan3.0.tsv.gz"),
    read_tsv("../data/relative_abundance_tables/metatranscriptomics_motus2.5_relabun.tsv.gz"),
    read_tsv("../data/relative_abundance_tables/metaproteomics_relabun.tsv.gz") %>% 
      select(-accession_number)
  )
  
)

# Subset to only keep shared samples

all_samples <- datasets %>% 
  mutate(samples = map(counts, colnames)) %>% 
  pull(samples)
shared_samples <- Reduce(intersect, all_samples)

datasets <- datasets %>% 
  mutate(shared_samples_long = map(counts, ~ .x %>% 
                                     pivot_longer(-species_name, names_to = "sample", values_to = "relabun") %>% 
                                     filter(sample %in% shared_samples)
                                )
  )

pdf("Supp_Figure_3.pdf", width=12, height=12)

datasets %>% 
  select(omic, counts) %>% 
  unnest(counts) %>% 
  pivot_longer(c(-omic, - species_name), names_to = "sample", values_to = "relabun") %>% 
  pivot_wider(names_from = omic, values_from = "relabun") %>% 
  select(-species_name, -sample) %>% 
  ggpairs() +
  theme(strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90)) +
  ggtitle("Comparison on species level")

dev.off()
