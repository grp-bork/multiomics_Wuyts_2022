pacman::p_load(tidyverse)
pacman::p_load(rlang)
pacman::p_load(patchwork)
pacman::p_load(ggpubr)
pacman::p_load(RColorBrewer)
pacman::p_load(vegan)
pacman::p_load(gghighlight)

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

# Number of taxa and samples

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

# Now finally, let's take the mean of the repeats to avoid overplotting

count_data <- count_data %>% 
  mutate(relabun_shared_samples = map(relabun_shared_samples,
                    ~ .x %>% 
                      separate(sample, into = c("run", "drug", "timepoint", "rep")) %>% 
                      group_by(species_name, run, drug, timepoint) %>% 
                      summarise(relabun = mean(relabun)) %>% 
                      ungroup() %>% 
                      mutate(sample = str_c(run, drug, timepoint, sep = "_")) %>% 
                      select(species_name, sample, relabun)))


# Relabun plots

## Technical repeats merged

relabun_plot_merged <- function(x, type = "overview") {
  
  # Extract data and merge technical repeats
  data <- count_data %>%
    filter(omic == x) %>%
    select(relabun) %>%
    unnest(relabun) %>% 
    mutate(sample = str_remove(sample, "_r[12]$")) %>% 
    group_by(species_name, sample) %>% 
    summarise(relabun = sum(relabun)) %>% 
    group_by(sample) %>% 
    mutate(relabun = relabun/sum(relabun)) %>%
    left_join(metadata %>%  select(-identifier,-rrna_pilot,-resequence, - batch, -replicate) %>% 
                mutate(sample = str_remove(sample, "_r[12]$")) %>%                 
                distinct()
              ) %>%
    filter(drug %in% c("control", "chlorpromazine", "niclosamide", "metformin")) 
  
  
  if (x == "metagenomics") {
    
    
    # Add timepoint 0 to all cases for metagenomics
    timepoint0 <- data %>% 
      filter(timepoint == 0,
             drug == 'control')
    
    for (z in  c("chlorpromazine", "niclosamide", "metformin")){
      data <- data %>% 
        bind_rows(
          timepoint0 %>% 
            mutate(drug = z)
        )
    }  
    
  }
  
  # Find most abundant species
  species_order <- data %>%
    group_by(species_name) %>%
    summarise(sum_relabun = sum(relabun)) %>%
    arrange(-sum_relabun) %>%
    pull(species_name)
  top12 <- species_order[1:12]
  
  if (type == "overview") { 
    
   plot <- data %>% 
      # Convert to factors for ordering of plot
      mutate(species = if_else(species_name %in% top12, species_name, "other"), 
             species = factor(species, levels = c(top12, "other")),
             species_name = factor(species_name, levels = species_order),
             drug = factor(drug, levels = c("control", "chlorpromazine", "metformin", "niclosamide")),
             run = factor(run, levels = c("A", "B"))
      )  %>% 
      # Plot
      ## Basics
      ggplot(aes(x = label, y = relabun, group = species_name)) +
      geom_area(aes(fill = species), colour = "black", size=0.2) +
      facet_grid(run ~ drug, scales = "free_x") +
      
      ## Make a scale of 13 colours based on a colourbrewer palette
      scale_fill_manual(values = 
                          c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                            '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                            '#cab2d6','#6a3d9a','#ffff99','#b15928',
                            "lightgrey")) +
      ## Only select 3 values, y-axis
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      
      ## Labels
      xlab("") +
      ylab("Relative abundance (%)\n") +
      ggtitle(str_c("Relative abundance (", x, ")")) +
      
      ## Customise lay-out
      theme_minimal() +
      theme(
        # Legend
        legend.position = "bottom",
        legend.title = element_blank(),
        
        # Axes labeling
        axis.text.y = element_text(size = 10, margin= margin(0,-15,0,0)),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9),
        
        # Gridlines
        panel.grid = element_blank(),
        
        # Spacing
        panel.spacing.x = unit(-2, "lines"),
        
        # Backgrounds
        panel.background = element_blank(),
        plot.background = element_blank(),
        
        # Facet titles
        strip.text.x = element_text(face = "bold", size = 14),
        strip.text.y = element_text(margin = margin(0,0,0,-1))
      )
  } else {
    
    plot <- data %>% 
      ggplot(aes(x = label, y = relabun, colour = drug)) +
      geom_line(aes(group = interaction(run, drug))) +
      geom_point() + 
      expand_limits(y = 0) +
      facet_wrap(~ species_name, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_color_brewer(palette = "Set2") +
      ylab("Relative abundance") +
      xlab("") +
      ggtitle(str_c("Relative abundance (", x, ")")) 
  }
  
  plot
}

colour_scheme_top12 <-
  colour_scheme %>% 
  group_by(colour) %>% 
  summarise(species_name = list(species_gtdb) %>% str_c(collapse = ";")) %>% 
  mutate(species_name = if_else(colour == "#d4d4d4", "other", species_name))

count_data_non_G <- filter(count_data, omic != "metagenomics")

pdf("Supp_Figure_6_C.pdf", height=6, width=12)

count_data_non_G %>% 
  select(omic) %>% 
  mutate(plots = map(omic,
                     ~relabun_plot_merged(.x) + geom_vline(aes(xintercept =  label), alpha = 0.4) +
                       scale_fill_manual(values = colour_scheme_top12$colour,
                                         breaks = colour_scheme_top12$species_name)
                     )
         ) %>% 
  pull(plots)

dev.off()

invsimpson_metaG <- count_data %>% 
  filter(omic == "metagenomics") %>%
  select(relabun_shared_samples) %>% 
  unnest(relabun_shared_samples) %>% 
  pivot_wider(names_from = sample, values_from = relabun) %>%
  column_to_rownames("species_name") %>% 
  # Calculate Bray-curtis distance
  as.matrix() %>% 
  t() %>% 
  diversity(index = "invsimpson")

pdf("Supp_Figure_6_AB.pdf", height=6, width=12)

tibble(
  sample = names(invsimpson_metaG),
  invsimpson = invsimpson_metaG
) %>% 
  separate(sample, into = c("run", "drug", "time", "tech_rep"), remove = T) %>% 
  ggplot(aes(x = time, y = invsimpson, group = drug, colour = drug)) +
  stat_summary(fun = "mean", geom = "line") +
  geom_point() +
  expand_limits(y = 0) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Alpha diversity",
    y = "Inverse simpson index",
    x = "",
    colour = "Drug"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~run, nrow = 2)

dev.off()
