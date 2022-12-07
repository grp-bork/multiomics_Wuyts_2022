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
species_metadata <- read_tsv("../data/metadata/species_annotations.tsv.gz") %>%
  rename(accession_number = NT_code) %>% 
  left_join(genomes)

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

## Scatter plots
#
# In this visualisation we will plot all species within shared samples of 2-omics on 2 axes. If the omics agree on their relative abundance, they should follow the diagonal perfectly. I choose metagenomics as reference omic and compared all omics to metaG.
#
# Here I will define a function that I will use later on to repeat the plotting:

generate_scatter_plot <- function(x, y, z, y_abundance_cutoff, custom_df = count_data, log = T){
  
  # Extract axis 1
  x_axis <- custom_df %>% 
    filter(omic == x) %>%
    select(omic, relabun_shared_samples) %>% 
    unnest(cols = relabun_shared_samples) 
  
  # Extract axis 2
  y_axis <- custom_df %>% 
    filter(omic == y) %>%
    select(omic, relabun_shared_samples) %>% 
    unnest(cols = relabun_shared_samples)
  
  # Construct dataframe for plotting
  df <- x_axis %>% 
    bind_rows(y_axis) %>% 
    mutate(condition = str_extract(sample, "(?<=_)[A-z]+(?=_)")) %>% 
    left_join(genomes) %>% 
    left_join(species_metadata) %>% 
    pivot_wider(names_from = omic, values_from = relabun)
  
  # Set links to allow usage in ggplot2
  x <- sym(x)
  y <- sym(y)
  
  # Add y_abundance for different grouping
  if(!missing(y_abundance_cutoff)){
    new_var <- str_c("abundance_", y)
    df <- df %>% 
      mutate({{new_var}} := if_else({{y}} < y_abundance_cutoff, str_c("< ", y_abundance_cutoff), str_c("> ", y_abundance_cutoff)))
  }

  if (log == T) {
    # Calculate Spearman corr on log transformed abundances
    df_r2 <- df %>%
      mutate(x_log10 = log10(!!x),
             y_log10 = log10(!!y)) %>%
      # Correlation does not work, so remove these
      filter(x_log10 != -Inf,
             y_log10 != -Inf)
    
    correlation <-
      cor(df_r2 %>% pull(x_log10), df_r2 %>% pull(y_log10), method="spearman")
  } else {
    # Remove missing values
    df_r2 <- df %>% 
      filter(!is.na(!!x), !is.na(!!y))
    
    correlation <-
      cor(df_r2 %>% pull({{y}}),df_r2 %>% pull({{x}}), method="spearman")
  }
  
  correlation <- correlation %>% round(digits = 3)

  # Check if a colouring condition is present
  if(missing(z)){
    p <- ggplot(df, aes(!!x, !!y))
  } else {
    z <- sym(z)
    p <- ggplot(df, aes(!!x, !!y, colour = !!z))
  }
  
  # Construct final plot
  p <- p +
    geom_point() +
    expand_limits(x = 0, y = 0) +
    geom_abline(intercept = 0, colour = "red") +
    scale_color_brewer(palette = "Paired")
  
  # Convert to log if wanted
  if (log == T){
    p <- p +
      scale_x_log10(limits = c(1e-06,1)) +
      scale_y_log10(limits = c(1e-06,1)) +
      xlab(str_c(x, " (log)")) +
      ylab(str_c(y, " (log)"))
  } else {
    p <- p +
      scale_x_continuous(limits = c(1e-06,1)) +
      scale_y_continuous(limits = c(1e-06,1))
    
  }

    p +
      labs(caption = str_c("Spearman corr: ", correlation))
  
}

## Non-coloured plots 

per_species_plots <-
  function(x, y){
    
    # Extract axis 1
    x_axis <- count_data %>% 
      filter(omic == x) %>%
      select(omic, relabun_shared_samples) %>% 
      unnest(cols = relabun_shared_samples) 
    
    # Extract axis 2
    y_axis <- count_data %>% 
      filter(omic == y) %>%
      select(omic, relabun_shared_samples) %>% 
      unnest(cols = relabun_shared_samples)
    
    # Construct dataframe for plotting
    df <- x_axis %>% 
      bind_rows(y_axis) %>% 
      mutate(condition = str_extract(sample, "(?<=_)[A-z]+(?=_)")) %>% 
      left_join(genomes) %>% 
      pivot_wider(names_from = omic, values_from = relabun)
    
    # Set links to allow usage in ggplot2
    x <- sym(x)
    y <- sym(y)
    
    df %>% 
      ggplot(aes(x = !!x, y = !!y, colour = species_name)) +
      geom_point() +
      scale_x_log10() +
      scale_y_log10() + 
      gghighlight() +
      facet_wrap(~species_name)
    
  }

scatterplots <-
  expand_grid(x_omic = count_data$omic, y_omic = rev(count_data$omic)) %>% 
  filter(x_omic == "metagenomics" & y_omic == "amplicon_sequencing") %>%
  mutate(plot_species = map2(x_omic, y_omic,~per_species_plots(.x, .y)))

scatterplots <- scatterplots %>% 
  mutate(plot_condition = map2(x_omic, y_omic, ~generate_scatter_plot(.x, .y, "condition", F, log = T)),
         plot_gram = map2(x_omic, y_omic, ~generate_scatter_plot(.x, .y, "Gram_stain", F, log = T))
  )

pdf("Supp_Figure_4_A.pdf", width=12, height=12)

scatterplots$plot_species

dev.off()

pdf("Supp_Figure_4_B.pdf")

scatterplots$plot_condition

dev.off()

pdf("Supp_Figure_4_C.pdf")

scatterplots$plot_gram[[1]] + scale_colour_brewer(palette = "Dark2")

dev.off()
