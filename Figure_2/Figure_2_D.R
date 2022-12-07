pacman::p_load(tidyverse)
pacman::p_load(ggridges)
pacman::p_load(cowplot)

kegg_pathways <- read_tsv("../data/kegg_pathways/kegg_pathways_filtered_by_community_members.csv.gz")

coverage_max <- tibble(omic = c("metaG", "metaT", "metaP", "metaB")) %>%
  mutate(data = map(omic, function(x){
                      if(x != "metaB"){
                        read_tsv(
                          str_c(
                            "../data/kegg_pathways/kegg_pathway_coverage_per_sample_",
                            x,
                            ".tsv.gz"
                          )
                        )} else {
                        read_tsv(
                          "../data/kegg_pathways/kegg_pathway_coverage_per_sample_metaB_specific_filteredann.tsv.gz"
                        )
                        }
    })) %>%
  unnest(data) %>%
  rename(pathway = ...1) %>%
  pivot_longer(c(-omic,-pathway), names_to = 'sample', values_to = "coverage") %>%
  remove_missing() %>%
  mutate(
    omic = case_when(
      omic == "metaG" ~ "metagenomics",
      omic == "metaT" ~ "metatranscriptomics",
      omic == "metaP" ~ "metaproteomics",
      omic == "metaB" ~ "metabolomics",
    )
  ) %>%
  filter(pathway %in% kegg_pathways$ID) %>%
  left_join(kegg_pathways %>% rename(pathway = ID))

coverage_pathway_class_species_ind <- tibble(
  omic = c("metaG", "metaT", "metaP", "metaB")
) %>% 
  mutate(data = map(omic, ~ if (.x == "metaB") {
    read_tsv("../data/kegg_pathways/kegg_pathwayCLASS_coverage_per_sample_metaB_specific.tsv.gz")
  } else {
    read_tsv(str_c("../data/kegg_pathways/kegg_pathwayCLASS_coverage_per_reaction_per_sample_", .x, ".tsv.gz"))
  }
  )
  ) %>% 
  unnest(data) %>% 
  rename(pathway = ...1) %>% 
  pivot_longer(c(-omic, - pathway), names_to = 'sample', values_to = "coverage") %>% 
  remove_missing() %>% 
  mutate(omic = case_when(
    omic == "metaG" ~ "metagenomics",
    omic == "metaT" ~ "metatranscriptomics",
    omic == "metaP" ~ "metaproteomics",
    omic == "metaB" ~ "metabolomics",
  ))

order <- coverage_pathway_class_species_ind %>% 
  filter(omic == "metabolomics") %>% 
  group_by(pathway) %>% 
  summarise(mean = mean(coverage)) %>% 
  arrange(mean) %>% 
  pull(pathway)

coverage <- coverage_pathway_class_species_ind %>%
  mutate(pathway = factor(pathway, levels = order)) %>% 
  filter(!is.na(pathway)) %>% 
  ggplot(aes(x = coverage, y = pathway, fill = omic)) +
  geom_density_ridges(colour = "white", alpha = .7, from = 0, to = 1, scale = 1) +
  theme_ridges() +
  ylab("") +
  xlab("Coverage") +
  scale_fill_manual(values = c("#045a8d", "#b30000", "#006d2c", "#984ea3"), breaks = c("metagenomics", "metatranscriptomics", "metaproteomics", "metabolomics"))+
  ggtitle("Pathway classes species-independent") +
  theme(legend.position = "none")

metadata <- read_csv("../data/kegg_pathways/keggCLASS_nummet.tsv.gz", col_names = c("pathway", "nummet"), skip = 1) %>% 
  left_join(read_csv("../data/kegg_pathways/keggCLASS_num_reactions.tsv.gz", col_names = c("pathway", "nreactions"), skip = 1)) %>% 
  mutate(pathway = factor(pathway, levels = order)) %>% 
  filter(!is.na(pathway))

nummet <- ggplot(metadata, aes(x = pathway, y = nummet, label = nummet)) +
  geom_col(fill = "lightblue", width = 0.5) +
  coord_flip() + 
  geom_text(hjust = 1) +
  scale_y_log10()+ 
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.5))


nreactions <- ggplot(metadata, aes(x = pathway, y = nreactions, label = nreactions)) +
  geom_col(fill = "lightblue", width = 0.5) +
  coord_flip() + 
  geom_text(hjust = 1) +
  scale_y_log10() +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.5)) 


svg("Figure_2_D.svg")

plot_grid(coverage, nummet, nreactions, align = "h", axis = "ltb",  rel_widths = c(10, 1, 1), nrow =1)

dev.off()


