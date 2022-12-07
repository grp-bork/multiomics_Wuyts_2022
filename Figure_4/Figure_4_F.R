pacman::p_load(tidyverse)
pacman::p_load(patchwork)

# Counts
metaB_count <- read_tsv("../data/metabolomics/metabolomics_jointdata_quantnorm_ann_filtered.tsv.gz")%>% 
  select(-`m/z`, -aqmode) %>% 
  pivot_longer(-`Ion_index`, names_to = "sample", values_to = "count") %>% 
  separate(sample, into = c("run", "drug", "timepoint", "repeat"), extra = "merge") %>% 
  mutate(type = if_else(str_detect(`repeat`, "nobact"), "nobact", "bact")) %>% 
  mutate(timepoint = factor(timepoint, levels = c("T0", "T15", "T30", "T1h", "T3h", "T48h", "T96h")))

# Get all metabolites that are differential abundant
## get p-value from over all timepoints
metaB_differential_abundant_pvalue <- read_tsv("../data/ranova/ranova.separate.metaB.log10.runAB.tsv.gz") %>% 
  filter(ihw.fdr <= 0.05) %>% 
  select(name, drug, ihw.fdr)

## Get log2foldchange from per timepoint analysis
## Only use the first few timepoints
metaB_differential_abundant_log2fold <- read_tsv("../data/ranova/ranova.metaB.runAB.post_hoc_ttest.tsv.gz") %>% 
  select(-...1) %>% 
  rename(drug = `drug...4`) %>% 
  filter(time %in% c("T15","T30", "T1h", "T3h")) %>% 
  group_by(`m/z`, drug) %>% 
  summarise(max_log2fc = max(abs(log2fc), na.rm = T)) %>% 
  filter(max_log2fc >= log2(1.5)) %>% 
  rename(name = `m/z`) %>% 
  select(name, drug, max_log2fc)

metaB_differential_abundant <- inner_join(metaB_differential_abundant_pvalue, metaB_differential_abundant_log2fold)

# Curated list of pathways
kegg_pathways_curated <- read_tsv("../data/kegg_pathways/kegg_pathways_filtered_by_community_members.csv.gz") %>% 
  mutate(ID = str_replace(ID, "path:", "")) %>% 
  # Filter out the Gloval and overview maps
  filter(! (str_detect(ID, "map011") | str_detect(ID, "map012"))) %>% 
  rename(KEGG_Pathway = ID)

# Link KEGG compound ID to KEGG pathway
compound_pathways <- read_tsv("../data/kegg_pathways/kegg_pathway_compound.tsv.gz", col_names = c("KEGG_Pathway", "KEGG_ID")) %>% 
  mutate(KEGG_Pathway = str_replace(KEGG_Pathway, "path:", ""))

# Read in complete annotation and add pathway data
metaB_annotation <- read_tsv("../data/metabolomics/metabolomics_joint_ann_filtered.tsv.gz")%>% 
  select(Ion_index, Ion_mz, KEGG_ID, Metabolite_name_HMDB, SMILES) %>% 
  left_join(compound_pathways) %>% 
  left_join(kegg_pathways_curated)%>% 
  filter(!is.na(Name))

metabolites_of_interest <- metaB_annotation %>% 
  filter(Name %in% c(
    "Galactose metabolism"
  )) %>% 
  filter(Ion_mz %in% metaB_differential_abundant$name)

concat_names <- metabolites_of_interest %>% 
  group_by(Ion_mz) %>% 
  summarise(names = str_c(unique(Metabolite_name_HMDB), collapse = "; ")) 

metabolite_plots <- metaB_count %>% 
  filter(Ion_index %in% metabolites_of_interest$Ion_index) %>% 
  left_join(metabolites_of_interest) %>%
  left_join(concat_names) %>% 
  mutate(count = if_else(count == 0, NaN, count)) %>% 
  group_by(Name) %>% 
  nest() %>% 
  mutate(plot = map2(data, Name,
                     ~.x %>% 
                       ggplot(aes(x = timepoint, y = count, colour = drug, group = interaction(run, Ion_index, drug, `repeat`))) +
                       geom_line(aes(linetype = type)) +
                       facet_wrap(~names, scales = "free_y", nrow = 1) +
                       ggtitle(.y) +
                       scale_colour_brewer(palette = "Dark2")
                    )
  ) 

metabolites_of_interest_chlor <- metaB_annotation %>% 
  filter(Name %in% c(
    "Galactose metabolism",
    "Fructose and mannose metabolism",
    "Starch and sucrose metabolism"
  )) %>% 
  filter(Ion_mz %in% metaB_differential_abundant$name)

#Select di-, tri- and tetrasaccharides for plotting.  
metaB_selected_norm <- metaB_count %>% 
  filter(Ion_index %in% metabolites_of_interest$Ion_index) %>% 
  left_join(metabolites_of_interest_chlor) %>%
  left_join(concat_names) %>% 
  select(-KEGG_Pathway, -Name) %>% 
  distinct() %>% 
  mutate(count = if_else(count == 0, NaN, count)) %>%
  mutate(names = str_replace_all(names, ";", "\n"))


#Normalize each time course in each sample by T0.
for (run in unique(metaB_selected_norm$run)){
  for (drug in unique(metaB_selected_norm$drug)){
    for (rep in unique(metaB_selected_norm$`repeat`))
      for (type in unique(metaB_selected_norm$type)){
        for (iidx in unique(metaB_selected_norm$Ion_index)){
          temp_select <- metaB_selected_norm[(metaB_selected_norm$run == run) &
               (metaB_selected_norm$drug == drug) &
               (metaB_selected_norm$`repeat` == rep) &
               (metaB_selected_norm$type == type) &
               (metaB_selected_norm$Ion_index) == iidx,]
          metaB_selected_norm[(metaB_selected_norm$run == run) &
               (metaB_selected_norm$drug == drug) &
               (metaB_selected_norm$`repeat` == rep) &
               (metaB_selected_norm$type == type) &
               (metaB_selected_norm$Ion_index) == iidx,]$count <- metaB_selected_norm[(metaB_selected_norm$run == run) &
               (metaB_selected_norm$drug == drug) &
               (metaB_selected_norm$`repeat` == rep) &
               (metaB_selected_norm$type == type) &
               (metaB_selected_norm$Ion_index) == iidx,]$count / 
            mean(temp_select[temp_select$timepoint == 'T0', ]$count)
        }
      }
  }
}

#Plot sugars normalized by T0.
pdf("Figure_4_F.pdf", width=12, height=8)

metaB_selected_norm %>% 
  ggplot(aes(x = timepoint, y = count, colour = drug, group = interaction(drug, run, `repeat`, type))) +
  geom_line(aes(linetype = type)) +
  facet_wrap(~names, nrow = 1) + #scales = "free_y", 
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 0)) +
  #ylim(0, 1) + 
  ggtitle("Chlorpromazine selected sugar metabolites normalized to T0") 

dev.off()
