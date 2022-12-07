pacman::p_load(tidyverse)
pacman::p_load(patchwork)
# library(ggh4x)
# library(scales)

data <- tibble(
  omic = c("metaT", "metaP", "metaB")
) %>% 
  # Read in data and reformat
  mutate(data = map(omic, function(x){
    df <- read_csv(str_c("../data/kegg_pathways/keggPTWenrichment_fischer_", x, ".csv.gz")) %>% 
      select(PTWNAME,
             pvalue_chlor, pFDR_chlor, meanlog2FC_chlor, num_change_ptw_chlor, num_nochange_ptw_chlor,
             pvalue_metfo, pFDR_metfo, meanlog2FC_metfo, num_change_ptw_metfo, num_nochange_ptw_metfo,
             pvalue_niclo, pFDR_niclo, meanlog2FC_niclo, num_change_ptw_niclo, num_nochange_ptw_niclo) %>% 
      filter(!is.na(PTWNAME)) %>% 
      mutate(PTWNAME = str_replace(PTWNAME, "\\[\\'", ""),
             PTWNAME = str_replace(PTWNAME, "\\'\\]", ""))
    
    pvals <- df %>%
      select(PTWNAME, pvalue_chlor, pvalue_metfo, pvalue_niclo) %>% 
      pivot_longer(-PTWNAME, names_to = "drug", values_to = "pval") %>% 
      mutate(drug = str_extract(drug, "(?<=_).+")) 
    
    padj <- df %>%
      select(PTWNAME, pFDR_chlor, pFDR_metfo, pFDR_niclo) %>% 
      pivot_longer(-PTWNAME, names_to = "drug", values_to = "padj") %>% 
      mutate(drug = str_extract(drug, "(?<=_).+")) 
    
    log2fc <- df %>%
      select(PTWNAME, meanlog2FC_chlor, meanlog2FC_metfo, meanlog2FC_niclo) %>% 
      pivot_longer(-PTWNAME, names_to = "drug", values_to = "meanlog2fc") %>% 
      mutate(drug = str_extract(drug, "(?<=_).+")) 
    
    num_change_ptw <- df %>%
      select(PTWNAME, num_change_ptw_chlor, num_change_ptw_metfo, num_change_ptw_niclo) %>% 
      pivot_longer(-PTWNAME, names_to = "drug", values_to = "num_change_ptw") %>% 
      mutate(drug = str_extract(drug, "(?<=ptw_).+$")) 
    
    num_nochange_ptw <- df %>%
      select(PTWNAME, num_nochange_ptw_chlor, num_nochange_ptw_metfo, num_nochange_ptw_niclo) %>% 
      pivot_longer(-PTWNAME, names_to = "drug", values_to = "num_nochange_ptw") %>% 
      mutate(drug = str_extract(drug, "(?<=ptw_).+$")) 
    
    pvals %>% 
      left_join(padj) %>% 
      left_join(log2fc) %>% 
      left_join(num_change_ptw) %>% 
      left_join(num_nochange_ptw) %>% 
      mutate(num_total_ptw = num_change_ptw + num_nochange_ptw)
    
  }
  )
  )  %>% 
  unnest(data) %>% 
  rename(pathway = PTWNAME)

kegg_pathways_curated <- read_tsv("../data/kegg_pathways/kegg_pathways_filtered_by_community_members.csv.gz") %>% 
  mutate(ID = str_replace(ID, "path:", "")) %>% 
  # Filter out the Global and overview maps
  filter(! (str_detect(ID, "map011") | str_detect(ID, "map012"))) %>% 
  rename(KEGG_Pathway = ID) %>% 
  mutate(Name = str_replace(Name, "\t", ","))

data <- data %>% 
  filter(pathway %in% kegg_pathways_curated$Name)

data <- data %>% 
  filter(num_total_ptw > 3)

data_pval <- data %>%  
    mutate(
    significance = case_when(
      pval > 0.05 ~ "",
      pval <= 0.05 & pval > 0.01 ~ "*",
      pval <= 0.01 & pval > 0.001 ~ "**",
      pval <= 0.001 ~ "***"
    )
  )  

n_significant <- data_pval %>% 
  group_by(pathway, significance) %>% 
  summarise(count = n()) %>% 
  filter(!is.na(significance),
         !significance == "")

pathway_order <- n_significant %>% 
  group_by(pathway) %>% 
  summarise(total_count = sum(count)) %>% 
  mutate(pathway = factor(pathway) %>% fct_reorder(total_count) %>% fct_rev()) %>%
  pull(pathway)

n_per_drug <- data_pval %>%  
  filter(!is.na(significance),
         !significance == "") %>% 
  group_by(drug, pathway) %>%  
  summarise(count = n()) %>% 
  group_by(pathway) %>% 
  mutate(relabun = count/sum(count))

pdf("Figure_4_B.pdf", height=8, width=12)

data_pval %>% 
  complete(omic, drug, pathway) %>% 
  mutate(pathway = factor(pathway, levels = levels(pathway_order))) %>% 
  mutate(omic = factor(omic, levels = c("metaB", "metaP", "metaT"))) %>% 
  mutate(drug = factor(drug, levels = c("niclo", "metfo", "chlor"))) %>% 
  # Filter pathways that do not show any significance at all
  filter(!is.na(pathway)) %>% {
  ggplot(., aes(x = pathway, y = omic, fill = meanlog2fc)) +
  geom_tile(colour = "grey") +
  geom_text(aes(label = significance)) + 
  facet_wrap(~ drug, ncol = 1) + 
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.spacing = unit(0.2, "lines"),
        panel.background = element_blank()) +
  scale_fill_distiller(palette = "RdBu", 
                       limits = max(abs(.$meanlog2fc), na.rm = T) * c(-1, 1))
  }

dev.off()
