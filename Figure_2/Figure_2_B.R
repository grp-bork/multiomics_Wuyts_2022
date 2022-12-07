pacman::p_load(tidyverse)
pacman::p_load(ggthemes)
pacman::p_load(ggrastr)

prepare_df <- function(g, t, p, pseudvalue){
  
  ## sample names shared in all datasets
  dup <- table(c(colnames(g), colnames(t), colnames(p)))
  shared_sample <- names(dup[dup == 3])
  
  ## use only the shared samples
  ext_shared_sample <- function(d){
    keep <- which(colnames(d) %in% shared_sample)
    if(sum(grep("gene_name", colnames(d))) > 0){
      rowname <- d$gene_name
      print("gene_name")
    }else if(sum(grep("species", colnames(d))) > 0){
      rowname <- d$species
      print("species")
    }else if(sum(grep("m/z", colnames(d))) > 0){
      rowname <- d$`m/z`
      print("m/z")      
    }else if(sum(grep("X1", colnames(d))) > 0){
      rowname <- d$X1
      print("KO")      
    }
    d <- d[, keep] %>% data.frame()
    rownames(d) <- rowname
    return(d)
  }
  g <- ext_shared_sample(g)
  t <- ext_shared_sample(t)
  p <- ext_shared_sample(p)

  g <- g[, c(-1)]
  t <- t[, c(-1)]
  p <- p[, c(-1)]
  
  sort.order <- function(df){
    df <- df[, order(colnames(df))]
    df <- df[order(rownames(df)), ]  
    return(df)
  }
  g <- sort.order(g)
  t <- sort.order(t)
  p <- sort.order(p)
  data.frame(colnames(g), colnames(t), colnames(p))
  
  ## transform data into relative abundance
  g <- prop.table(as.matrix(g), margin = 2)
  t <- prop.table(as.matrix(t), margin = 2)
  p <- prop.table(as.matrix(p), margin = 2)
  
  ## calc mean
  g.mean <- apply(g, 1, mean)
  t.mean <- apply(t, 1, mean)
  p.mean <- apply(p, 1, mean)  

  ## add 1E-8 as pseudvalue
  g.mean <- g.mean + pseudvalue
  t.mean <- t.mean + pseudvalue
  p.mean <- p.mean + pseudvalue
  
  ## prepare df
  df <- data.frame(gene = names(g.mean), metagenomics = log10(g.mean), metatranscriptomics = log10(t.mean), metaproteomics = log10(p.mean))
  return(df)
}
## read data
g.ko <- read_tsv("../data/gene_count_tables/metagenomics_gene_counts_normalised_ko.tsv.gz") %>% select(-contains("Transfer")) %>% select(-contains("Inoculum"))
t.ko <- read_tsv("../data/gene_count_tables/metatranscriptomics_gene_counts_normalised_ko.tsv.gz") %>% select(-contains("Transfer"))
p.ko <- read_tsv("../data/gene_count_tables/metaproteomics_gene_counts_normalised_ko.tsv.gz") %>% select(-contains("n_genes_PG"))
df <- prepare_df(g.ko, t.ko, p.ko, 1E-8)

cor <- c()
cor[1] <- cor(df$metagenomics, df$metatranscriptomics, method = "spearman")
cor[2] <- cor(df$metagenomics, df$metaproteomics, method = "spearman")
cor[3] <- cor(df$metatranscriptomics, df$metaproteomics, method = "spearman")
cor <- paste0("rho = ", round(cor, digits = 2))

pdf("Figure_2_B.pdf", width=3, height=3)

## metaG vs metaT
ggplot(df, aes(x = metagenomics, y = metatranscriptomics)) +
  geom_abline(intercept = 0, slope = 1, col = "gray80") +
  theme_few() +
  geom_point_rast(size = 0.1) +  
  annotate(geom = 'text', x = -3.5, y = -7, label = cor[1], hjust = 0, vjust = 1.5) +
  scale_y_continuous(limits = c(-8, -1.5)) +
  scale_x_continuous(limits = c(-8, -1.5)) 

## metaG vs metaP
ggplot(df, aes(x = metagenomics, y = metaproteomics)) +
  geom_abline(intercept = 0, slope = 1, col = "gray80") +
  theme_few() +
  geom_point_rast(size = 0.1) +  
  annotate(geom = 'text', x = -3.5, y = -7, label = cor[2], hjust = 0, vjust = 1.5) +  
  scale_y_continuous(limits = c(-8, -1.5)) +
  scale_x_continuous(limits = c(-8, -1.5)) 

## metaT vs metaP
ggplot(df, aes(x = metatranscriptomics, y = metaproteomics)) +
  geom_abline(intercept = 0, slope = 1, col = "gray80") +
  theme_few() +
  geom_point_rast(size = 0.1) +  
  annotate(geom = 'text', x = -3.5, y = -7, label = cor[3], hjust = 0, vjust = 1.5) +  
  scale_y_continuous(limits = c(-8, -1.5)) +
  scale_x_continuous(limits = c(-8, -1.5)) 

dev.off()
