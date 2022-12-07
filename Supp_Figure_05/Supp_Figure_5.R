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
  dim(g)
  dim(t)
  dim(p)
  
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
  
  length(g.mean)
  length(t.mean)
  length(p.mean)
  
  sum(g.mean == 0)
  sum(t.mean == 0)
  sum(p.mean == 0)
  
  ## add 1E-8 as pseudvalue
  g.mean <- g.mean + pseudvalue
  t.mean <- t.mean + pseudvalue
  p.mean <- p.mean + pseudvalue
  
  ## prepare df
  df <- data.frame(gene = names(g.mean), metagenomics = log10(g.mean), metatranscriptomics = log10(t.mean), metaproteomics = log10(p.mean))
  return(df)
}

## read data (gene)
g.gene <- read_tsv("../data/gene_count_tables/metagenomics_gene_counts_raw.tsv.gz") %>% select(-contains("species_name"))
t.gene <- read_tsv("../data/gene_count_tables/metatranscriptomics_gene_counts_raw.tsv.gz") %>% select(-contains("species_name"))
p.gene <- read_tsv("../data/gene_count_tables/metaproteomics_gene_counts_raw.tsv.gz") %>% select(-contains("species_name"))
df <- prepare_df(g.gene, t.gene, p.gene, 1E-8)

## correlation using all genes
cor <- c()
cor[1] <- cor(df$metagenomics, df$metatranscriptomics, method = "spearman")
cor[2] <- cor(df$metagenomics, df$metaproteomics, method = "spearman")
cor[3] <- cor(df$metatranscriptomics, df$metaproteomics, method = "spearman")

mylab <- c()
mylab[1] <- parse(text = sprintf('rho=="%.2f (all genes)"', cor[1]))
mylab[2] <- c(parse(text = sprintf('rho=="%.2f (all genes)"', cor[2])))
mylab[3] <- c(parse(text = sprintf('rho=="%.2f (all genes)"', cor[3])))

## correlatin using non-zero genes
cor.non_zero <- c()
df.gt <- df %>% filter(metagenomics != -8 & metatranscriptomics != -8)
df.gp <- df %>% filter(metagenomics != -8 & metaproteomics != -8)
df.tp <- df %>% filter(metatranscriptomics != -8 & metaproteomics != -8)

cor.non_zero[1] <- cor(df.gt$metagenomics, df.gt$metatranscriptomics, method = "spearman")
cor.non_zero[2] <- cor(df.gp$metagenomics, df.gp$metaproteomics, method = "spearman")
cor.non_zero[3] <- cor(df.tp$metatranscriptomics, df.tp$metaproteomics, method = "spearman")

cor(df.gt$metagenomics, df.gt$metatranscriptomics)
cor(df.gp$metagenomics, df.gp$metaproteomics)
cor(df.tp$metatranscriptomics, df.tp$metaproteomics)
cor.non_zero
#cor.non_zero <- paste0("rho = ", round(cor.non_zero, digits = 2), " (non-zero genes)")

mylab2 <- c()
mylab2[1] <- parse(text = sprintf('rho=="%.2f (non-zero genes)"', cor.non_zero[1]))
mylab2[2] <- c(parse(text = sprintf('rho=="%.2f (non-zero genes)"', cor.non_zero[2])))
mylab2[3] <- c(parse(text = sprintf('rho=="%.2f (non-zero genes)"', cor.non_zero[3])))

pdf("Supp_Figure_5.pdf", width=3, height=3)

## plot (gene)
ggplot(df, aes(x = metagenomics, y = metatranscriptomics)) +
  geom_abline(intercept = 0, slope = 1, col = "gray80") +
  theme_few() +
  geom_point_rast(size = 0.1, alpha = 0.1) +  
  annotate(geom = 'text', x = -8, y = -2, label = mylab[1], hjust = 0, vjust = 1.5, size = 2.5) +
  annotate(geom = 'text', x = -8, y = -2.3, label = mylab2[1], hjust = 0, vjust = 1.5, size = 2.5) +  
  scale_y_continuous(limits = c(-8, -2)) +
  scale_x_continuous(limits = c(-8, -2))

ggplot(df, aes(x = metagenomics, y = metaproteomics)) +
  geom_abline(intercept = 0, slope = 1, col = "gray80") +
  theme_few() +
  geom_point_rast(size = 0.1, alpha = 0.1) +  
  annotate(geom = 'text', x = -8, y = -2, label = mylab[2], hjust = 0, vjust = 1.5, size = 2.5) +
  annotate(geom = 'text', x = -8, y = -2.3, label = mylab2[2], hjust = 0, vjust = 1.5, size = 2.5) +  
  scale_y_continuous(limits = c(-8, -2)) +
  scale_x_continuous(limits = c(-8, -2))

ggplot(df, aes(x = metatranscriptomics, y = metaproteomics)) +
  geom_abline(intercept = 0, slope = 1, col = "gray80") +
  theme_few() +
  geom_point_rast(size = 0.1, alpha = 0.1) +  
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour = "red") +    
  annotate(geom = 'text', x = -8, y = -2, label = mylab[3], hjust = 0, vjust = 1.5, size = 2.5) +
  annotate(geom = 'text', x = -8, y = -2.3, label = mylab2[3], hjust = 0, vjust = 1.5, size = 2.5) +    
  scale_y_continuous(limits = c(-8, -2)) +
  scale_x_continuous(limits = c(-8, -2))

dev.off()
