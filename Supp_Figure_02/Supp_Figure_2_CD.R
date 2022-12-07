pacman::p_load(tidyverse)
pacman::p_load(RColorBrewer)
pacman::p_load(vegan)
pacman::p_load(pheatmap)

mantel_cor_technical_cor <- function(name){
  res <- list()
  
  ## sample names shared in all datasets
  dup <- table(c(colnames(g.sp), colnames(t.sp), colnames(a.full), colnames(p.sp), colnames(b.both)))
  shared_sample <- names(dup[dup == 5])
  shared_sample
  
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
    }
    d <- d[, keep] %>% data.frame()
    rownames(d) <- rowname
    return(d)
  }
  g.sp <- ext_shared_sample(g.sp)
  g.gene <- ext_shared_sample(g.gene)
  t.sp <- ext_shared_sample(t.sp)
  t.gene <- ext_shared_sample(t.gene)
  p.sp <- ext_shared_sample(p.sp)
  p.gene <- ext_shared_sample(p.gene)
  a.full <- ext_shared_sample(a.full)
  b.both <- ext_shared_sample(b.both)
  
  ## transform data into relative abundance
  g.sp <- g.sp %>% as.matrix() %>% prop.table(margin = 2)
  g.gene <- g.gene %>% as.matrix() %>% prop.table(margin = 2)
  t.sp <- t.sp %>% as.matrix() %>% prop.table(margin = 2)
  t.gene <- t.gene %>% as.matrix() %>% prop.table(margin = 2)
  p.sp <- p.sp %>% as.matrix() %>% prop.table(margin = 2)
  p.gene <- p.gene %>% as.matrix() %>% prop.table(margin = 2)
  a.full <- a.full %>% as.matrix() %>% prop.table(margin = 2)
  b.both <- b.both %>% as.matrix() %>% prop.table(margin = 2)
  
  ## excludes features less than threshold
  filter_threshold <- function(d, num){
    keep <- apply(d, 1, mean) > threshold.vec[num]
    d <- d[keep, ]
    return(d)
  }
  g.sp <- filter_threshold(g.sp, 1)
  g.gene <- filter_threshold(g.gene, 2)
  t.sp <- filter_threshold(t.sp, 3)
  t.gene <- filter_threshold(t.gene, 4)
  p.sp <- filter_threshold(p.sp, 5)
  p.gene <- filter_threshold(p.gene, 6)
  
  # replace 0 with non-zero minimum value (metaB)
  for(i in 1:nrow(b.both)){
    min <- min(b.both[i, b.both[i, ] != 0])
    b.both[i, b.both[i, ] == 0] <- min
  }
  
  ## count number of features 
  num.vec <- c(nrow(g.sp), nrow(g.gene), nrow(t.sp), nrow(t.gene), nrow(p.sp), nrow(p.gene), nrow(a.full), nrow(b.both))
  
  ## modify name
  name <- paste0(name, " (n=", num.vec, ")")
  name <- gsub("\\) \\(", "\\, ", name)
  
  ## put all data into list
  dat <- list()
  dat[[1]] <- g.sp
  dat[[2]] <- g.gene
  dat[[3]] <- t.sp
  dat[[4]] <- t.gene
  dat[[5]] <- p.sp
  dat[[6]] <- p.gene
  dat[[7]] <- a.full
  dat[[8]] <- b.both
  
  ## check size and NAs of the datasets
  dim(dat[[1]])
  dim(dat[[2]])
  dim(dat[[3]])
  dim(dat[[4]])
  dim(dat[[5]])
  dim(dat[[6]])
  dim(dat[[7]])
  dim(dat[[8]])
  
  sum(is.na(dat[[1]]))
  sum(is.na(dat[[2]]))
  sum(is.na(dat[[3]]))
  sum(is.na(dat[[4]]))
  sum(is.na(dat[[5]]))
  sum(is.na(dat[[6]]))
  sum(is.na(dat[[7]]))
  sum(is.na(dat[[8]]))
  
  ## sort rows and columns
  dat <- map(dat, function(df){
    df <- df[, order(colnames(df))]
    #df <- df[order(rownames(df)), ]  
    return(df)
  }
  )
  
  ## calc distance
  cor.spe <- map(dat, function(df){1 - cor(df, method = "spearman")})
  cor.pea <- map(dat, function(df){1 - cor(df, method = "pearson")})
  
  ## plot mantel correlation using pheatmap
  breaks <- seq(0, 1, by = 0.01)
  col <- colorRampPalette(c("white", "#ffffbf", "firebrick3"))(length(breaks))
  
  ## calc correlation between technical replicates
  ext.list <- c(5, 16, 21, 22, 35, 46, 59) # ex samples without corresponding replicate
  pair.spe <- list()
  pair.pea <- list()  
  pair.spe <- map(cor.spe, function(df){df[-ext.list, -ext.list]})
  pair.pea <- map(cor.pea, function(df){df[-ext.list, -ext.list]})  
  
  seq1 <- seq(1, ncol(pair.spe[[1]])-1, by = 2) # replicate1
  seq2 <- seq(2, ncol(pair.spe[[1]]), by = 2) # replicate2
  
  diag.spe <- map(pair.spe, function(df){diag(df[seq1, seq2])})
  diag.pea <- map(pair.pea, function(df){diag(df[seq1, seq2])})  
  
  # spearman correlation
  df <- data.frame(diag.spe)
  colnames(df) <- name
  df$pair1 <- colnames(pair.spe[[1]])[seq1]
  df$pair2 <- colnames(pair.spe[[1]])[seq2]
  
  df <- gather(df, 1:(ncol(df)-2), key = "data", value = "var")
  df$data <- factor(df$data, levels = name)
  df$run <- ifelse(str_detect(df$pair1, "^A") & str_detect(df$pair2, "^A"), "A", "B")
  
  plot.spe <- ggplot(df, aes(x = data, y = 1 - var, fill = run)) +
    geom_boxplot(show.legend = F) +
    ylab("Spearman correlation") +
    xlab("") +
    ylim(c(0, 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ## pearson correlation
  df <- data.frame(diag.pea)
  colnames(df) <- name
  df$pair1 <- colnames(pair.pea[[1]])[seq1]
  df$pair2 <- colnames(pair.pea[[1]])[seq2]
  df <- gather(df, 1:(ncol(df)-2), key = "data", value = "var")
  df$data <- factor(df$data, levels = name)
  df$run <- ifelse(str_detect(df$pair1, "^A") & str_detect(df$pair2, "^A"), "A", "B")

  plot.pea <- ggplot(df, aes(x = data, y = 1 - var, fill = run)) +
    geom_boxplot(show.legend = F) +
    ylab("Pearson correlation") +
    xlab("") +
    ylim(c(0, 1)) +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  res[[1]] <- plot.pea
  res[[2]] <- plot.spe
  return(res)
}

g.sp <- read_tsv("../data/relative_abundance_tables/metagenomics_relabun.tsv.gz")
g.gene <- read_tsv("../data/gene_count_tables/metagenomics_gene_counts_raw.tsv.gz")
t.sp <- read_tsv("../data/relative_abundance_tables/metatranscriptomics_relabun.tsv.gz")
t.gene <- read_tsv("../data/gene_count_tables/metatranscriptomics_gene_counts_raw.tsv.gz")
p.sp <- read_tsv("../data/relative_abundance_tables/metaproteomics_relabun.tsv.gz")
p.gene <- read_tsv("../data/gene_count_tables/metaproteomics_gene_counts_raw.tsv.gz")
a.full <- read_tsv("../data/relative_abundance_tables/amplicon_sequencing_relabun.tsv.gz")

b.both <- read_tsv("../data/metabolomics/metabolomics_jointdata_quantnorm_ann_filtered.tsv.gz") %>%
  select(-contains("nobact")) %>%
  select(-contains("aqmode")) %>%
  select(-contains("Ion_index"))

threshold.vec <- c(0, 1E-7, 0, 1E-7, 0, 1E-5, 0, 0)
name <- c("metaG (species)", "metaG (genes)", "metaT (species)", "metaT (genes)", "metaP (species)", "metaP (genes)", "16S rRNA profiling", "metaB")

res <- mantel_cor_technical_cor(name)

pdf("Supp_Figure_2_CD.pdf", width=6, height=4)

res[[1]]
res[[2]]

dev.off()
