pacman::p_load(tidyverse)

## read source file
d <- read.delim("../data/seq_stats/seq_stat.tsv.gz", header = T)

keep.g <- grepl("metaG", d$file)
keep.t <- grepl("metaT", d$file)

g <- d[keep.g, ]
t <- d[keep.t, ]

name.list <- unique(gsub("-metaG.*", "", g$file))
num_seqs <- c()
sum_len <- c()

## merge forward and reverse seq files (metaG)
for(i in seq(1, nrow(g), by = 2)){
  ii <- (i + 1)
  num_seqs.temp <- g$num_seqs[i] + g$num_seqs[ii]
  sum_len.temp <- g$sum_len[i] + g$sum_len[ii]
  
  num_seqs <- c(num_seqs, num_seqs.temp)
  sum_len <- c(sum_len, sum_len.temp)  
}
g <- data.frame(file = name.list, num_seqs = num_seqs, sum_len = sum_len, data = "Metagenomics")

t <- t[, c(1, 4, 5)]
t <- data.frame(t, data = "Metatranscriptomics")
df <- rbind(g, t)

pdf("Supp_Figure_2_AB.pdf", width=3, height=4)

ggplot(df, aes(x = data, y = sum_len)) +
  theme_classic() +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(width = 0.1, size = 0.1) +
  ylab("Total bases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title.x = element_blank())

ggplot(df, aes(x = data, y = num_seqs)) +
  theme_classic() +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(width = 0.1, size = 0.1) +
  ylab("Total reads") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title.x = element_blank())

dev.off()
