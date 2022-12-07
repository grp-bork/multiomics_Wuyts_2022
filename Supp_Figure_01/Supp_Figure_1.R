pacman::p_load(tidyverse)
pacman::p_load(vegan)
pacman::p_load(ggrepel)
pacman::p_load(stringi)

# Read in abundance
abundance <- tibble(
  omic = c("amplicon_sequencing", "metagenomics", "metatranscriptomics")
)  %>% 
  mutate(abundance = map(omic, ~read_tsv(str_c("../data/relative_abundance_tables/", .x, "_relabun.tsv.gz")) %>% 
                           select(accession_number, species_name, contains("Inoculum"), contains("Transfer")) %>% 
                           pivot_longer(c(-accession_number, -species_name), names_to = "sample", values_to = "abundance") %>% 
                           separate(sample, into = c("run", "time", "tech_rep"), remove = F)
                           )
         ) %>% 
  unnest(abundance)


mean_abundance <- abundance %>% 
  filter(omic == "amplicon_sequencing") %>% 
  group_by(run, time, species_name) %>% 
  summarise(mean_abundance = mean(abundance))

# Get distance matrix
distmat <- mean_abundance %>% 
  ungroup() %>% 
  mutate(sample = str_c(run, time, sep = "_")) %>% 
  select(-run, -time) %>% 
  pivot_wider(names_from = sample, values_from = mean_abundance) %>%
  column_to_rownames("species_name") %>% 
  # Calculate Bray-curtis distance
  as.matrix() %>% 
  t() %>% 
  vegdist(method = "bray", upper = T, diag = T)

pdf("Supp_Figure_1_C.pdf")

# Sample order according to time
sample_order <- mean_abundance %>% 
  mutate(sample = str_c(run, time, sep = "_")) %>% 
  select(sample, time, run) %>% 
  distinct() %>% 
  arrange(time) %>% 
  pull(sample)


distmat %>%
  as.matrix() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "sample1")  %>% 
  pivot_longer(-sample1, names_to = "sample2", values_to = "dist")%>% 
  separate(sample1, into = c("run1", "time1"), remove = F) %>%
  separate(sample2, into = c("run2", "time2"), remove = F) %>% 
  mutate(sample1 = factor(sample1, levels = sample_order),
         sample2 = factor(sample2, levels = sample_order)) %>% 
  # Calculate similarity
  ggplot(aes(x = sample1, y = sample2, fill = dist, label = round(dist, 2))) +
  geom_tile(colour = "lightgrey", size = 0.5) +
  geom_text(size = 3, alpha = 0.8) + 
  expand_limits(fill = 1) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  labs(
    x = "",
    y = "",
    fill = "Bray-Curtis dissimilarity index"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    
    # Gridlines
    panel.grid = element_blank(),
    
    # Backgrounds
    panel.background = element_blank(),
    plot.background = element_blank(),
    
    # Facet titles
    strip.text.x = element_text(face = "bold", size = 14),
    strip.background = element_blank()
    
    )

dev.off()

pdf("Supp_Figure_1_B.pdf")

invsimpson <- abundance %>% 
  filter(omic == "amplicon_sequencing") %>% 
  select(-run, -time, -tech_rep, -omic, - accession_number) %>% 
  pivot_wider(names_from = sample, values_from = abundance) %>%
  column_to_rownames("species_name") %>% 
  # Calculate Bray-curtis distance
  as.matrix() %>% 
  t() %>% 
  diversity(index = "invsimpson")

tibble(
  sample = names(invsimpson),
  invsimpson = invsimpson
) %>% 
  separate(sample, into = c("run", "time", "tech_rep"), remove = T) %>% 
  ggplot(aes(x = time, y = invsimpson, group = run, colour = run)) +
  stat_summary(fun = "mean", geom = "line") +
  geom_point() +
  expand_limits(y = 0) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Alpha diversity",
    y = "Inverse simpson index",
    x = "",
    colour = "Run"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank()
  )

dev.off()

# Find top 12 species
species_order <- mean_abundance %>%
    group_by(species_name) %>%
    summarise(sum_relabun = sum(mean_abundance)) %>%
    arrange(-sum_relabun) %>%
    pull(species_name)
top12 <- species_order[1:12]

pdf("Supp_Figure_1_A.pdf")

mean_abundance %>% 
  mutate(species_name = if_else(species_name %in% top12, species_name, "other"), 
         species_name = factor(species_name, levels = c(top12, "other"))
  ) %>% 
  ggplot(aes(x = time, y = mean_abundance, fill = species_name)) +
  geom_col() +
  ## Make a scale of 13 colours based on a colourbrewer palette
  scale_fill_manual(values = 
                      c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                        '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                        '#cab2d6','#6a3d9a','#ffff99','#b15928',
                        "lightgrey")) +
  facet_wrap(~run, nrow = 1, scales = "free_x") +
  labs(
    y = "mean relative abundance",
    x = ""
  ) +
  ## Customise lay-out
  theme_minimal() +
  theme(
    # Legend
    legend.position = "right",
    legend.title = element_blank(),
    
    # Axes labeling
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    
    # Gridlines
    panel.grid = element_blank(),
    
    # Backgrounds
    panel.background = element_blank(),
    plot.background = element_blank(),
    
    # Facet titles
    strip.text.x = element_text(face = "bold", size = 14),
    
    # Axis ticks
    axis.ticks.y = element_line()
  )

dev.off()

d <- read_tsv("../data/growth_curve/values.txt.gz", col_types = cols(
  Day = col_character(),
  time = col_double(),
  Sample = col_character(),
  OD = col_double()
))

d <- d %>% mutate(
  kind = stri_sub(Sample, 1, -3),
  measurement = stri_sub(Sample, -1, -1)
)

X <- seq(0, 12, 0.1)

fitSigmoid <- function(x) {
  for (timescale_start in 1:5) {
    try( {
      n <- nls( OD ~ lo + (hi-lo)/(1+exp(-(time-delta)*timescale ) ), data = x,
                start = list( lo = min(x$OD), hi = max(x$OD), delta = median(x$time), timescale = timescale_start ) )
      
      return(
        tibble(
          time = X,
          OD = predict(n, data.frame(time = X))
        )
      )
      return(as.data.frame(t(coef(n))))
    }, silent = T
    )
  }
  
  data.frame( lo = NA, hi = NA, delta = NA, timescale = NA)
}

md <- d %>% filter(kind != "medium_con") %>% group_by(Day, time) %>%
  summarise( OD = median(OD, na.rm=T) ) %>% mutate(kind = "measured")

fd <- md %>% group_by(Day) %>% do(fitSigmoid(.)) %>% mutate(kind = "fit")

dd <- bind_rows(md, fd) %>% group_by(Day) %>% 
  mutate(ODc0 = OD - min(OD[ kind == "fit" ]),
         ODc01 = ODc0 / max(ODc0[ kind == "fit" ])
  )

dd <- dd %>% mutate(day = as.factor(as.integer(stri_extract_first_regex(Day, "[1-9][0-9]*"))))

pdf("Supp_Figure_1_D.pdf", width=6, height=6)

dd %>% filter( kind == "fit", Day %in% c("SE07", "SE09", "SE11") ) %>%
  ggplot(aes(time, ODc01, color=day)) +
  annotate("rect", xmin = 5, xmax = 8, ymin = -Inf, ymax = Inf, alpha = .2, size = 0) +
  geom_vline(xintercept = 5) +
  geom_line() +
  geom_point(data = dd %>% filter( kind == "measured", Day %in% c("SE07", "SE09", "SE11") )) +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = seq(0,12,4)) +
  scale_y_continuous(name = "normalised OD") + 
  theme_minimal() +
  theme(legend.position = "bottom")

dev.off()
