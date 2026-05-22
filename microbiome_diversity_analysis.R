############################################################
# Microbiome Diversity Analysis in R
# Author: Sufia Akter Neha
#
# Description:
# Clean and reproducible workflow for importing microbiome
# count data, taxonomy, metadata, and a phylogenetic tree;
# creating a phyloseq object; estimating alpha diversity;
# visualizing community composition; and testing beta diversity.
#
# Required input files:
#   zotutab.txt   = OTU/ZOTU count table; taxa as rows, samples as columns
#   zotus.tax     = taxonomy table from SINTAX or similar classifier
#   met.txt       = sample metadata; samples as rows
#   tree.tre      = phylogenetic tree
############################################################

############################
# 1. Load libraries
############################

library(phyloseq)
library(ape)
library(vegan)
library(tidyverse)
library(scales)

############################
# 2. User-defined input files
############################

otu_file  <- "zotutab.txt"
tax_file  <- "zotus.tax"
meta_file <- "met.txt"
tree_file <- "tree.tre"

output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)

############################
# 3. Import OTU/ZOTU table
############################

otu <- read.delim(otu_file, row.names = 1, check.names = FALSE)
otu <- as.data.frame(otu)
otu[] <- lapply(otu, as.numeric)

############################
# 4. Import and clean taxonomy
############################

tax_raw <- read.delim(tax_file, row.names = 1, header = FALSE, check.names = FALSE)

tax_string <- tax_raw[, 1]

tax_split <- strsplit(as.character(tax_string), ",")

tax <- do.call(rbind, lapply(tax_split, function(x) {
  length(x) <- 7
  x
}))

tax <- as.data.frame(tax)
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(tax) <- rownames(tax_raw)

tax <- tax %>%
  mutate(across(everything(), ~ gsub("^[a-zA-Z]:", "", .))) %>%
  mutate(across(everything(), ~ gsub("\\(.*\\)", "", .))) %>%
  mutate(across(everything(), ~ trimws(.)))

############################
# 5. Import metadata
############################

metadata <- read.delim(meta_file, row.names = 1, check.names = FALSE)
metadata$Sample_ID <- rownames(metadata)

############################
# 6. Match samples across OTU table and metadata
############################

shared_samples <- intersect(colnames(otu), rownames(metadata))

otu <- otu[, shared_samples, drop = FALSE]
metadata <- metadata[shared_samples, , drop = FALSE]

############################
# 7. Import phylogenetic tree
############################

tree <- read.tree(tree_file)

shared_taxa <- Reduce(intersect, list(rownames(otu), rownames(tax), tree$tip.label))

otu <- otu[shared_taxa, , drop = FALSE]
tax <- tax[shared_taxa, , drop = FALSE]
tree <- keep.tip(tree, shared_taxa)

############################
# 8. Create phyloseq object
############################

ps <- phyloseq(
  otu_table(as.matrix(otu), taxa_are_rows = TRUE),
  tax_table(as.matrix(tax)),
  sample_data(metadata),
  phy_tree(tree)
)

ps

############################
# 9. Rarefy samples
############################

set.seed(123)

min_depth <- min(sample_sums(ps))

ps_rare <- rarefy_even_depth(
  ps,
  sample.size = min_depth,
  rngseed = 123,
  replace = FALSE,
  verbose = FALSE
)

############################
# 10. Alpha diversity
############################

alpha_df <- estimate_richness(ps_rare, measures = c("Observed", "Shannon", "Simpson")) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(metadata %>% rownames_to_column("Sample_ID"), by = "Sample_ID")

write.csv(alpha_df, file.path(output_dir, "alpha_diversity.csv"), row.names = FALSE)

# Change this variable to match your metadata column
group_var <- "Sites"

p_alpha <- ggplot(alpha_df, aes(x = .data[[group_var]], y = Shannon, color = .data[[group_var]])) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.8, size = 2) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = group_var, y = "Shannon diversity")

ggsave(
  filename = file.path(output_dir, "alpha_diversity_shannon.png"),
  plot = p_alpha,
  width = 7,
  height = 5,
  dpi = 600
)

############################
# 11. Relative abundance plot
############################

ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

genus_df <- psmelt(ps_genus_rel)

genus_df$Genus <- as.character(genus_df$Genus)
genus_df$Genus[is.na(genus_df$Genus) | genus_df$Genus == ""] <- "Unclassified"

top_genera <- genus_df %>%
  group_by(Genus) %>%
  summarize(mean_abundance = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

genus_df <- genus_df %>%
  mutate(Genus_plot = ifelse(Genus %in% top_genera, Genus, "Other"))

genus_summary <- genus_df %>%
  group_by(Sample, Genus_plot, .data[[group_var]]) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

p_bar <- ggplot(genus_summary, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
  geom_bar(stat = "identity") +
  facet_wrap(as.formula(paste("~", group_var)), scales = "free_x") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = "Relative abundance", fill = "Genus")

ggsave(
  filename = file.path(output_dir, "genus_relative_abundance.png"),
  plot = p_bar,
  width = 12,
  height = 7,
  dpi = 600
)

############################
# 12. Beta diversity: Bray-Curtis PCoA
############################

bray_dist <- phyloseq::distance(ps_rare, method = "bray")

ordination <- ordinate(ps_rare, method = "PCoA", distance = bray_dist)

p_pcoa <- plot_ordination(ps_rare, ordination, color = group_var) +
  geom_point(size = 4, alpha = 0.8) +
  theme_bw(base_size = 14) +
  labs(color = group_var)

ggsave(
  filename = file.path(output_dir, "bray_curtis_pcoa.png"),
  plot = p_pcoa,
  width = 7,
  height = 5,
  dpi = 600
)

############################
# 13. PERMANOVA
############################

metadata_rare <- data.frame(sample_data(ps_rare))

adonis_formula <- as.formula(paste("bray_dist ~", group_var))

permanova_res <- adonis2(
  formula = adonis_formula,
  data = metadata_rare,
  permutations = 999
)

capture.output(
  permanova_res,
  file = file.path(output_dir, "permanova_bray_curtis.txt")
)

############################
# 14. Beta dispersion test
############################

dispersion <- betadisper(bray_dist, metadata_rare[[group_var]])
dispersion_res <- permutest(dispersion, permutations = 999)

capture.output(
  dispersion_res,
  file = file.path(output_dir, "beta_dispersion_test.txt")
)

