############################################################
## RNA-seq DESeq2 analysis — publication-ready (fixed)
## Outputs all plots to Plots_Generated/
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(org.At.tair.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(pheatmap)
  library(ggplot2)
})

theme_set(theme_bw(base_size = 14))

############################
## Create output dirs
############################
dir.create("Plots_Generated", showWarnings = FALSE)
dir.create("Plots_Generated/PCA", showWarnings = FALSE)
dir.create("Plots_Generated/heatmaps", showWarnings = FALSE)
dir.create("Plots_Generated/boxplots", showWarnings = FALSE)
dir.create("Plots_Generated/GO", showWarnings = FALSE)

############################
## 1. Metadata
############################
meta <- read.table(
  "tabularmeta.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  mutate(
    sample = factor(sample),
    treatment = factor(treatment, levels = c("mock", "DC3000")),
    temperature = factor(temperature, levels = c("23", "30"))
  )

############################
## 2. featureCounts files -> count matrix
############################
files <- list.files(
  "featureCounts_collection",
  pattern = "\\.tabular$",
  full.names = TRUE
)

sample_names <- gsub("\\.tabular$", "", basename(files))
names(files) <- sample_names

read_fc <- function(fpath) {
  df <- read.delim(fpath, header = TRUE, stringsAsFactors = FALSE)
  # Expect two columns: Geneid, count
  if (ncol(df) < 2) stop("featureCounts file has <2 columns: ", fpath)
  # second column might be named differently; keep as counts col
  counts_col <- colnames(df)[2]
  out <- df[, c("Geneid", counts_col)]
  colnames(out) <- c("Geneid", basename(fpath))
  out
}

count_list <- lapply(files, read_fc)

count_matrix_df <- Reduce(function(x, y) full_join(x, y, by = "Geneid"), count_list)
rownames(count_matrix_df) <- count_matrix_df$Geneid
count_matrix <- as.matrix(count_matrix_df[, -1])
count_matrix[is.na(count_matrix)] <- 0
colnames(count_matrix) <- sample_names

stopifnot(all(colnames(count_matrix) == meta$sample))

############################
## 3. DESeq object (no interaction)
############################
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = meta,
  design = ~ temperature + treatment
)

dds <- dds[rowSums(counts(dds)) > 10, ]

############################
## 4. Run contrasts (re-fit per subset)
############################
run_contrast <- function(dds, subset_idx, design_formula, contrast) {
  d <- dds[, subset_idx]
  design(d) <- design_formula
  d <- DESeq(d)
  results(d, contrast = contrast)
}

res_list <- list(
  DC3000_vs_mock_23C =
    run_contrast(dds, meta$temperature == "23",
                 ~ treatment,
                 c("treatment", "DC3000", "mock")),
  DC3000_vs_mock_30C =
    run_contrast(dds, meta$temperature == "30",
                 ~ treatment,
                 c("treatment", "DC3000", "mock")),
  Temp_30_vs_23_mock =
    run_contrast(dds, meta$treatment == "mock",
                 ~ temperature,
                 c("temperature", "30", "23")),
  Temp_30_vs_23_DC3000 =
    run_contrast(dds, meta$treatment == "DC3000",
                 ~ temperature,
                 c("temperature", "30", "23"))
)

# human-readable titles for outputs (keeps filenames short but titles informative)
contrast_titles <- c(
  DC3000_vs_mock_23C = "Infection (DC3000 vs mock) at 23°C",
  DC3000_vs_mock_30C = "Infection (DC3000 vs mock) at 30°C",
  Temp_30_vs_23_mock = "Temperature (30°C vs 23°C) in mock",
  Temp_30_vs_23_DC3000 = "Temperature (30°C vs 23°C) in DC3000"
)

############################
## 5. VST for visualization
############################
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

############################################################
## PCA (all samples)
############################################################
pca_df <- plotPCA(vsd, intgroup = c("treatment", "temperature"), returnData = TRUE)
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = treatment, shape = temperature)) +
  geom_point(size = 4) +
  labs(title = "PCA — all samples")

ggsave("Plots_Generated/PCA/PCA_all_samples.pdf", p_pca, width = 6, height = 5)
ggsave("Plots_Generated/PCA/PCA_all_samples.png", p_pca, dpi = 300, width = 6, height = 5)

############################################################
## Heatmaps for all contrasts (top 50 DE by padj)
############################################################
save_pheatmap_png <- function(mat, annotation_col, filename_png, main = NULL) {
  # pheatmap with png wrapper
  png(filename_png, width = 1600, height = 1200, res = 150)
  pheatmap(mat, annotation_col = annotation_col, show_rownames = FALSE, main = main)
  dev.off()
}

for (nm in names(res_list)) {
  res <- res_list[[nm]]
  samples <- switch(nm,
                    DC3000_vs_mock_23C = meta$sample[meta$temperature == "23"],
                    DC3000_vs_mock_30C = meta$sample[meta$temperature == "30"],
                    Temp_30_vs_23_mock = meta$sample[meta$treatment == "mock"],
                    Temp_30_vs_23_DC3000 = meta$sample[meta$treatment == "DC3000"])
  # select top DE genes (padj)
  genes <- res %>%
    as.data.frame() %>%
    rownames_to_column("Geneid") %>%
    filter(!is.na(padj), padj < 0.05) %>%
    arrange(padj) %>%
    pull(Geneid) %>%
    head(50)
  if (length(genes) == 0) {
    message("No significant DE genes for heatmap: ", nm)
    next
  }
  mat <- vsd_mat[genes, samples, drop = FALSE]
  mat <- mat - rowMeans(mat)
  ann <- meta %>% filter(sample %in% samples) %>% select(sample, treatment, temperature) %>% column_to_rownames("sample")
  # pdf
  pheatmap(mat, annotation_col = ann, show_rownames = FALSE, main = contrast_titles[nm],
           filename = paste0("Plots_Generated/heatmaps/heatmap_", nm, ".pdf"))
  # png
  save_pheatmap_png(mat, ann, paste0("Plots_Generated/heatmaps/heatmap_", nm, ".png"), main = contrast_titles[nm])
}

############################################################
## Defense gene boxplots (exact 4-condition panels; max 6 genes)
############################################################
defense_symbols <- c("PR1", "PR2", "PR5", "ICS1", "ALD1", "FMO1")

# map SYMBOL -> TAIR, ensure one TAIR per SYMBOL and present in matrix
gene_map_raw <- AnnotationDbi::select(org.At.tair.db, keys = defense_symbols, keytype = "SYMBOL", columns = "TAIR")
gene_map_clean <- gene_map_raw %>%
  filter(!is.na(TAIR)) %>%
  filter(TAIR %in% rownames(vsd_mat)) %>%
  group_by(SYMBOL) %>%
  slice(1) %>%         # deterministic pick if multiple TAIRs (first)
  ungroup()

if (nrow(gene_map_clean) == 0) {
  warning("No defense genes found in expression matrix. Skipping defense plots.")
} else {
  if (nrow(gene_map_clean) < length(defense_symbols)) {
    message("Some defense symbols not found or not expressed: ", paste(setdiff(defense_symbols, gene_map_clean$SYMBOL), collapse = ", "))
  }
  stopifnot(nrow(gene_map_clean) <= 6)
  df_def <- vsd_mat[gene_map_clean$TAIR, , drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("TAIR") %>%
    pivot_longer(-TAIR, names_to = "sample", values_to = "expression") %>%
    left_join(meta, by = "sample") %>%
    left_join(gene_map_clean, by = c("TAIR" = "TAIR")) %>%
    mutate(SYMBOL = ifelse(is.na(SYMBOL), TAIR, SYMBOL),
           treatment = factor(treatment, levels = c("mock", "DC3000")),
           temperature = factor(temperature, levels = c("23", "30")))
  p_def <- ggplot(df_def, aes(x = treatment, y = expression, fill = temperature)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                size = 1.2, alpha = 0.8) +
    facet_wrap(~ SYMBOL, scales = "free_y", ncol = 3) +
    labs(title = "Defense-related genes (4-condition panels)", x = "Treatment", y = "VST expression", fill = "Temperature")
  ggsave("Plots_Generated/boxplots/defense_genes.pdf", p_def, width = 10, height = 6)
  ggsave("Plots_Generated/boxplots/defense_genes.png", p_def, dpi = 300, width = 10, height = 6)
}

############################################################
## Top expressed genes per contrast (max 8 each) — now for ALL contrasts
############################################################
plot_top_genes_contrast <- function(samples, xvar, title, out_prefix, topN = 8) {
  mat <- vsd_mat[, samples, drop = FALSE]
  if (nrow(mat) == 0) {
    message("No rows in vsd_mat for these samples: ", out_prefix); return(NULL)
  }
  topN <- min(topN, nrow(mat))
  top_genes <- rownames(mat)[order(rowMeans(mat), decreasing = TRUE)][1:topN]
  df <- mat[top_genes, , drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("Geneid") %>%
    pivot_longer(-Geneid, names_to = "sample", values_to = "expression") %>%
    left_join(meta, by = "sample")
  p <- ggplot(df, aes_string(x = xvar, y = "expression", fill = xvar)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.2, alpha = 0.8) +
    facet_wrap(~ Geneid, scales = "free_y", ncol = min(4, topN)) +
    labs(title = title, x = xvar, y = "VST expression")
  ggsave(paste0(out_prefix, ".pdf"), p, width = 12, height = 6)
  ggsave(paste0(out_prefix, ".png"), p, dpi = 300, width = 12, height = 6)
}

# Generate top gene boxplots for all contrasts
plot_top_genes_contrast(meta$sample[meta$temperature == "23"],
                        "treatment",
                        "Top expressed genes — Infection at 23°C",
                        "Plots_Generated/boxplots/top_genes_23C",
                        topN = 8)

plot_top_genes_contrast(meta$sample[meta$temperature == "30"],
                        "treatment",
                        "Top expressed genes — Infection at 30°C",
                        "Plots_Generated/boxplots/top_genes_30C",
                        topN = 8)

plot_top_genes_contrast(meta$sample[meta$treatment == "mock"],
                        "temperature",
                        "Top expressed genes — Temperature effect in mock",
                        "Plots_Generated/boxplots/top_genes_mock_temp",
                        topN = 8)

plot_top_genes_contrast(meta$sample[meta$treatment == "DC3000"],
                        "temperature",
                        "Top expressed genes — Temperature effect in DC3000",
                        "Plots_Generated/boxplots/top_genes_DC3000_temp",
                        topN = 8)

############################################################
## GO enrichment — clearer filenames and titles
############################################################
tair2entrez <- AnnotationDbi::select(
  org.At.tair.db,
  keys = rownames(dds),
  keytype = "TAIR",
  columns = "ENTREZID"
) %>% filter(!is.na(ENTREZID))

universe_entrez <- unique(tair2entrez$ENTREZID)

for (nm in names(res_list)) {
  
  res <- res_list[[nm]]
  sig <- res %>%
    as.data.frame() %>%
    rownames_to_column("Geneid") %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    pull(Geneid)
  
  if (length(sig) == 0) {
    message("No sig genes for GO: ", nm); next
  }
  
  entrez_genes <- tair2entrez$ENTREZID[tair2entrez$TAIR %in% sig] %>% unique()
  if (length(entrez_genes) == 0) {
    message("No ENTREZ mappings for sig genes in contrast: ", nm); next
  }
  
  ego <- enrichGO(
    gene = entrez_genes,
    OrgDb = org.At.tair.db,
    keyType = "ENTREZID",
    universe = universe_entrez,
    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  if (is.null(ego) || nrow(ego@result) == 0) {
    message("No enriched GO terms for contrast: ", nm); next
  }
  
  # use human title and safe filename
  human_title <- contrast_titles[nm]
  safe_name <- gsub("[^A-Za-z0-9_]", "_", nm)
  
  p_go <- dotplot(ego, showCategory = 10) + ggtitle(human_title)
  ggsave(paste0("Plots_Generated/GO/GO_", safe_name, ".pdf"), p_go, width = 8, height = 6)
  ggsave(paste0("Plots_Generated/GO/GO_", safe_name, ".png"), p_go, dpi = 300, width = 8, height = 6)
}

############################################################
## Done
############################################################
