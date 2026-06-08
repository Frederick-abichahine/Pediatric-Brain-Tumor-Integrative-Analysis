##################################################################
# Name:                     Frederick Abi Chahine
# Masters:                  Quantitative and Computational Biology
# Matriculation Number:     256293
# Course:                   Network-based Data Analysis
# Dataset:                  GSE50161
##################################################################

##################################################
# Loading libraries and creating necessary folders
##################################################

library(GEOquery)
library(Biobase)
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(randomForest)
library(MASS)
library(glmnet)
library(caret)
library(gprofiler2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(pathfindR)
library(STRINGdb)
library(rScudo)
library(igraph)
library(RCy3)

set.seed(123)

if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
  if (!is.null(current_file) && nzchar(current_file)) {
    setwd(dirname(current_file))
  }
}

GEO_ACCESSION <- "GSE50161"
results_dir <- file.path(getwd(), "results", GEO_ACCESSION)
fig_dir <- file.path(results_dir, "figures")
table_dir <- file.path(results_dir, "tables")
object_dir <- file.path(results_dir, "objects")
cytoscape_dir <- file.path(results_dir, "cytoscape")
geo_cache_dir <- file.path(results_dir, "geo_cache")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cytoscape_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(geo_cache_dir, recursive = TRUE, showWarnings = FALSE)

################################################
# Creating helper functions to reduce redundancy
################################################

write_csv <- function(x, file_name, row.names = FALSE) {
  if (is.data.frame(x)) {
    list_columns <- vapply(x, is.list, logical(1))
    if (any(list_columns)) {
      for (nm in names(x)[list_columns]) {
        x[[nm]] <- vapply(x[[nm]], function(value) {
          value <- unlist(value, recursive = TRUE, use.names = FALSE)
          value <- as.character(value)
          value <- value[nzchar(value)]
          paste(value, collapse = ";")
        }, character(1))
      }
    }
  }
  write.csv(x, file.path(table_dir, file_name), row.names = row.names)
}

open_png <- function(file_name, width = 1800, height = 1400, res = 180) {
  png(file.path(fig_dir, file_name), width = width, height = height, res = res)
}

close_png <- function() {
  invisible(dev.off())
}

row_var <- function(x) {
  apply(x, 1, var, na.rm = TRUE)
}

row_zscore <- function(x) {
  z <- t(scale(t(x)))
  z[!is.finite(z)] <- 0
  z
}

first_token <- function(x, sep = "///") {
  x <- as.character(x)
  x[is.na(x)] <- ""
  out <- vapply(strsplit(x, sep, fixed = TRUE), function(parts) {
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0) NA_character_ else parts[[1]]
  }, character(1))
  out
}

clean_symbol <- function(x) {
  out <- first_token(x)
  out <- gsub("\\s+", "", out)
  out[out == ""] <- NA_character_
  out
}

feature_column <- function(feature_data, name) {
  if (name %in% colnames(feature_data)) {
    feature_data[[name]]
  } else {
    rep(NA_character_, nrow(feature_data))
  }
}

parse_group <- function(sample_title) {
  prefix <- toupper(sub("-.*$", "", sample_title))
  group <- ifelse(prefix == "NORMAL", "Normal", prefix)
  ifelse(group %in% c("Normal", "EPN", "GBM", "MED", "PA"), group, NA)
}

classification_metrics <- function(predicted, truth) {
  truth <- factor(truth)
  predicted <- factor(predicted, levels = levels(truth))
  cm <- table(Predicted = predicted, Actual = truth)
  accuracy <- sum(diag(cm)) / sum(cm)

  sensitivity <- numeric(length(levels(truth)))
  specificity <- numeric(length(levels(truth)))
  names(sensitivity) <- levels(truth)
  names(specificity) <- levels(truth)

  for (lvl in levels(truth)) {
    tp <- cm[lvl, lvl]
    fn <- sum(cm[, lvl]) - tp
    fp <- sum(cm[lvl, ]) - tp
    tn <- sum(cm) - tp - fn - fp
    sensitivity[lvl] <- ifelse(tp + fn == 0, NA, tp / (tp + fn))
    specificity[lvl] <- ifelse(tn + fp == 0, NA, tn / (tn + fp))
  }

  list(
    confusion = cm,
    accuracy = accuracy,
    balanced_accuracy = mean(sensitivity, na.rm = TRUE),
    by_class = data.frame(
      Class = levels(truth),
      Sensitivity = as.numeric(sensitivity),
      Specificity = as.numeric(specificity)
    )
  )
}

make_folds <- function(y, k = 5) {
  y <- factor(y)
  folds <- vector("list", k)
  for (lvl in levels(y)) {
    idx <- sample(which(y == lvl))
    fold_id <- rep(seq_len(k), length.out = length(idx))
    for (i in seq_len(k)) {
      folds[[i]] <- c(folds[[i]], idx[fold_id == i])
    }
  }
  lapply(folds, sort)
}

feature_f_score <- function(x, y) {
  y <- factor(y)
  overall_mean <- colMeans(x)
  between <- rep(0, ncol(x))
  within <- rep(0, ncol(x))

  for (lvl in levels(y)) {
    group_x <- x[y == lvl, , drop = FALSE]
    group_mean <- colMeans(group_x)
    between <- between + nrow(group_x) * (group_mean - overall_mean)^2
    within <- within + colSums(sweep(group_x, 2, group_mean, "-")^2)
  }

  between / pmax(within, .Machine$double.eps)
}

count_gene_list <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(strsplit(x, ",", fixed = TRUE), function(z) {
    z <- trimws(z)
    z <- z[nzchar(z)]
    length(unique(z))
  }, integer(1))
}

plot_confusion_matrix <- function(cm, title, file_name) {
  cm_df <- as.data.frame(cm)
  colnames(cm_df) <- c("Predicted", "Actual", "N")

  p <- ggplot(cm_df, aes(x = Actual, y = Predicted, fill = N)) +
    geom_tile(color = "white") +
    geom_text(aes(label = N), size = 4) +
    scale_fill_gradient(low = "#f7fbff", high = "#08519c") +
    theme_minimal(base_size = 12) +
    labs(title = title, x = "Actual class", y = "Predicted class")

  ggsave(file.path(fig_dir, file_name), p, width = 7, height = 5, dpi = 180)
}

##############################################
# Loading GSE50161 dataset and sample metadata
##############################################

gse_list <- getGEO(
  GEO_ACCESSION,
  GSEMatrix = TRUE,
  AnnotGPL = TRUE,
  destdir = geo_cache_dir
)

gset <- gse_list[[1]]
expr_raw <- exprs(gset)
pheno <- pData(gset)
feature_data <- fData(gset)
pheno <- pheno[colnames(expr_raw), , drop = FALSE]
pheno$Sample <- rownames(pheno)
pheno$Title <- as.character(pheno$title)
pheno$Group <- parse_group(pheno$Title)

if (any(is.na(pheno$Group))) {
  stop("Some samples could not be assigned to Normal, EPN, GBM, MED, or PA.")
}

pheno$Group <- factor(pheno$Group, levels = c("Normal", "EPN", "GBM", "MED", "PA"))
pheno$TumorStatus <- factor(
  ifelse(pheno$Group == "Normal", "Normal", "Tumor"),
  levels = c("Normal", "Tumor")
)

sample_counts <- as.data.frame(table(pheno$Group))
colnames(sample_counts) <- c("Group", "N")
write_csv(sample_counts, "01_sample_group_counts.csv")

metadata_export <- pheno[, c(
  "Sample",
  "Title",
  "Group",
  "TumorStatus",
  "disease state:ch1",
  "tissue:ch1"
)]
write_csv(metadata_export, "01_sample_metadata.csv")

##############################################
# Performing preprocessing and quality control
##############################################

expr_raw <- as.matrix(expr_raw)
storage.mode(expr_raw) <- "numeric"

raw_quantiles <- quantile(
  expr_raw,
  probs = c(0, 0.25, 0.50, 0.75, 0.99, 1),
  na.rm = TRUE
)

needs_log2 <- (raw_quantiles[[5]] > 100) ||
  ((raw_quantiles[[6]] - raw_quantiles[[1]] > 50) && raw_quantiles[[2]] > 0)

expr_processed <- expr_raw
log2_applied <- FALSE

if (needs_log2) {
  expr_processed[expr_processed <= 0] <- NA
  expr_processed <- log2(expr_processed)
  log2_applied <- TRUE
}

complete_rows <- apply(expr_processed, 1, function(x) all(is.finite(x)))
expr_processed <- expr_processed[complete_rows, , drop = FALSE]
gset <- gset[rownames(expr_processed), ]
feature_data <- fData(gset)

median_range <- diff(range(apply(expr_processed, 2, median)))
iqr_range <- diff(range(apply(expr_processed, 2, IQR)))
normalization_applied <- median_range > 1 || iqr_range > 1.5

if (normalization_applied) {
  expr_processed <- normalizeBetweenArrays(expr_processed, method = "quantile")
}

preprocess_summary <- data.frame(
  Metric = c(
    "Raw min",
    "Raw Q1",
    "Raw median",
    "Raw Q3",
    "Raw 99th percentile",
    "Raw max",
    "Log2 transform applied",
    "Between-array normalization applied",
    "Rows retained after missing-value filtering",
    "Samples"
  ),
  Value = c(
    as.character(round(raw_quantiles, 4)),
    as.character(log2_applied),
    as.character(normalization_applied),
    as.character(nrow(expr_processed)),
    as.character(ncol(expr_processed))
  )
)
write_csv(preprocess_summary, "02_preprocessing_summary.csv")

open_png("01_boxplot_processed_expression.png", width = 2200, height = 1400)
boxplot(
  expr_processed,
  outline = FALSE,
  las = 2,
  col = as.integer(pheno$Group),
  main = "GSE50161 processed expression distributions",
  ylab = "Expression"
)
legend("topright", legend = levels(pheno$Group), fill = seq_along(levels(pheno$Group)), cex = 0.8)
close_png()

open_png("02_density_processed_expression.png", width = 1800, height = 1300)
plotDensities(expr_processed, group = pheno$Group, main = "Expression density by group")
close_png()

########################################
# Annotating probe and gene level matrix
########################################

probe_annotation <- data.frame(
  ProbeID = rownames(expr_processed),
  GeneSymbol = clean_symbol(feature_column(feature_data, "Gene symbol")),
  EntrezID = first_token(feature_column(feature_data, "Gene ID")),
  GeneTitle = first_token(feature_column(feature_data, "Gene title")),
  GOProcess = feature_column(feature_data, "GO:Process"),
  GOProcessID = feature_column(feature_data, "GO:Process ID"),
  GOFunction = feature_column(feature_data, "GO:Function"),
  GOFunctionID = feature_column(feature_data, "GO:Function ID"),
  GOComponent = feature_column(feature_data, "GO:Component"),
  GOComponentID = feature_column(feature_data, "GO:Component ID")
)

annotated <- !is.na(probe_annotation$GeneSymbol) & nzchar(probe_annotation$GeneSymbol)
expr_annotated <- expr_processed[annotated, , drop = FALSE]
probe_annotation <- probe_annotation[annotated, , drop = FALSE]

probe_variance <- row_var(expr_annotated)
best_probe_order <- order(probe_annotation$GeneSymbol, -probe_variance)
best_probe_idx <- best_probe_order[!duplicated(probe_annotation$GeneSymbol[best_probe_order])]

expr_gene <- expr_annotated[best_probe_idx, , drop = FALSE]
gene_annotation <- probe_annotation[best_probe_idx, , drop = FALSE]
rownames(expr_gene) <- gene_annotation$GeneSymbol
rownames(gene_annotation) <- gene_annotation$GeneSymbol
gene_annotation$ProbeVariance <- probe_variance[best_probe_idx]

write_csv(probe_annotation, "03_probe_annotation_all.csv")
write_csv(gene_annotation, "03_gene_annotation_selected_probe_per_gene.csv")

save(
  expr_raw,
  expr_processed,
  expr_gene,
  pheno,
  gene_annotation,
  file = file.path(object_dir, "01_expression_objects.RData")
)

####################################
# Performing differential expression
####################################

FDR_CUTOFF <- 0.05
LOGFC_CUTOFF <- 1

make_deg_table <- function(fit, coef_name) {
  tab <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH")
  tab$GeneSymbol <- rownames(tab)
  tab$EntrezID <- gene_annotation[tab$GeneSymbol, "EntrezID"]
  tab$GeneTitle <- gene_annotation[tab$GeneSymbol, "GeneTitle"]
  tab$Status <- "Not significant"
  tab$Status[tab$adj.P.Val < FDR_CUTOFF & tab$logFC >= LOGFC_CUTOFF] <- "Up"
  tab$Status[tab$adj.P.Val < FDR_CUTOFF & tab$logFC <= -LOGFC_CUTOFF] <- "Down"
  tab
}

plot_volcano <- function(tab, contrast_name) {
  tab$MinusLog10FDR <- -log10(pmax(tab$adj.P.Val, .Machine$double.xmin))
  tab$Label <- ""

  up_idx <- which(tab$Status == "Up")
  down_idx <- which(tab$Status == "Down")
  top_up <- up_idx[order(tab$adj.P.Val[up_idx])][seq_len(min(10, length(up_idx)))]
  top_down <- down_idx[order(tab$adj.P.Val[down_idx])][seq_len(min(10, length(down_idx)))]
  label_idx <- unique(c(top_up, top_down))
  tab$Label[label_idx] <- tab$GeneSymbol[label_idx]

  p <- ggplot(tab, aes(logFC, MinusLog10FDR, color = Status)) +
    geom_point(alpha = 0.70, size = 1.25) +
    geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = tab[nzchar(tab$Label), ],
      aes(label = Label),
      size = 3,
      min.segment.length = 0,
      max.overlaps = 40
    ) +
    scale_color_manual(values = c("Down" = "#2C7FB8", "Not significant" = "grey70", "Up" = "#D95F02")) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Volcano plot:", contrast_name),
      x = "log2 fold change",
      y = "-log10 adjusted P value"
    )

  open_png(paste0("03_volcano_", contrast_name, ".png"), width = 1350, height = 1080)
  print(p)
  close_png()
}

status_design <- model.matrix(~ 0 + pheno$TumorStatus)
colnames(status_design) <- levels(pheno$TumorStatus)
status_fit <- lmFit(expr_gene, status_design)
status_contrast <- makeContrasts(TumorVsNormal = Tumor - Normal, levels = status_design)
status_fit <- eBayes(contrasts.fit(status_fit, status_contrast))

tumor_vs_normal <- make_deg_table(status_fit, "TumorVsNormal")
write_csv(tumor_vs_normal, "04_DEG_Tumor_vs_Normal_limma.csv")
plot_volcano(tumor_vs_normal, "TumorVsNormal")

group_design <- model.matrix(~ 0 + pheno$Group)
colnames(group_design) <- levels(pheno$Group)
group_fit <- lmFit(expr_gene, group_design)
group_contrast <- makeContrasts(
  EPNvsNormal = EPN - Normal,
  GBMvsNormal = GBM - Normal,
  MEDvsNormal = MED - Normal,
  PAvsNormal = PA - Normal,
  levels = group_design
)
group_fit <- eBayes(contrasts.fit(group_fit, group_contrast))

deg_tables <- list()
for (contrast_name in colnames(group_contrast)) {
  deg_tables[[contrast_name]] <- make_deg_table(group_fit, contrast_name)
  write_csv(deg_tables[[contrast_name]], paste0("04_DEG_", contrast_name, "_limma.csv"))
  plot_volcano(deg_tables[[contrast_name]], contrast_name)
}

deg_counts <- data.frame(
  Contrast = c("TumorVsNormal", names(deg_tables)),
  Up = c(
    sum(tumor_vs_normal$Status == "Up"),
    vapply(deg_tables, function(x) sum(x$Status == "Up"), integer(1))
  ),
  Down = c(
    sum(tumor_vs_normal$Status == "Down"),
    vapply(deg_tables, function(x) sum(x$Status == "Down"), integer(1))
  ),
  Total = c(
    sum(tumor_vs_normal$Status %in% c("Up", "Down")),
    vapply(deg_tables, function(x) sum(x$Status %in% c("Up", "Down")), integer(1))
  )
)
write_csv(deg_counts, "04_DEG_counts.csv")

annotation_col <- data.frame(Group = pheno$Group, TumorStatus = pheno$TumorStatus)
rownames(annotation_col) <- pheno$Sample

top_heatmap_genes <- head(tumor_vs_normal$GeneSymbol[order(tumor_vs_normal$adj.P.Val)], 50)
open_png("04_heatmap_top50_DEG_TumorVsNormal.png", width = 1800, height = 1600)
pheatmap(
  row_zscore(expr_gene[top_heatmap_genes, , drop = FALSE]),
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 6,
  main = "Top 50 tumor vs normal genes"
)
close_png()

############################################
# Performing PCA and unsupervised clustering
############################################

gene_variance <- row_var(expr_gene)
top_pca_genes <- names(sort(gene_variance, decreasing = TRUE))[1:2000]

pca <- prcomp(t(expr_gene[top_pca_genes, , drop = FALSE]), center = TRUE, scale. = TRUE)
pca_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = pheno$Group,
  TumorStatus = pheno$TumorStatus
)
write_csv(pca_df, "05_PCA_coordinates.csv")
write_csv(data.frame(PC = seq_along(pca_var), VariancePercent = pca_var), "05_PCA_scree_variance.csv")

group_colors <- c(
  "Normal" = "#3C3C3C",
  "EPN" = "#2C7FB8",
  "GBM" = "#D95F02",
  "MED" = "#1B9E77",
  "PA" = "#7570B3"
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.90) +
  scale_color_manual(values = group_colors) +
  theme_minimal(base_size = 12) +
  labs(
    title = "PCA of GSE50161",
    x = paste0("PC1: ", round(pca_var[1], 2), "%"),
    y = paste0("PC2: ", round(pca_var[2], 2), "%")
  )
ggsave(file.path(fig_dir, "05_PCA_group_PC1_PC2.png"), p_pca, width = 7, height = 5.5, dpi = 180)

p_scree <- ggplot(data.frame(PC = 1:25, VariancePercent = pca_var[1:25]),
                  aes(PC, VariancePercent)) +
  geom_line(color = "#555555") +
  geom_point(color = "#1B9E77") +
  theme_minimal(base_size = 12) +
  labs(title = "PCA scree plot", x = "Principal component", y = "Variance explained (%)")
ggsave(file.path(fig_dir, "06_PCA_scree_plot.png"), p_scree, width = 7, height = 5, dpi = 180)

kmeans_input <- pca$x[, 1:10, drop = FALSE]
for (k in c(2, 5)) {
  km <- kmeans(kmeans_input, centers = k, nstart = 50)
  cluster <- factor(paste0("Cluster_", km$cluster))

  write_csv(as.data.frame.matrix(table(Cluster = cluster, Group = pheno$Group)),
            paste0("06_kmeans_k", k, "_vs_group.csv"), row.names = TRUE)
  write_csv(as.data.frame.matrix(table(Cluster = cluster, TumorStatus = pheno$TumorStatus)),
            paste0("06_kmeans_k", k, "_vs_tumor_status.csv"), row.names = TRUE)

  km_df <- pca_df
  km_df$Cluster <- cluster
  p_km <- ggplot(km_df, aes(PC1, PC2, color = Cluster, shape = Group)) +
    geom_point(size = 3, alpha = 0.90) +
    theme_minimal(base_size = 12) +
    labs(title = paste("K-means clustering on PCA space, k =", k))
  ggsave(file.path(fig_dir, paste0("07_kmeans_k", k, "_on_PCA.png")),
         p_km, width = 7.5, height = 5.5, dpi = 180)
}

sample_cor <- cor(expr_gene[top_pca_genes, , drop = FALSE], method = "pearson")
sample_hc <- hclust(as.dist(1 - sample_cor), method = "ward.D2")

open_png("08_hierarchical_clustering_dendrogram.png", width = 2200, height = 1300)
plot(sample_hc, labels = pheno$Title, main = "Hierarchical clustering of GSE50161 samples", xlab = "", sub = "", cex = 0.55)
close_png()

open_png("09_sample_correlation_heatmap.png", width = 1700, height = 1600)
pheatmap(
  sample_cor,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Sample-sample Pearson correlation"
)
close_png()

######################################
# Performing supervised classification
######################################

top_ml_genes <- names(sort(gene_variance, decreasing = TRUE))[1:1000]
feature_map <- data.frame(
  GeneSymbol = top_ml_genes,
  FeatureName = make.names(top_ml_genes, unique = TRUE)
)

x_all <- as.data.frame(t(expr_gene[top_ml_genes, , drop = FALSE]), check.names = FALSE)
colnames(x_all) <- feature_map$FeatureName
y_group <- pheno$Group
folds <- make_folds(y_group, k = 5)

# Random Forest
rf_pred <- rep(NA_character_, length(y_group))
for (i in seq_along(folds)) {
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_along(y_group), test_idx)
  rf_fit <- randomForest(
    x = x_all[train_idx, , drop = FALSE],
    y = y_group[train_idx],
    ntree = 1000
  )
  rf_pred[test_idx] <- as.character(predict(rf_fit, x_all[test_idx, , drop = FALSE]))
}

rf_metrics <- classification_metrics(rf_pred, y_group)
write_csv(as.data.frame.matrix(rf_metrics$confusion), "07_RF_crossvalidated_confusion_matrix.csv", row.names = TRUE)
write_csv(rf_metrics$by_class, "07_RF_crossvalidated_by_class_metrics.csv")
write_csv(data.frame(Model = "Random Forest", Accuracy = rf_metrics$accuracy,
                     BalancedAccuracy = rf_metrics$balanced_accuracy),
          "07_RF_crossvalidated_overall_metrics.csv")
plot_confusion_matrix(rf_metrics$confusion, "Random Forest cross-validated confusion matrix", "10_RF_crossvalidated_confusion_matrix.png")

rf_final <- randomForest(x = x_all, y = y_group, ntree = 1000, importance = TRUE, keep.forest = TRUE)
saveRDS(rf_final, file.path(object_dir, "02_random_forest_final_model.rds"))

rf_oob <- data.frame(Trees = seq_len(nrow(rf_final$err.rate)), OOBError = rf_final$err.rate[, "OOB"])
write_csv(rf_oob, "07_RF_OOB_error_by_tree.csv")

p_oob <- ggplot(rf_oob, aes(Trees, OOBError)) +
  geom_line(color = "#2C7FB8", linewidth = 0.8) +
  theme_minimal(base_size = 12) +
  labs(title = "Random Forest OOB error", x = "Trees", y = "OOB error")
ggsave(file.path(fig_dir, "11_RF_OOB_error.png"), p_oob, width = 7, height = 5, dpi = 180)

rf_importance <- data.frame(FeatureName = rownames(importance(rf_final)),
                            importance(rf_final),
                            row.names = NULL,
                            check.names = FALSE)
rf_importance <- merge(rf_importance, feature_map, by = "FeatureName", all.x = TRUE)
rf_importance <- rf_importance[order(rf_importance$MeanDecreaseGini, decreasing = TRUE), ]
write_csv(rf_importance, "07_RF_feature_importance.csv")

top_rf <- head(rf_importance, 25)
p_rf <- ggplot(top_rf, aes(reorder(GeneSymbol, MeanDecreaseGini), MeanDecreaseGini)) +
  geom_col(fill = "#1B9E77") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(title = "Top Random Forest genes", x = "Gene", y = "Mean decrease Gini")
ggsave(file.path(fig_dir, "12_RF_top25_feature_importance.png"), p_rf, width = 7, height = 6, dpi = 180)

open_png("13_RF_top25_gene_heatmap.png", width = 1700, height = 1500)
pheatmap(
  row_zscore(expr_gene[top_rf$GeneSymbol, , drop = FALSE]),
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 8,
  main = "Top 25 Random Forest genes"
)
close_png()

# LDA
lda_pred <- rep(NA_character_, length(y_group))
lda_feature_count <- min(50, ncol(x_all), nrow(x_all) - length(levels(y_group)) - 2)

for (i in seq_along(folds)) {
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_along(y_group), test_idx)
  f_scores <- feature_f_score(as.matrix(x_all[train_idx, , drop = FALSE]), y_group[train_idx])
  selected <- names(sort(f_scores, decreasing = TRUE))[1:lda_feature_count]
  lda_fit <- lda(x = x_all[train_idx, selected, drop = FALSE], grouping = y_group[train_idx])
  lda_pred[test_idx] <- as.character(predict(lda_fit, x_all[test_idx, selected, drop = FALSE])$class)
}

lda_metrics <- classification_metrics(lda_pred, y_group)
write_csv(as.data.frame.matrix(lda_metrics$confusion), "08_LDA_crossvalidated_confusion_matrix.csv", row.names = TRUE)
write_csv(lda_metrics$by_class, "08_LDA_crossvalidated_by_class_metrics.csv")
write_csv(data.frame(Model = "LDA", Accuracy = lda_metrics$accuracy,
                     BalancedAccuracy = lda_metrics$balanced_accuracy),
          "08_LDA_crossvalidated_overall_metrics.csv")
plot_confusion_matrix(lda_metrics$confusion, "LDA cross-validated confusion matrix", "14_LDA_crossvalidated_confusion_matrix.png")

lda_scores <- feature_f_score(as.matrix(x_all), y_group)
lda_features <- names(sort(lda_scores, decreasing = TRUE))[1:lda_feature_count]
lda_final <- lda(x = x_all[, lda_features, drop = FALSE], grouping = y_group)
saveRDS(lda_final, file.path(object_dir, "03_LDA_final_model.rds"))

lda_projection <- predict(lda_final, x_all[, lda_features, drop = FALSE])
lda_df <- data.frame(
  Sample = rownames(x_all),
  LD1 = lda_projection$x[, 1],
  LD2 = lda_projection$x[, 2],
  Group = y_group
)
write_csv(lda_df, "08_LDA_coordinates_final_model.csv")

p_lda <- ggplot(lda_df, aes(LD1, LD2, color = Group)) +
  geom_point(size = 3, alpha = 0.90) +
  scale_color_manual(values = group_colors) +
  theme_minimal(base_size = 12) +
  labs(title = "LDA projection", x = "LD1", y = "LD2")
ggsave(file.path(fig_dir, "15_LDA_projection.png"), p_lda, width = 7, height = 5.5, dpi = 180)

lda_coef <- data.frame(FeatureName = rownames(lda_final$scaling),
                       lda_final$scaling,
                       row.names = NULL,
                       check.names = FALSE)
lda_coef <- merge(lda_coef, feature_map, by = "FeatureName", all.x = TRUE)
write_csv(lda_coef, "08_LDA_feature_coefficients.csv")

# LASSO and Ridge
x_mat <- as.matrix(x_all)
run_glmnet_model <- function(alpha_value, model_name) {
  pred <- rep(NA_character_, length(y_group))
  lambdas <- numeric(length(folds))

  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(y_group), test_idx)
    cv_fit <- cv.glmnet(
      x = x_mat[train_idx, , drop = FALSE],
      y = y_group[train_idx],
      family = "multinomial",
      alpha = alpha_value,
      nfolds = min(5, min(table(y_group[train_idx]))),
      type.measure = "class"
    )
    lambdas[i] <- cv_fit$lambda.min
    pred[test_idx] <- as.character(predict(cv_fit, x_mat[test_idx, , drop = FALSE], s = "lambda.min", type = "class"))
  }

  metrics <- classification_metrics(pred, y_group)
  write_csv(as.data.frame.matrix(metrics$confusion),
            paste0("09_", model_name, "_crossvalidated_confusion_matrix.csv"), row.names = TRUE)
  write_csv(metrics$by_class, paste0("09_", model_name, "_crossvalidated_by_class_metrics.csv"))
  plot_confusion_matrix(metrics$confusion, paste(model_name, "cross-validated confusion matrix"),
                        paste0("16_", model_name, "_crossvalidated_confusion_matrix.png"))

  final_fit <- cv.glmnet(
    x = x_mat,
    y = y_group,
    family = "multinomial",
    alpha = alpha_value,
    nfolds = min(5, min(table(y_group))),
    type.measure = "class"
  )
  saveRDS(final_fit, file.path(object_dir, paste0("04_", model_name, "_glmnet_model.rds")))

  coef_df <- do.call(rbind, lapply(names(coef(final_fit, s = "lambda.min")), function(class_name) {
    mat <- as.matrix(coef(final_fit, s = "lambda.min")[[class_name]])
    data.frame(Class = class_name, FeatureName = rownames(mat), Coefficient = mat[, 1])
  }))
  coef_df <- subset(coef_df, FeatureName != "(Intercept)")
  coef_df <- merge(coef_df, feature_map, by = "FeatureName", all.x = TRUE)
  coef_df <- coef_df[order(abs(coef_df$Coefficient), decreasing = TRUE), ]
  write_csv(coef_df, paste0("09_", model_name, "_coefficients.csv"))

  data.frame(
    Model = model_name,
    Accuracy = metrics$accuracy,
    BalancedAccuracy = metrics$balanced_accuracy,
    LambdaMedian = median(lambdas)
  )
}

glmnet_results <- rbind(
  run_glmnet_model(alpha_value = 1, model_name = "LASSO"),
  run_glmnet_model(alpha_value = 0, model_name = "Ridge")
)

model_comparison <- rbind(
  data.frame(Model = "Random Forest", Accuracy = rf_metrics$accuracy, BalancedAccuracy = rf_metrics$balanced_accuracy),
  data.frame(Model = "LDA", Accuracy = lda_metrics$accuracy, BalancedAccuracy = lda_metrics$balanced_accuracy),
  glmnet_results[, c("Model", "Accuracy", "BalancedAccuracy")]
)
write_csv(model_comparison, "10_model_comparison_crossvalidated.csv")

p_models <- ggplot(model_comparison, aes(Model, Accuracy, group = 1)) +
  geom_point(size = 3, color = "#D95F02") +
  geom_line(color = "#555555") +
  ylim(0, 1) +
  theme_minimal(base_size = 12) +
  labs(title = "Cross-validated model comparison", y = "Accuracy", x = "")
ggsave(file.path(fig_dir, "17_model_comparison_accuracy.png"), p_models, width = 7, height = 5, dpi = 180)

# Repeated CV using caret
train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 3)
caret_models <- list(
  rf = train(x_all, y_group, method = "rf", trControl = train_control, ntree = 1000, tuneLength = 3),
  lda = train(x_all[, lda_features, drop = FALSE], y_group, method = "lda", trControl = train_control),
  glmnet = train(x_all, y_group, method = "glmnet", family = "multinomial",
                 trControl = train_control, tuneLength = 5)
)
saveRDS(caret_models, file.path(object_dir, "05_caret_model_comparison.rds"))
caret_acc <- do.call(rbind, lapply(names(caret_models), function(nm) {
  data.frame(Model = nm, Accuracy = max(caret_models[[nm]]$results$Accuracy, na.rm = TRUE))
}))
write_csv(caret_acc, "10_caret_repeated_cv_model_comparison.csv")

#################################
# Obtaining rank-based signatures
#################################

centered_expr <- sweep(expr_gene, 1, rowMeans(expr_gene), "-")
signature_length <- 100

up_signatures <- apply(centered_expr, 2, function(x) {
  names(sort(x, decreasing = TRUE))[1:signature_length]
})

down_signatures <- apply(centered_expr, 2, function(x) {
  names(sort(x, decreasing = FALSE))[1:signature_length]
})

signature_frequency <- function(signature_matrix, groups) {
  out <- do.call(rbind, lapply(levels(groups), function(grp) {
    genes <- as.vector(signature_matrix[, groups == grp, drop = FALSE])
    tab <- sort(table(genes), decreasing = TRUE)
    data.frame(Group = grp, GeneSymbol = names(tab), Frequency = as.integer(tab))
  }))
  rownames(out) <- NULL
  out
}

write_csv(signature_frequency(up_signatures, pheno$Group), "11_rank_signature_consensus_up_genes.csv")
write_csv(signature_frequency(down_signatures, pheno$Group), "11_rank_signature_consensus_down_genes.csv")

signature_similarity <- matrix(
  0,
  nrow = ncol(expr_gene),
  ncol = ncol(expr_gene),
  dimnames = list(colnames(expr_gene), colnames(expr_gene))
)

for (i in seq_len(ncol(expr_gene))) {
  sig_i <- c(up_signatures[, i], down_signatures[, i])
  for (j in i:ncol(expr_gene)) {
    sig_j <- c(up_signatures[, j], down_signatures[, j])
    sim <- length(intersect(sig_i, sig_j)) / length(union(sig_i, sig_j))
    signature_similarity[i, j] <- sim
    signature_similarity[j, i] <- sim
  }
}

write_csv(signature_similarity, "11_rank_signature_sample_jaccard_similarity.csv", row.names = TRUE)

open_png("18_rank_signature_similarity_heatmap.png", width = 1600, height = 1500)
pheatmap(
  signature_similarity,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Rank-signature sample similarity"
)
close_png()

try({
  scudo_genes <- names(sort(gene_variance, decreasing = TRUE))[1:5000]
  scudo_expr <- expr_gene[scudo_genes, , drop = FALSE]
  scudo_eset <- ExpressionSet(assayData = as.matrix(scudo_expr))

  train_idx <- as.vector(createDataPartition(y_group, p = 0.70, list = FALSE))
  test_idx <- setdiff(seq_along(y_group), train_idx)

  scudo_train <- scudoTrain(
    scudo_eset[, train_idx],
    groups = y_group[train_idx],
    nTop = 100,
    nBottom = 100,
    alpha = 0.05
  )

  scudo_graph <- scudoNetwork(scudo_train, N = 0.25)
  saveRDS(scudo_train, file.path(object_dir, "06_rScudo_train_results.rds"))
  write_graph(scudo_graph, file.path(cytoscape_dir, "rScudo_training_network.graphml"), format = "graphml")

  open_png("19_rScudo_training_network.png", width = 1600, height = 1300)
  scudoPlot(scudo_graph, vertex.label = NA)
  close_png()

  scudo_class <- scudoClassify(
    scudo_eset[, train_idx],
    scudo_eset[, test_idx],
    N = 0.25,
    nTop = 100,
    nBottom = 100,
    trainGroups = y_group[train_idx],
    alpha = 0.05
  )

  scudo_metrics <- confusionMatrix(
    factor(scudo_class$predicted, levels = levels(y_group)),
    y_group[test_idx]
  )
  capture.output(scudo_metrics, file = file.path(table_dir, "11_rScudo_caret_confusion_matrix.txt"))

  if (exists("consensusUpSignatures", where = asNamespace("rScudo"), mode = "function")) {
    write_csv(consensusUpSignatures(scudo_train), "11_rScudo_consensus_up_signatures.csv", row.names = TRUE)
    write_csv(consensusDownSignatures(scudo_train), "11_rScudo_consensus_down_signatures.csv", row.names = TRUE)
  }
}, silent = TRUE)

##################################
# Performing functional enrichment
##################################

primary_deg <- subset(tumor_vs_normal, adj.P.Val < FDR_CUTOFF & abs(logFC) >= LOGFC_CUTOFF)
primary_deg <- primary_deg[!duplicated(primary_deg$GeneSymbol), ]
sig_symbols <- primary_deg$GeneSymbol
write_csv(primary_deg, "12_primary_DEG_list_for_enrichment.csv")

# g:Profiler
gp <- gost(
  query = sig_symbols,
  organism = "hsapiens",
  correction_method = "fdr",
  significant = TRUE
)
gprof_results <- gp$result
write_csv(gprof_results, "12_gProfiler_enrichment.csv")

gp_plot <- gprof_results
gp_plot$p_value <- as.numeric(gp_plot$p_value)
gp_plot$intersection_size <- as.numeric(gp_plot$intersection_size)
gp_plot$source <- as.character(gp_plot$source)
gp_plot$source <- ifelse(gp_plot$source == "REAC", "Reactome", gp_plot$source)
gp_plot <- gp_plot[order(gp_plot$p_value), ]
gp_plot <- do.call(rbind, lapply(split(gp_plot, gp_plot$source), function(x) head(x, 5)))
gp_plot <- gp_plot[order(gp_plot$source, gp_plot$p_value), ]
write_csv(gp_plot, "12_gProfiler_enrichment_plot_data_by_source.csv")

p_gp <- ggplot(gp_plot, aes(source, -log10(p_value), size = intersection_size, color = source)) +
  geom_jitter(width = 0.20, alpha = 0.85) +
  theme_minimal(base_size = 12) +
  labs(
    title = "gProfiler enrichment summary by annotation source",
    x = "Database",
    y = "-log10 P value",
    size = "Genes"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

open_png("20_gProfiler_enrichment_summary.png", width = 1440, height = 990)
print(p_gp)
close_png()

# clusterProfiler and ReactomePA
gene_df <- bitr(sig_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- unique(gene_df$ENTREZID)

ego_bp <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)
write_csv(as.data.frame(ego_bp), "12_clusterProfiler_GO_BP.csv")

if (nrow(as.data.frame(ego_bp)) > 0) {
  p_go <- dotplot(ego_bp, showCategory = 20) + ggtitle("GO Biological Process enrichment")
  ggsave(file.path(fig_dir, "21_clusterProfiler_GO_BP_dotplot.png"), p_go, width = 8, height = 6, dpi = 180)
}

ekegg <- enrichKEGG(gene = entrez_ids, organism = "hsa", pvalueCutoff = 0.05)
write_csv(as.data.frame(ekegg), "12_clusterProfiler_KEGG.csv")

reactome_results <- enrichPathway(
  gene = entrez_ids,
  organism = "human",
  pvalueCutoff = 0.05,
  readable = TRUE
)
write_csv(as.data.frame(reactome_results), "12_ReactomePA_enrichment.csv")

# Simple GPL570 GO fallback
go_ids <- strsplit(as.character(gene_annotation$GOProcessID), "///", fixed = TRUE)
go_terms <- strsplit(as.character(gene_annotation$GOProcess), "///", fixed = TRUE)

go_map <- do.call(rbind, lapply(seq_len(nrow(gene_annotation)), function(i) {
  ids <- trimws(go_ids[[i]])
  terms <- trimws(go_terms[[i]])
  ids <- ids[nzchar(ids)]
  terms <- terms[nzchar(terms)]
  if (length(ids) == 0) return(NULL)
  if (length(terms) < length(ids)) terms <- c(terms, rep(NA_character_, length(ids) - length(terms)))
  data.frame(GeneSymbol = gene_annotation$GeneSymbol[i], GO_ID = ids, GO_Term = terms[seq_along(ids)])
}))

go_map <- unique(go_map)
universe <- rownames(expr_gene)
selected <- intersect(sig_symbols, universe)
go_split <- split(go_map$GeneSymbol, go_map$GO_ID)

go_results <- do.call(rbind, lapply(names(go_split), function(go_id) {
  term_genes <- unique(intersect(go_split[[go_id]], universe))
  if (length(term_genes) < 5 || length(term_genes) > 500) return(NULL)

  a <- length(intersect(selected, term_genes))
  b <- length(setdiff(selected, term_genes))
  c <- length(setdiff(term_genes, selected))
  d <- length(setdiff(universe, union(selected, term_genes)))

  p_value <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
  data.frame(GO_ID = go_id, SelectedGenesInTerm = a, TermSize = length(term_genes), P.Value = p_value)
}))

go_names <- go_map[!duplicated(go_map$GO_ID), c("GO_ID", "GO_Term")]
go_results <- merge(go_results, go_names, by = "GO_ID", all.x = TRUE)
go_results$adj.P.Val <- p.adjust(go_results$P.Value, method = "BH")
go_results <- go_results[order(go_results$adj.P.Val, go_results$P.Value), ]
write_csv(go_results, "12_fallback_GPL570_GO_BP_overrepresentation.csv")

top_go_fallback <- head(go_results, 20)
p_fallback <- ggplot(top_go_fallback, aes(-log10(adj.P.Val), reorder(GO_Term, -log10(adj.P.Val)))) +
  geom_col(fill = "#7570B3") +
  theme_minimal(base_size = 11) +
  labs(
    title = "Fallback GO-BP enrichment from GPL570 annotation",
    x = "-log10 adjusted P value",
    y = ""
  )
ggsave(file.path(fig_dir, "22_fallback_GO_BP_overrepresentation.png"), p_fallback, width = 8, height = 6, dpi = 180)

###########################################
# Performing network-based pathway analysis
###########################################

pathfindr_input <- primary_deg[, c("GeneSymbol", "logFC", "adj.P.Val")]
colnames(pathfindr_input) <- c("Gene.symbol", "logFC", "adj.P.Val")
pathfindr_input <- pathfindr_input[!duplicated(pathfindr_input$Gene.symbol), ]
write_csv(pathfindr_input, "13_pathfindR_input.csv")

pathfindr_results <- run_pathfindR(
  input = pathfindr_input,
  gene_sets = "Reactome",
  pin_name_path = "Biogrid",
  p_val_threshold = FDR_CUTOFF
)
write_csv(pathfindr_results, "13_pathfindR_Reactome_results.csv")

path_plot <- pathfindr_results
path_plot$lowest_p <- as.numeric(path_plot$lowest_p)
path_plot$Fold_Enrichment <- as.numeric(path_plot$Fold_Enrichment)
path_plot$UpGeneCount <- count_gene_list(path_plot$Up_regulated)
path_plot$DownGeneCount <- count_gene_list(path_plot$Down_regulated)
path_plot$TotalGeneCount <- path_plot$UpGeneCount + path_plot$DownGeneCount
path_plot$MinusLog10P <- -log10(pmax(path_plot$lowest_p, .Machine$double.xmin))
write_csv(path_plot, "13_pathfindR_Reactome_plot_data.csv")

top_path <- head(path_plot[order(path_plot$lowest_p), ], 20)
top_path$TermLabel <- vapply(top_path$Term_Description, function(x) {
  paste(strwrap(x, width = 55), collapse = "\n")
}, character(1))

p_path <- ggplot(top_path, aes(Fold_Enrichment, reorder(TermLabel, Fold_Enrichment))) +
  geom_point(aes(color = MinusLog10P, size = TotalGeneCount), alpha = 0.85) +
  scale_color_gradient(low = "#2C7FB8", high = "#D95F02") +
  scale_size_continuous(range = c(2.5, 9)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 8), legend.position = "right", plot.title = element_text(face = "bold")) +
  labs(
    title = "pathfindR Reactome active-subnetwork enrichment",
    x = "Fold enrichment",
    y = "",
    color = "-log10 P",
    size = "Genes"
  )

open_png("23_pathfindR_Reactome_dotplot.png", width = 1710, height = 1260)
print(p_path)
close_png()

open_png("23_pathfindR_Reactome_dotplot_clean.png", width = 1710, height = 1260)
print(p_path)
close_png()

############################################
# Obtaining STRING and co-expression network
############################################

string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400, input_directory = object_dir)
string_input <- data.frame(GeneSymbol = sig_symbols)
string_mapped <- string_db$map(string_input, "GeneSymbol", removeUnmappedRows = TRUE)
string_interactions <- string_db$get_interactions(string_mapped$STRING_id)
write_csv(string_mapped, "14_STRING_mapped_genes.csv")
write_csv(string_interactions, "14_STRING_interactions.csv")

network_genes <- head(intersect(primary_deg$GeneSymbol, rownames(expr_gene)), 220)
expr_network <- expr_gene[network_genes, , drop = FALSE]

gene_cor <- cor(t(expr_network), method = "pearson")
n_samples <- ncol(expr_network)
t_stat <- gene_cor * sqrt((n_samples - 2) / pmax(1 - gene_cor^2, .Machine$double.eps))
p_mat <- 2 * pt(-abs(t_stat), df = n_samples - 2)
diag(p_mat) <- 1

upper_idx <- which(upper.tri(gene_cor), arr.ind = TRUE)
edge_pool <- data.frame(
  From = rownames(gene_cor)[upper_idx[, 1]],
  To = colnames(gene_cor)[upper_idx[, 2]],
  Correlation = gene_cor[upper_idx],
  P.Value = p_mat[upper_idx]
)
edge_pool$adj.P.Val <- p.adjust(edge_pool$P.Value, method = "BH")

selected_edges <- NULL
for (cutoff in c(0.85, 0.80, 0.75, 0.70, 0.65)) {
  candidate <- subset(edge_pool, abs(Correlation) >= cutoff & adj.P.Val < 0.05)
  if (nrow(candidate) >= 30) {
    selected_edges <- candidate
    break
  }
}

if (is.null(selected_edges)) {
  selected_edges <- head(edge_pool[order(edge_pool$adj.P.Val, -abs(edge_pool$Correlation)), ], 30)
}

selected_edges <- head(selected_edges[order(-abs(selected_edges$Correlation)), ], 1500)

network_status <- tumor_vs_normal[match(network_genes, tumor_vs_normal$GeneSymbol), ]
nodes <- data.frame(
  name = network_genes,
  GeneSymbol = network_genes,
  logFC = network_status$logFC,
  adj.P.Val = network_status$adj.P.Val,
  Status = network_status$Status
)
nodes$Status[is.na(nodes$Status)] <- "Not significant"

coexpression_graph <- graph_from_data_frame(selected_edges, directed = FALSE, vertices = nodes)
E(coexpression_graph)$Weight <- abs(E(coexpression_graph)$Correlation)
E(coexpression_graph)$Sign <- ifelse(E(coexpression_graph)$Correlation >= 0, "Positive", "Negative")

centrality <- data.frame(
  GeneSymbol = V(coexpression_graph)$name,
  Degree = degree(coexpression_graph),
  Strength = strength(coexpression_graph, weights = E(coexpression_graph)$Weight),
  Closeness = closeness(coexpression_graph, weights = 1 / E(coexpression_graph)$Weight, normalized = TRUE),
  Betweenness = betweenness(coexpression_graph, weights = 1 / E(coexpression_graph)$Weight, normalized = TRUE),
  EigenCentrality = eigen_centrality(coexpression_graph, weights = E(coexpression_graph)$Weight)$vector
)

community <- cluster_louvain(coexpression_graph, weights = E(coexpression_graph)$Weight)
centrality$Community <- membership(community)[centrality$GeneSymbol]
V(coexpression_graph)$Community <- centrality$Community[match(V(coexpression_graph)$name, centrality$GeneSymbol)]

centrality <- merge(centrality, nodes, by = "GeneSymbol", all.x = TRUE)
centrality <- centrality[order(centrality$Degree, decreasing = TRUE), ]

write_csv(centrality, "15_coexpression_network_centrality.csv")
write_csv(selected_edges, "15_coexpression_network_edges.csv")
write_csv(nodes, "15_coexpression_network_nodes.csv")

write.csv(centrality, file.path(cytoscape_dir, "coexpression_nodes.csv"), row.names = FALSE)
write.csv(selected_edges, file.path(cytoscape_dir, "coexpression_edges.csv"), row.names = FALSE)
write_graph(coexpression_graph, file.path(cytoscape_dir, "coexpression_network.graphml"), format = "graphml")

set.seed(123)
layout_fr <- layout_with_fr(coexpression_graph, weights = E(coexpression_graph)$Weight)
node_colors <- c("Up" = "#D95F02", "Down" = "#2C7FB8", "Not significant" = "grey70")
node_size <- 4 + 10 * degree(coexpression_graph) / max(degree(coexpression_graph))

open_png("24_coexpression_network.png", width = 1700, height = 1500)
plot(
  coexpression_graph,
  layout = layout_fr,
  vertex.label = NA,
  vertex.size = node_size,
  vertex.color = node_colors[V(coexpression_graph)$Status],
  edge.width = 0.5 + 2 * E(coexpression_graph)$Weight,
  edge.color = ifelse(E(coexpression_graph)$Sign == "Positive", "#77777780", "#D95F0280"),
  main = "GSE50161 co-expression network"
)
legend("topleft", legend = names(node_colors), col = node_colors, pch = 19, bty = "n", cex = 0.8)
close_png()

try({
  createNetworkFromIgraph(
    coexpression_graph,
    title = paste0(GEO_ACCESSION, " coexpression network"),
    collection = GEO_ACCESSION
  )
}, silent = TRUE)

#####################################
# Performing a brief analysis summary
#####################################

summary_lines <- c(
  paste0("# ", GEO_ACCESSION, " dataset explanation and pipeline summary"),
  "",
  "## Dataset",
  paste0(
    GEO_ACCESSION,
    " contains Affymetrix HG-U133 Plus 2.0 gene-expression profiles from ependymoma, ",
    "glioblastoma, medulloblastoma, pilocytic astrocytoma, and normal brain tissue."
  ),
  paste0("Processed probe matrix: ", nrow(expr_processed), " probes x ", ncol(expr_processed), " samples."),
  paste0("Gene-level matrix: ", nrow(expr_gene), " genes x ", ncol(expr_gene), " samples."),
  "",
  "## Sample division",
  paste(apply(sample_counts, 1, function(x) paste0("- ", x[["Group"]], ": ", x[["N"]], " samples")), collapse = "\n"),
  "",
  "## Preprocessing",
  paste0("- Log2 transform applied: ", log2_applied),
  paste0("- Between-array normalization applied: ", normalization_applied),
  "- One representative probe per gene was retained using highest probe variance.",
  "",
  "## Main analyses",
  "- limma differential expression",
  "- PCA, k-means, hierarchical clustering",
  "- Random Forest, LDA, LASSO, Ridge, caret repeated CV",
  "- rank-based signatures and rScudo",
  "- gProfiler, clusterProfiler, ReactomePA, pathfindR, STRING",
  "- co-expression network centrality and Cytoscape exports"
)

writeLines(summary_lines, file.path(results_dir, "GSE50161_dataset_explanation_and_pipeline_summary.md"))
capture.output(sessionInfo(), file = file.path(results_dir, "sessionInfo.txt"))