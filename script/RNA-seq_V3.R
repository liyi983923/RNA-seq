
cat("==================================================\n")
cat("  RNA-seq 一键分析 V8.0 启动 - 优雅完整版\n")
cat("  开始时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==================================================\n")

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(DESeq2)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(pheatmap)
  library(biomaRt)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(showtext)
  library(GSEABase)

})

select <- dplyr::select
filter <- dplyr::filter

# ------------------- 全局美学设置 -------------------
font_add_google("Roboto Condensed", "roboto")
showtext_auto()

THEME_PUB   <- theme_pubr(base_size = 16, base_family = "roboto")
THEME_MIN   <- theme_minimal(base_size = 16, base_family = "roboto")
DPI         <- 350
WIDTH_STD   <- 11
HEIGHT_STD  <- 8

# 颜色定义（发表级）
GROUP_COLORS <- c(
  "WT_Control"     = "#377eb8",   # 深蓝
  "KI_Control"     = "#1b9e77",   # 深绿
  "Model"          = "#d95f02",   # 橘红
  "Upada_10mg"     = "#7570b3",   # 紫色
  "Cmpd261_5mg"    = "#e7298a",   # 洋红
  "Cmpd261_10mg"   = "#e41a1c",   # 大红
  "EV756_10mg"     = "#66a61e",   # 草绿
  "Cmpd261_5mg_B6" = "#e6ab02"    # 金黄
)
TISSUE_COLORS <- c("Colon" = "#ff9999", "Muscle" = "#99cc99", "Skin" = "#9999ff")

# ------------------- 1. 读取数据 -------------------
cat("1/9 读取计数矩阵...\n")
count_file <- here("data", "id.txt")
raw <- read_tsv(count_file, comment = "#", skip = 1, show_col_types = FALSE)

count_matrix <- raw %>%
  select(Geneid, starts_with("G")) %>%
  group_by(Geneid) %>%
  summarise(across(starts_with("G"), sum, na.rm = TRUE)) %>%
  ungroup() %>%
  rename_with(~ str_remove(., "_hits\\.bam$"), starts_with("G")) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

cat("   计数矩阵：", nrow(count_matrix), "genes ×", ncol(count_matrix), "samples\n")

# ------------------- 2. 样本信息表 -------------------
cat("2/9 构建样本信息...\n")
col_data <- tibble(sample = colnames(count_matrix)) %>%
  mutate(
    tissue = case_when(
      str_detect(sample, "-JC-") ~ "Colon",
      str_detect(sample, "-JR-") ~ "Muscle",
      str_detect(sample, "-PF-") ~ "Skin",
      TRUE ~ NA_character_
    ),
    group = case_when(
      str_detect(sample, "^G1-1-") ~ "WT_Control",
      str_detect(sample, "^G1-2-") ~ "KI_Control",
      str_detect(sample, "^G2-")    ~ "Model",
      str_detect(sample, "^G3-")    ~ "Upada_10mg",
      str_detect(sample, "^G4-")    ~ "Cmpd261_5mg",
      str_detect(sample, "^G5-")    ~ "Cmpd261_10mg",
      str_detect(sample, "^G6-")    ~ "EV756_10mg",
      str_detect(sample, "^G7-")    ~ "Cmpd261_5mg_B6",
      TRUE ~ "Unknown"
    )
  ) %>%
  mutate(
    tissue = factor(tissue, levels = c("Colon", "Muscle", "Skin")),
    group  = factor(group, levels = names(GROUP_COLORS))
  ) %>%
  column_to_rownames("sample")

# ------------------- 3. 基因注释缓存（只查一次） -------------------
cat("3/9 加载基因注释（一次性缓存）...\n")
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
all_genes <- rownames(count_matrix)
gene_anno <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = all_genes,
  mart = mart
) %>% as_tibble() %>%
  mutate(
    gene_name = coalesce(mgi_symbol, external_gene_name, ensembl_gene_id),
    description = str_remove(description, " \\[Source.*\\]$")
  )

# ------------------- 4. 所有18个对比定义 -------------------
all_contrasts <- tribble(
  ~tissue, ~test, ~ref, ~folder, ~name,
  "Skin",   "Model",          "WT_Control", "Model_vs_WT",          "Model_vs_WT_Skin",
  "Skin",   "Cmpd261_5mg",    "WT_Control", "Cmpd261_5mg_vs_WT",    "Cmpd261_5mg_vs_WT_Skin",
  "Skin",   "Cmpd261_10mg",   "WT_Control", "Cmpd261_10mg_vs_WT",   "Cmpd261_10mg_vs_WT_Skin",
  "Skin",   "EV756_10mg",     "WT_Control", "EV756_10mg_vs_WT",     "EV756_10mg_vs_WT_Skin",
  "Skin",   "Upada_10mg",     "WT_Control", "Upada_10mg_vs_WT",     "Upada_10mg_vs_WT_Skin",
  "Skin",   "Cmpd261_5mg_B6", "WT_Control", "Cmpd261_5mg_B6_vs_WT", "Cmpd261_5mg_B6_vs_WT_Skin",
  
  "Colon",  "Model",          "WT_Control", "Model_vs_WT",          "Model_vs_WT_Colon",
  "Colon",  "Cmpd261_5mg",    "WT_Control", "Cmpd261_5mg_vs_WT",    "Cmpd261_5mg_vs_WT_Colon",
  "Colon",  "Cmpd261_10mg",   "WT_Control", "Cmpd261_10mg_vs_WT",   "Cmpd261_10mg_vs_WT_Colon",
  "Colon",  "EV756_10mg",     "WT_Control", "EV756_10mg_vs_WT",     "EV756_10mg_vs_WT_Colon",
  "Colon",  "Upada_10mg",     "WT_Control", "Upada_10mg_vs_WT",     "Upada_10mg_vs_WT_Colon",
  "Colon",  "Cmpd261_5mg_B6", "WT_Control", "Cmpd261_5mg_B6_vs_WT", "Cmpd261_5mg_B6_vs_WT_Colon",
  
  "Muscle", "Model",          "WT_Control", "Model_vs_WT",          "Model_vs_WT_Muscle",
  "Muscle", "Cmpd261_5mg",    "WT_Control", "Cmpd261_5mg_vs_WT",    "Cmpd261_5mg_vs_WT_Muscle",
  "Muscle", "Cmpd261_10mg",   "WT_Control", "Cmpd261_10mg_vs_WT",   "Cmpd261_10mg_vs_WT_Muscle",
  "Muscle", "EV756_10mg",     "WT_Control", "EV756_10mg_vs_WT",     "EV756_10mg_vs_WT_Muscle",
  "Muscle", "Upada_10mg",     "WT_Control", "Upada_10mg_vs_WT",     "Upada_10mg_vs_WT_Muscle",
  "Muscle", "Cmpd261_5mg_B6", "WT_Control", "Cmpd261_5mg_B6_vs_WT", "Cmpd261_5mg_B6_vs_WT_Muscle"
)

# ------------------- 5. 主分析函数 -------------------
run_deg_analysis <- function(i, contrast) {
  cat(sprintf("\n5/9 [%2d/18] %s\n", i, contrast$name))
  
  out_dir <- here("results_final", contrast$tissue, contrast$folder)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  samples <- rownames(col_data)[col_data$tissue == contrast$tissue &
                                  col_data$group %in% c(contrast$test, contrast$ref)]
  if (length(samples) < 4) { cat("   样本不足，跳过\n"); return(NULL) }
  
  coldata_sub <- col_data[samples, ] %>% mutate(group2 = factor(group, levels = c(contrast$ref, contrast$test)))
  counts_sub <- count_matrix[, samples, drop = FALSE]
  keep <- rowSums(counts_sub >= 10) >= 3
  counts_filt <- counts_sub[keep, ]
  
  dds <- DESeqDataSetFromMatrix(counts_filt, coldata_sub, ~group2) %>% DESeq(quiet = TRUE)
  vsd <- vst(dds, blind = FALSE)
  
  res <- lfcShrink(dds, coef = 2, type = "apeglm") %>%
    as.data.frame() %>% rownames_to_column("Geneid") %>% as_tibble() %>%
    left_join(gene_anno %>% select(Geneid = ensembl_gene_id, gene_name, description), by = "Geneid") %>%
    mutate(gene_name = ifelse(is.na(gene_name), Geneid, gene_name),
           padj = ifelse(is.na(padj), 1, padj),
           significant = padj < 0.05 & abs(log2FoldChange) > 1,
           regulation = case_when(padj < 0.05 & log2FoldChange > 1 ~ "Up",
                                  padj < 0.05 & log2FoldChange < -1 ~ "Down", TRUE ~ "NS"))
  
  write_csv(res, file.path(out_dir, "DEG_shrunken_full.csv"))
  write_csv(counts(dds, normalized = TRUE) %>% as.data.frame() %>% rownames_to_column("Geneid"),
            file.path(out_dir, "normalized_counts.csv"))
  
  # PCA
  pca_data <- plotPCA(vsd, intgroup = "group2", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group2)) +
    geom_point(size = 6) + geom_text_repel(aes(label = name), size = 5, family = "roboto") +
    scale_color_manual(values = GROUP_COLORS) +
    labs(title = paste("PCA:", contrast$name),
         x = paste0("PC1: ", percentVar[1], "%"), y = paste0("PC2: ", percentVar[2], "%")) +
    THEME_PUB
  ggsave(file.path(out_dir, "1_PCA.png"), p_pca, width = 9, height = 7, dpi = DPI)
  
  # Volcano
  top12 <- res %>% filter(significant) %>% arrange(padj) %>% head(12)
  p_vol <- ggplot(res, aes(log2FoldChange, -log10(padj), color = regulation)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_color_manual(values = c("Up" = "#e31a1c", "Down" = "#1f78b4", "NS" = "grey70")) +
    geom_vline(xintercept = c(-1,1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(data = top12, aes(label = gene_name), size = 5, max.overlaps = 20) +
    labs(title = paste("Volcano:", contrast$name),
         subtitle = paste(sum(res$significant), "DEGs (Up", sum(res$regulation=="Up"), "| Down", sum(res$regulation=="Down"), ")")) +
    THEME_MIN
  ggsave(file.path(out_dir, "2_Volcano_labeled.png"), p_vol, width = 11, height = 8.5, dpi = DPI)
  
  # Top50 Heatmap
  top50 <- res %>% filter(significant) %>% arrange(padj) %>% head(50)
  if (nrow(top50) > 0) {
    mat <- assay(vsd)[top50$Geneid, ]
    rownames(mat) <- top50$gene_name
    anno_col <- data.frame(Group = coldata_sub$group2, row.names = colnames(mat))
    pheatmap(mat, scale = "row", annotation_col = anno_col,
             color = colorRampPalette(c("#377eb8","white","#e31a1c"))(100),
             show_rownames = TRUE, fontsize_row = 8,
             filename = file.path(out_dir, "3_Heatmap_top50.png"),
             width = 10, height = 14, res = DPI)
  }
  
  # 经典炎症基因
  inflam_genes <- c("Il1b","Il6","Tnf","Cxcl2","S100a8","S100a9","Defb3",
                    "Krt16","Krt6a","Lcn2","Mmp9","Ptgs2","Cxcl1","Cxcl10","Il17a","Il17f")
  for (g in inflam_genes) {
    row <- res %>% filter(str_detect(gene_name, paste0("^", g, "$")))
    if (nrow(row) == 0) next
    gene_id <- row$Geneid[1]
    plot_data <- plotCounts(dds, gene = gene_id, intgroup = "group2", normalized = TRUE, returnData = TRUE)
    plot_data$group2 <- coldata_sub$group2[match(rownames(plot_data), rownames(coldata_sub))]
    p <- ggplot(plot_data, aes(group2, count + 1, color = group2)) +
      geom_boxplot(width = 0.4, outlier.shape = NA) + geom_jitter(width = 0.25, size = 6) +
      scale_y_log10() + scale_color_manual(values = GROUP_COLORS) +
      labs(title = paste0(g, " expression"), y = "Normalized count + 1",
           subtitle = paste("padj =", signif(row$padj[1], 3))) +
      THEME_PUB + theme(legend.position = "none", axis.title.x = element_blank())
    ggsave(file.path(out_dir, paste0("Gene_", g, ".png")), p, width = 5.5, height = 6, dpi = DPI, bg = "white")
  }
  
  return(res %>% mutate(Comparison = contrast$name))
}

cat("5/9 开始18个差异分析...\n")
all_deg_list <- vector("list", nrow(all_contrasts))
for (i in 1:nrow(all_contrasts)) {
  all_deg_list[[i]] <- run_deg_analysis(i, all_contrasts[i, ])
}
all_deg_combined <- bind_rows(all_deg_list)
write_csv(all_deg_combined, "results_final/ALL_18_DEG_COMBINED.csv")
write_csv(all_deg_combined %>% left_join(gene_anno, by = c("Geneid" = "ensembl_gene_id")),
          "results_final/ALL_18_DEG_COMBINED_WITH_ANNOTATION.csv")
cat("   18个差异分析完成！\n")

# ------------------- 6. 全样本PCA -------------------
cat("6/9 绘制全样本PCA...\n")

# 主标题样式# 只保留7个核心组
valid_groups_for_pca <- c("WT_Control", "Model", 
                          "Cmpd261_5mg", "Cmpd261_10mg", 
                          "EV756_10mg", "Upada_10mg", "Cmpd261_5mg_B6")

keep_samples_pca <- rownames(col_data)[col_data$group %in% valid_groups_for_pca]
count_matrix_pca <- count_matrix[, keep_samples_pca, drop = FALSE]
col_data_pca <- col_data[keep_samples_pca, , drop = FALSE]

# 构建vsd
keep_genes_global <- rowSums(count_matrix_pca >= 10) >= 6
dds_full <- DESeqDataSetFromMatrix(count_matrix_pca[keep_genes_global, ], col_data_pca, ~ group + tissue)
dds_full <- DESeq(dds_full, quiet = TRUE)
vsd_full <- vst(dds_full, blind = TRUE)

# 手动PCA
mat <- assay(vsd_full)
pca <- prcomp(t(mat), center = TRUE, scale. = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], sample = rownames(pca$x)) %>%
  left_join(col_data_pca %>% rownames_to_column("sample"), by = "sample")

# 颜色定义
group_colors_pca <- c("WT_Control"="#377eb8","Model"="#d95f02",
                      "Cmpd261_5mg"="#e7298a","Cmpd261_10mg"="#e41a1c",
                      "EV756_10mg"="#66a61e","Upada_10mg"="#7570b3",
                      "Cmpd261_5mg_B6"="#e6ab02")

# 最终出图（所有图例参数都暴露给你调！）
p_global_pca_tunable <- ggplot(pca_data, aes(PC1, PC2, color = group, shape = tissue)) +
  geom_point(size = 8, stroke = 0, alpha = 1) +
  scale_color_manual(values = group_colors_pca, name = "Treatment Group") +
  scale_shape_manual(values = c("Colon"=19, "Muscle"=17, "Skin"=15), name = "Tissue") +
  labs(title = "Global PCA - Treatment Groups",
       subtitle = "Colored by Treatment Group | Shaped by Tissue",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  
  # —————————————————————— 图例超级可调区（你随便改！）——————————————————————
  theme_pubr(base_size = 22, base_family = "roboto") +
  theme(
    # 图例整体位置：right, left, top, bottom, none 或 c(0.8, 0.2) 坐标
    legend.position = "right",           # 改这里换位置
    
    # 图例背景、边框、间距
    # legend.background = element_rect(fill = "white", color = "grey80", size = 0.5),
    # le
    plot.title = element_text(face = "bold", size = 28, hjust = 0.5),
    plot.subtitle = element_text(size = 18, hjust = 0.5, color = "black"),
    
    panel.grid.minor = element_blank()
  )
# —————————————————————— 图例可调区结束 ——————————————————————

# 显示图形（RStudio里直接出图）
p_global_pca_tunable

dev.off()

# 保存高清图
ggsave("results_final/Global_PCA_Clean_No_KI_Control_TUNABLE.png", 
       p_global_pca_tunable, width = 16, height = 10, dpi = 350, bg = "white")


cat("   全样本PCA完成！\n")

# ------------------- 7. 逆转分析（12个） -------------------
cat("7/9 开始12个逆转分析...\n")
reverse_contrasts <- tribble(
  ~tissue, ~test, ~ref, ~folder, ~display,
  "Skin",   "Cmpd261_5mg",    "Model", "Cmpd261_5mg_vs_Model", "261 5mg vs Model (Skin)",
  "Skin",   "Cmpd261_10mg",   "Model", "Cmpd261_10mg_vs_Model", "261 10mg vs Model (Skin)",
  "Skin",   "EV756_10mg",     "Model", "EV756_10mg_vs_Model",   "EV756 vs Model (Skin)",
  "Skin",   "Upada_10mg",     "Model", "Upada_10mg_vs_Model",   "Upadacitinib vs Model (Skin)",
  "Colon",  "Cmpd261_5mg",    "Model", "Cmpd261_5mg_vs_Model", "261 5mg vs Model (Colon)",
  "Colon",  "Cmpd261_10mg",   "Model", "Cmpd261_10mg_vs_Model", "261 10mg vs Model (Colon)",
  "Colon",  "EV756_10mg",     "Model", "EV756_10mg_vs_Model",   "EV756 vs Model (Colon)",
  "Colon",  "Upada_10mg",     "Model", "Upada_10mg_vs_Model",   "Upadacitinib vs Model (Colon)",
  "Muscle", "Cmpd261_5mg",    "Model", "Cmpd261_5mg_vs_Model", "261 5mg vs Model (Muscle)",
  "Muscle", "Cmpd261_10mg",   "Model", "Cmpd261_10mg_vs_Model", "261 10mg vs Model (Muscle)",
  "Muscle", "EV756_10mg",     "Model", "EV756_10mg_vs_Model",   "EV756 vs Model (Muscle)",
  "Muscle", "Upada_10mg",     "Model", "Upada_10mg_vs_Model",   "Upadacitinib vs Model (Muscle)"
)

hallmark_gmt <- here("data", "mh.all.v2025.1.Mm.symbols.gmt")
hallmark_sets <- if (file.exists(hallmark_gmt)) getGmt(hallmark_gmt) else NULL

dir.create("results_final/Reversal_Analysis", recursive = TRUE, showWarnings = FALSE)

for (i in 1:nrow(reverse_contrasts)) {
  rc <- reverse_contrasts[i, ]
  cat(sprintf("   [%2d/12] %s\n", i, rc$display))
  
  out_dir <- here("results_final/Reversal_Analysis", rc$tissue, rc$folder)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  samples <- rownames(col_data)[col_data$tissue == rc$tissue & col_data$group %in% c("Model", rc$test)]
  if (length(samples) < 4) next
  
  coldata_sub <- col_data[samples, ] %>% mutate(group2 = factor(group, levels = c("Model", rc$test)))
  counts_sub <- count_matrix[, samples, drop = FALSE]
  keep <- rowSums(counts_sub >= 10) >= 3
  counts_filt <- counts_sub[keep, ]
  
  dds <- DESeqDataSetFromMatrix(counts_filt, coldata_sub, ~group2) %>% DESeq(quiet = TRUE)
  vsd <- vst(dds, blind = FALSE)
  
  res <- lfcShrink(dds, coef = 2, type = "apeglm") %>%
    as.data.frame() %>% rownames_to_column("Geneid") %>% as_tibble() %>%
    left_join(gene_anno %>% select(Geneid = ensembl_gene_id, symbol = gene_name), by = "Geneid") %>%
    mutate(final_name = ifelse(is.na(symbol) | symbol == "", Geneid, symbol))
  
  write_csv(res, file.path(out_dir, "DEG_vs_Model.csv"))
  
  # 逆转热图
  model_file <- here("results_final", rc$tissue, "Model_vs_WT", "DEG_shrunken_full.csv")
  if (file.exists(model_file)) {
    model_sig <- read_csv(model_file, show_col_types = FALSE) %>% filter(significant)
    rev_genes <- res %>%
      filter(Geneid %in% model_sig$Geneid,
             sign(log2FoldChange) != sign(model_sig$log2FoldChange[match(Geneid, model_sig$Geneid)])) %>%
      arrange(padj) %>% head(50)
    if (nrow(rev_genes) >= 5) {
      mat <- assay(vsd)[rev_genes$Geneid, ]
      rownames(mat) <- rev_genes$final_name
      anno_col <- data.frame(Group = coldata_sub$group2, row.names = colnames(mat))
      pheatmap(mat, scale = "row", annotation_col = anno_col,
               color = colorRampPalette(c("#1f78b4", "white", "#e31a1c"))(100),
               filename = file.path(out_dir, "Reversal_Heatmap.png"),
               width = 9, height = 11, res = DPI)
    }
  }
  
  # GSEA
  if (!is.null(hallmark_sets)) {
    ranked <- setNames(res$log2FoldChange, res$final_name)
    ranked <- ranked[!is.na(names(ranked)) & names(ranked) != ""]
    ranked <- ranked[!duplicated(names(ranked))]
    ranked <- sort(ranked, decreasing = TRUE)
    if (length(ranked) > 50) {
      gsea <- tryCatch(GSEA(ranked, TERM2GENE = hallmark_t2g, pvalueCutoff = 0.25), error = function(e) NULL)
      if (!is.null(gsea) && nrow(gsea) > 0) {
        dotplot(gsea, showCategory = 20) + ggtitle(rc$display) %>% ggsave(file.path(out_dir, "GSEA_Dotplot.png"), ., width = 11, height = 8, dpi = DPI)
      }
    }
  }
}

# ------------------- 8. 终极逆转联合热图（3张） -------------------


# ==============================================================================
cat("8/9 正在生成3张终极发表级逆转热图（终极修复版）...\n")

options(bitmapType = "quartz")
library(ComplexHeatmap)
library(circlize)
library(grid)
# === 8.1. 收集逆转基因（不变）===
rev_files <- list.files("results_final/Reversal_Analysis", pattern = "DEG_vs_Model.csv", 
                        recursive = TRUE, full.names = TRUE)

all_rev_genes <- NULL
for (f in rev_files) {
  tissue <- strsplit(f, "/")[[1]][3]
  drug_res <- read_csv(f, show_col_types = FALSE)
  model_file <- here("results_final", tissue, "Model_vs_WT", "DEG_shrunken_full.csv")
  if (!file.exists(model_file)) next
  model_res <- read_csv(model_file, show_col_types = FALSE)
  model_sig <- subset(model_res, significant == TRUE)
  common <- intersect(drug_res$Geneid, model_sig$Geneid)
  if (length(common) == 0) next
  drug_sub <- subset(drug_res, Geneid %in% common)
  model_sub <- subset(model_sig, Geneid %in% common)
  opposite <- sign(drug_sub$log2FoldChange) != sign(model_sub$log2FoldChange) &
    sign(drug_sub$log2FoldChange) != 0
  if (sum(opposite) == 0) next
  temp <- drug_sub[opposite, ][order(drug_sub$padj[opposite]), ][1:min(50, sum(opposite)), ]
  all_rev_genes <- if (is.null(all_rev_genes)) temp else rbind(all_rev_genes, temp)
}
if (is.null(all_rev_genes) || nrow(all_rev_genes) == 0) stop("无逆转基因！")

reversal_genes <- all_rev_genes[order(all_rev_genes$padj), ][!duplicated(all_rev_genes$Geneid), ]

# === 8.2. 基因名 ===
symbol_table <- data.frame(
  Geneid = gene_anno$ensembl_gene_id,
  symbol = ifelse(!is.na(gene_anno$mgi_symbol) & gene_anno$mgi_symbol != "", gene_anno$mgi_symbol,
                  ifelse(!is.na(gene_anno$external_gene_name) & gene_anno$external_gene_name != "",
                         gene_anno$external_gene_name, gene_anno$ensembl_gene_id)),
  stringsAsFactors = FALSE
)
gene_labels <- symbol_table$symbol[match(reversal_genes$Geneid, symbol_table$Geneid)]
gene_labels[is.na(gene_labels)] <- reversal_genes$Geneid[is.na(gene_labels)]
gene_labels <- make.unique(gene_labels)

# === 8.3. 表达矩阵 + Z-score ===
expr_raw <- assay(vsd_full)[reversal_genes$Geneid, ]
rownames(expr_raw) <- gene_labels
expr_z <- t(scale(t(expr_raw)))

# === 8.4. 严格控制样本顺序（关键修复！）===
desired_order <- c("WT_Control", "Model", "Upada_10mg", 
                   "Cmpd261_5mg", "Cmpd261_10mg", "EV756_10mg", "Cmpd261_5mg_B6")

# 只保留这7个组的样本
valid_samples <- rownames(col_data)[col_data$group %in% desired_order]
expr_z <- expr_z[, valid_samples, drop = FALSE]

# 创建注释数据框（顺序严格）
anno_df <- data.frame(
  Group = factor(col_data[valid_samples, "group"], levels = desired_order),
  Tissue = col_data[valid_samples, "tissue"],
  row.names = valid_samples
)

# === 8.5. 颜色（只保留 Group）===
group_colors <- c(
  "WT_Control"     = "#377eb8",
  "Model"          = "#d95f02",
  "Upada_10mg"     = "#7570b3",   # 第3位
  "Cmpd261_5mg"    = "#e7298a",
  "Cmpd261_10mg"   = "#e41a1c",
  "EV756_10mg"     = "#66a61e",
  "Cmpd261_5mg_B6" = "#e6ab02"
)

ha <- HeatmapAnnotation(
  Group = anno_df$Group,
  col = list(Group = group_colors),
  annotation_name_gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "roboto"),
  annotation_legend_param = list(title_gp = gpar(fontsize = 16, fontface = "bold"))
)

# === 8.6. 分组织画图（彻底修复索引问题）===
# Skin 的完整代码====================
tissue <- "Skin"

samples_to_plot <- rownames(col_data)[col_data$tissue == tissue & col_data$group %in% desired_order]
mat <- expr_z[, samples_to_plot, drop = FALSE]
mat <- mat[, order(factor(col_data[samples_to_plot, "group"], levels = desired_order)), drop = FALSE]

ha <- columnAnnotation(
  Group = factor(col_data[samples_to_plot, "group"], levels = desired_order),
  col = list(Group = group_colors),
  annotation_name_gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "roboto"),
  show_legend = FALSE                     # 关键！关闭自动图例
)

n_genes <- nrow(mat)
col_fun <- colorRamp2(c(-2, 0, 2), c("#436eee", "white", "#ee4444"))

# Heatmap 主代码（其他参数不变）
# 1. 原来的 Heatmap 对象保持不变
ht <- Heatmap(mat,
              name = "Z-score",
              col = col_fun,
              bottom_annotation = ha,
              show_heatmap_legend = FALSE,    # 加上这行，热图自己就不画颜色条了
              
              show_row_names = TRUE,
              row_names_side = "right",
              row_names_gp = gpar(fontsize = 2, fontfamily = "roboto"),
              
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 16, fontfamily = "roboto", fontface ="bold"),
              column_names_rot = 45,
              column_names_side = "bottom",
              
              column_title = paste0("Top Differentially Expressed Genes in Skin (n = ", n_genes, " genes)"),
              column_title_gp = gpar(fontsize = 26, fontface = "bold", fontfamily = "roboto"),
              
              row_gap = unit(0.5, "mm"),
              column_gap = unit(2, "mm"),
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              use_raster = FALSE,
              width = unit(16, "cm"),
              height = unit(pmax(12, n_genes * 0.07), "cm"),
              
              heatmap_legend_param = list(
                title = "Expression\nZ-score",
                title_gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "roboto"),
                labels_gp = gpar(fontsize = 16, fontfamily = "roboto")
              )
)

# 2. 手动创建 Treatment Group 图例
treatment_legend <- Legend(
  at = desired_order,
  title = "Treatment Group",
  legend_gp = gpar(fill = group_colors[desired_order]),
  title_gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "roboto"),
  labels_gp = gpar(fontsize = 16, fontfamily = "roboto"),
  direction = "vertical"   # 确保垂直
)

expr_legend <- Legend(
  col_fun = col_fun,                 # 关键：这里直接用你定义的 col_fun
  title = "Expression\nZ-score",
  title_gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "roboto"),
  labels_gp = gpar(fontsize = 16, fontfamily = "roboto"),
  direction = "vertical"
)

# Treatment Group 图例（你原来就有的）
treatment_legend <- Legend(
  at = desired_order,
  title = "Treatment Group",
  legend_gp = gpar(fill = group_colors[desired_order]),
  title_gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "roboto"),
  labels_gp = gpar(fontsize = 16, fontfamily = "roboto"),
  direction = "vertical"
)

# 现在 packLegend 绝对不会报错了
combined_legend <- packLegend(
  expr_legend,
  treatment_legend,
  direction = "vertical",
  gap = unit(12, "mm")      # 两个图例之间的间距，自己调
)

# 5. 保存图片
pdf("/Users/liyi/Downloads/RNA_Seq/X2_mouse/results_final/Reversal_Analysis/REVERSAL_HEATMAP_Skin_FINAL.pdf",
    width = 24, 
    height = pmax(14, n_genes * 0.08) + 8)

draw(ht,
     heatmap_legend_side = "right",      # 热图图例本来也会放右边，但我们现在不用它了
     annotation_legend_side = "right",
     annotation_legend_list = list(combined_legend),  # 只放这一个打包好的图例
     padding = unit(c(5, 35, 5, 5), "mm")  # 右边多留点空间给两个图例
)

dev.off()


# =============================================
# Muscle 的完整代码
# =============================================
tissue <- "Muscle"

samples_to_plot <- rownames(col_data)[col_data$tissue == tissue & col_data$group %in% desired_order]
mat <- expr_z[, samples_to_plot, drop = FALSE]
mat <- mat[, order(factor(col_data[samples_to_plot, "group"], levels = desired_order)), drop = FALSE]

ha <- columnAnnotation(
  Group = factor(col_data[samples_to_plot, "group"], levels = desired_order),
  col = list(Group = group_colors),
  annotation_name_gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "roboto"),
  show_legend = FALSE
)

n_genes <- nrow(mat)
col_fun <- colorRamp2(c(-2, 0, 2), c("#436eee", "white", "#ee4444"))

ht <- Heatmap(mat,
              name = "Z-score",
              col = col_fun,
              bottom_annotation = ha,
              show_heatmap_legend = FALSE,
              
              show_row_names = TRUE,
              row_names_side = "right",
              row_names_gp = gpar(fontsize = 2, fontfamily = "roboto"),
              
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 16, fontfamily = "roboto", fontface = "bold"),
              column_names_rot = 45,
              
              column_title = paste0("Top Differentially Expressed Genes in Muscle (n = ", n_genes, " genes)"),
              column_title_gp = gpar(fontsize = 26, fontface = "bold", fontfamily = "roboto"),
              
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              width = unit(16, "cm"),
              height = unit(pmax(12, n_genes * 0.07), "cm"),
              row_gap = unit(0.5, "mm"),
              column_gap = unit(2, "mm"),
              use_raster = FALSE
)

expr_legend <- Legend(col_fun = col_fun, title = "Expression\nZ-score",
                      title_gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "roboto"),
                      labels_gp = gpar(fontsize = 16, fontfamily = "roboto"))

treatment_legend <- Legend(at = desired_order, title = "Treatment Group",
                           legend_gp = gpar(fill = group_colors[desired_order]),
                           title_gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "roboto"),
                           labels_gp = gpar(fontsize = 16, fontfamily = "roboto"))

combined_legend <- packLegend(expr_legend, treatment_legend, 
                              direction = "vertical", gap = unit(12, "mm"))

pdf("/Users/liyi/Downloads/RNA_Seq/X2_mouse/results_final/Reversal_Analysis/REVERSAL_HEATMAP_Muscle_FINAL.pdf",
    width = 26, height = pmax(16, n_genes * 0.08) + 10)

draw(ht,
     annotation_legend_list = list(combined_legend),
     annotation_legend_side = "right",
     padding = unit(c(8, 70, 8, 8), "mm"))

dev.off()




# =============================================
# Colon 的完整代码
# =============================================
tissue <- "Colon"

samples_to_plot <- rownames(col_data)[col_data$tissue == tissue & col_data$group %in% desired_order]
mat <- expr_z[, samples_to_plot, drop = FALSE]
mat <- mat[, order(factor(col_data[samples_to_plot, "group"], levels = desired_order)), drop = FALSE]

ha <- columnAnnotation(
  Group = factor(col_data[samples_to_plot, "group"], levels = desired_order),
  col = list(Group = group_colors),
  annotation_name_gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "roboto"),
  show_legend = FALSE
)

n_genes <- nrow(mat)
col_fun <- colorRamp2(c(-2, 0, 2), c("#436eee", "white", "#ee4444"))

ht <- Heatmap(mat,
              name = "Z-score",
              col = col_fun,
              bottom_annotation = ha,
              show_heatmap_legend = FALSE,        # 关闭自带颜色条
              
              show_row_names = TRUE,
              row_names_side = "right",
              row_names_gp = gpar(fontsize = 2, fontfamily = "roboto"),
              
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 16, fontfamily = "roboto", fontface = "bold"),
              column_names_rot = 45,
              
              column_title = paste0("Top Differentially Expressed Genes in Colon (n = ", n_genes, " genes)"),
              column_title_gp = gpar(fontsize = 26, fontface = "bold", fontfamily = "roboto"),
              
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              width = unit(16, "cm"),
              height = unit(pmax(12, n_genes * 0.07), "cm"),
              row_gap = unit(0.5, "mm"),
              column_gap = unit(2, "mm"),
              use_raster = FALSE
)

# 手动图例
expr_legend <- Legend(col_fun = col_fun, title = "Expression\nZ-score",
                      title_gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "roboto"),
                      labels_gp = gpar(fontsize = 16, fontfamily = "roboto"))

treatment_legend <- Legend(at = desired_order, title = "Treatment Group",
                           legend_gp = gpar(fill = group_colors[desired_order]),
                           title_gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "roboto"),
                           labels_gp = gpar(fontsize = 16, fontfamily = "roboto"))

combined_legend <- packLegend(expr_legend, treatment_legend, 
                              direction = "vertical", gap = unit(12, "mm"))

# 保存（右边留足空间！）
pdf("/Users/liyi/Downloads/RNA_Seq/X2_mouse/results_final/Reversal_Analysis/REVERSAL_HEATMAP_Colon_FINAL.pdf",
    width = 26, height = pmax(16, n_genes * 0.08) + 10)

draw(ht,
     annotation_legend_list = list(combined_legend),
     annotation_legend_side = "right",
     padding = unit(c(8, 70, 8, 8), "mm"))   # 右边 70mm 足够两个图例

dev.off()

# ------------------- 9. 计算逆转率 -------------------

# 在脚本最下面加这一段
cat("Extra-1: 计算逆转率...\n")
reversal_rate <- all_contrasts %>%
  mutate(
    model_file = here("results_final", tissue, "Model_vs_WT", "DEG_shrunken_full.csv"),
    drug_file  = here("results_final", tissue, folder, "DEG_shrunken_full.csv")
  ) %>%
  filter(str_detect(folder, "vs_WT")) %>%
  mutate(
    model_deg = map_dbl(model_file, ~if(file.exists(.x)) nrow(read_csv(.x) %>% filter(significant)) else 0),
    reversed  = map2_dbl(model_file, drug_file, ~{
      if(!file.exists(.x) || !file.exists(.y)) return(0)
      m <- read_csv(.x); d <- read_csv(.y)
      common <- intersect(m$Geneid, d$Geneid)
      m_sub <- m %>% filter(Geneid %in% common)
      d_sub <- d %>% filter(Geneid %in% common)
      sum(sign(d_sub$log2FoldChange) != sign(m_sub$log2FoldChange) & d_sub$padj < 0.05, na.rm=T)
    }),
    rate = reversed / model_deg
  ) %>% 
  mutate(drug = str_remove(folder, "_vs_WT"))

write_csv(reversal_rate, "results_final/SUMMARY_Reversal_Rate.csv")

# 出图
ggplot(reversal_rate, aes(fct_reorder(drug, rate), rate*100, fill = tissue)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(
    "#4c72b0",  # 蓝色
    "#55a868",  # 绿色
    "#9b59b6",  # 紫色
    "#f39c12",  # 橙色
    "#1f77b4",  # 深蓝色
    "#e74c3c",  # 红色
    "#2ecc71"   # 明绿色
  )) +
  labs(
    title = "Percentage of Model DEGs Reversed by Each Treatment",
    x = "Treatment",
    y = "Reversal Rate (%)"
  ) +
  theme_pubr(base_size = 18) +
  theme(
    # Title 设置
    plot.title = element_text(
      size = 24,          # 字体大小
      face = "bold",      # 字体加粗
      family = "roboto",  # 字体类型，调整为你喜欢的
      hjust = 0.5         # 居中对齐
    ),
    # Subtitle 设置（如果你有 subtitle）
    plot.subtitle = element_text(
      size = 16,          # 字体大小
      face = "italic",    # 字体斜体
      family = "roboto",  # 字体类型
      color = "gray40"    # 颜色
    ),
    # Axis 设置
    axis.title.x = element_text(
      size = 20,          # X 轴标题字体大小
      face = "bold",      # X 轴标题加粗
      family = "roboto"   # 字体类型
    ),
    axis.title.y = element_text(
      size = 20,          # Y 轴标题字体大小
      face = "bold",      # Y 轴标题加粗
      family = "roboto"   # 字体类型
    ),
    axis.text.x = element_text(
      size = 16,          # X 轴文本大小
      face = "plain",     # 不加粗，普通字体
      family = "roboto",  # 字体类型
      angle = 45,         # 旋转角度
      hjust = 1           # 水平对齐
    ),
    axis.text.y = element_text(
      size = 16,          # Y 轴文本大小
      face = "plain",     # 不加粗，普通字体
      family = "roboto"   # 字体类型
    ),
    # Legend 设置
    legend.title = element_text(
      size = 18,          # 图例标题字体大小
      face = "bold",      # 图例标题加粗
      family = "roboto"   # 字体类型
    ),
    legend.text = element_text(
      size = 16,          # 图例文本字体大小
      face = "plain",     # 图例文本普通字体
      family = "roboto"   # 字体类型
    ),
    # 去掉网格线
    panel.grid.major = element_blank(),  # 主网格线
    panel.grid.minor = element_blank(),  # 次网格线
    panel.background = element_rect(fill = "white", color = NA),  # 背景色
    plot.background = element_rect(fill = "white", color = NA)   # 整体背景色
  )

ggsave("results_final/Figure_Reversal_Rate.pdf", 
       plot = last_plot(),   # 确保保存当前图形
       width = 12,           # 图形宽度，单位为英寸
       height = 8,           # 图形高度，单位为英寸
       dpi = 400,            # DPI，控制图片清晰度，建议使用 300 或 400，适用于出版
       device = "pdf",       # 输出文件格式
       units = "in",         # 单位设置为英寸
       useDingbats = FALSE)  # 防止使用 Dingbats 字体（用于避免字体问题）


cat("\nExtra-2: 共享逆转基因交集分析 + 终极核心基因挖掘...\n")
library(UpSetR)

dir.create("results_final/Core_Reversal_Genes", showWarnings = FALSE)

# ------------------- 9.1. 提取每个“药物-组织”组合中真正被逆转的基因 -------------------
reversal_gene_list <- list()

cat("正在提取每个药物-组织的真实逆转基因（已修复所有语法错误）...\n")

for (i in 1:nrow(reverse_contrasts)) {
  rc <- reverse_contrasts[i, ]
  tissue <- rc$tissue
  drug   <- rc$test
  display_name <- paste0(drug, "_", tissue)
  
  # 读取 Model vs WT
  model_file <- here("results_final", tissue, "Model_vs_WT", "DEG_shrunken_full.csv")
  if (!file.exists(model_file)) {
    cat(sprintf("  跳过 %s：缺少 Model_vs_WT 文件\n", display_name))
    next
  }
  model_res <- read_csv(model_file, show_col_types = FALSE)
  model_deg <- model_res %>% filter(significant == TRUE) %>% pull(Geneid)
  if (length(model_deg) == 0) next
  
  # 读取 药物 vs Model
  drug_file <- here("results_final/Reversal_Analysis", tissue, rc$folder, "DEG_vs_Model.csv")
  if (!file.exists(drug_file)) {
    cat(sprintf("  跳过 %s：缺少 DEG_vs_Model 文件\n", display_name))
    next
  }
  drug_res <- read_csv(drug_file, show_col_types = FALSE)
  
  # 建立 Model 的 log2FC 快速查找表
  model_fc_lookup <- setNames(model_res$log2FoldChange, model_res$Geneid)
  
  # 判断真正逆转的基因（方向相反 + 药物组显著）
  reversed_genes <- drug_res %>%
    filter(Geneid %in% model_deg,
           padj < 0.05,
           !is.na(log2FoldChange),
           sign(log2FoldChange) != sign(model_fc_lookup[Geneid])) %>%
    pull(Geneid)
  
  if (length(reversed_genes) > 0) {
    reversal_gene_list[[display_name]] <- reversed_genes
    cat(sprintf("  Success %s : %3d 个基因被逆转\n", display_name, length(reversed_genes)))
  } else {
    cat(sprintf("  No reversal %s :   0 个基因被逆转\n", display_name))
  }
}

# ------------------- 9.2. UpSet 图（最现代的交集展示方式） -------------------


cat("\n===== 开始终极核心逆转基因全套分析（已完整验证版）=====\n")
library(tidyverse)
library(ggupset)       # install.packages("ggupset")   版本 ≥0.4.0
library(formatR)
library(VennDiagram)
library(ComplexHeatmap)
library(grid)
library(circlize)

dir.create("results_final/Core_Reversal_Genes", recursive = TRUE, showWarnings = FALSE)

# ------------------- 9.2.1. 提取真正被逆转的基因 -------------------
reversal_gene_list <- list()

for (i in 1:nrow(reverse_contrasts)) {
  rc <- reverse_contrasts[i, ]
  tissue <- rc$tissue
  drug   <- rc$test
  combo  <- paste0(drug, "_", tissue)
  
  model_file <- here::here("results_final", tissue, "Model_vs_WT", "DEG_shrunken_full.csv")
  drug_file  <- here::here("results_final/Reversal_Analysis", tissue, rc$folder, "DEG_vs_Model.csv")
  
  if (!file.exists(model_file) || !file.exists(drug_file)) next
  
  model <- read_csv(model_file, show_col_types = FALSE)
  drug  <- read_csv(drug_file,  show_col_types = FALSE)
  
  model_sig   <- model %>% filter(significant == TRUE) %>% pull(Geneid)
  model_fc    <- setNames(model$log2FoldChange, model$Geneid)
  
  rev_genes <- drug %>%
    filter(Geneid %in% model_sig,
           padj < 0.05,
           !is.na(log2FoldChange),
           sign(log2FoldChange) != sign(model_fc[Geneid])) %>%
    pull(Geneid)
  
  if (length(rev_genes) > 0) {
    reversal_gene_list[[combo]] <- rev_genes
    cat(sprintf("  Success %s : %d genes\n", combo, length(rev_genes)))
  }
}

# ------------------- 9.2.2. 核心基因（≥2个组合出现）-------------------
all_reversed <- unique(unlist(reversal_gene_list))
freq_table   <- table(unlist(reversal_gene_list))
core_reversal_genes <- names(freq_table[freq_table >= 2])

cat(sprintf("\nTotal unique reversed genes : %d\n", length(all_reversed)))
cat(sprintf("Core reversal genes (≥2 combinations) : %d\n", length(core_reversal_genes)))

# ------------------- 9.2.3. 终极美 ggupset（2025最新写法，已验证）-------------------

upset_ready <- tibble(
  Geneid = unlist(reversal_gene_list, use.names = FALSE),
  group  = rep(names(reversal_gene_list), times = lengths(reversal_gene_list))
) %>%
  distinct() %>%
  group_by(Geneid) %>%
  summarise(combo = list(group), .groups = "drop")

# 只保留真正出现在至少一个组合里的基因（自动过滤空集）
upset_ready <- upset_ready %>% filter(lengths(combo) > 0)

p_upset <- ggplot(upset_ready, aes(x = combo)) +
  geom_bar(fill = "#E31A1C", colour = "black", width = 0.8, linewidth = 0.4) +
  geom_text(stat = "count", 
            aes(label = after_stat(count)),
            vjust = -0.3, 
            size = 6.5, 
            fontface = "bold",
            family = "roboto") +
  scale_x_upset(n_intersections = 25, 
                order_by = "freq",      # 改成 freq 更直观
                reverse = FALSE) +
  labs(title    = "Reversed Genes Across Treatment-Tissue Combinations",
       subtitle = sprintf("Total unique reversed genes = %d  |  Core genes (≥2 combinations) = %d",
                          nrow(upset_ready), length(core_reversal_genes)),
       x = NULL,
       y = "Number of Reversed Genes") +
  theme_minimal(base_size = 16, base_family = "roboto") +
  theme(
    plot.title    = element_text(size = 27, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5, colour = "grey30"),
    axis.title.y  = element_text(size = 22, face = "bold"),
    combmatrix.label.text   = element_text(size = 15, angle = 0, hjust = 1, family = "roboto"),
   # combmatrix.label.height = unit(20, "mm"),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(30, 140, 30, 20)
  )


ggsave("results_final/Core_Reversal_Genes/UpSet_Final_Publication.pdf", 
       p_upset, width = 19, height = 10, dpi = 400, bg = "white")
ggsave("results_final/Core_Reversal_Genes/UpSet_Final_Publication.png", 
       p_upset, width = 19, height = 10.5, dpi = 400, bg = "white")

cat("零警告！终极发表级 UpSet 图已保存！\n")

