#!/usr/bin/env Rscript
# =============================================================================
# RNA-seq + ATAC-seq Analysis Pipeline
# Project : PLK1 inhibitors (BI2536, BI6727) in CAOV3 and OVCAR3
# Server  : rbgo-server2
# Data    : /userhome/natalia/Documents/PLK1_Kadie_Agne/
#
# Usage (RStudio):
#   source("analysis_pipeline.R")
#   run_pipeline()
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(DESeq2)
  library(msigdbr)
  library(fgsea)
  library(UpSetR)
})

# =============================================================================
# PATHS - all hardcoded to your exact directory structure on rbgo-server2
# =============================================================================

BASE_DIR <- "/userhome/natalia/Documents/PLK1_Kadie_Agne"
ATAC_DIR <- file.path(BASE_DIR, "ATAC")
RNA_DIR  <- file.path(BASE_DIR, "RNAseq")

# ── ATAC nf-core output (merged-replicate = primary analysis level) ──────────
ATAC_MREP_BASE <- file.path(ATAC_DIR, "bwa/merged_replicate")
ATAC_MLIB_BASE <- file.path(ATAC_DIR, "bwa/merged_library")

# featureCounts peak x sample matrix (18 samples, consensus peaks)
ATAC_COUNTS_FILE <- file.path(
  ATAC_MREP_BASE,
  "macs2/broad_peak/consensus/consensus_peaks.mRp.clN.featureCounts.txt"
)

# DESeq2 RDS produced by nf-core (DESeqDataSet, all 18 samples)
ATAC_DESEQ2_RDS <- file.path(
  ATAC_MREP_BASE,
  "macs2/broad_peak/consensus/deseq2/consensus_peaks.mRp.clN.rds"
)

# HOMER annotatePeaks on consensus peaks
ATAC_CONSENSUS_ANNOT <- file.path(
  ATAC_MREP_BASE,
  "macs2/broad_peak/consensus/consensus_peaks.mRp.clN.annotatePeaks.txt"
)

# Per-condition HOMER annotated peaks (merged-replicate)
ATAC_MREP_PEAK_DIR <- file.path(ATAC_MREP_BASE, "macs2/broad_peak")
ATAC_MACS2_ANNOT <- list(
  CAOV3_DMSO   = file.path(ATAC_MREP_PEAK_DIR, "CAOV3_DMSO.mRp.clN_peaks.annotatePeaks.txt"),
  CAOV3_BI2536 = file.path(ATAC_MREP_PEAK_DIR, "CAOV3_BI2536.mRp.clN_peaks.annotatePeaks.txt"),
  CAOV3_BI6727 = file.path(ATAC_MREP_PEAK_DIR, "CAOV3_BI6727.mRp.clN_peaks.annotatePeaks.txt"),
  OVCAR3_DMSO  = file.path(ATAC_MREP_PEAK_DIR, "OVCAR3_DMSO.mRp.clN_peaks.annotatePeaks.txt"),
  OVCAR3_BI2536= file.path(ATAC_MREP_PEAK_DIR, "OVCAR3_BI2536.mRp.clN_peaks.annotatePeaks.txt"),
  OVCAR3_BI6727= file.path(ATAC_MREP_PEAK_DIR, "OVCAR3_BI6727.mRp.clN_peaks.annotatePeaks.txt")
)

# Picard insert size metrics (per-replicate, merged-library)
ATAC_PICARD_DIR <- file.path(ATAC_MLIB_BASE, "picard_metrics")

# MultiQC aggregated insert size table
ATAC_INSERT_SIZE_MQC <- file.path(
  ATAC_DIR,
  "multiqc/broad_peak/multiqc_data/mqc_picard_insert_size_Counts.txt"
)

# FRiP QC directory
ATAC_FRIP_DIR <- file.path(ATAC_MREP_BASE, "macs2/broad_peak/qc")

# ── RNA-seq ──────────────────────────────────────────────────────────────────
RNA_BAM_DIR <- file.path(RNA_DIR, "processed_files")

# Generate this with featureCounts (see generate_rna_featurecounts_cmd())
RNA_COUNT_MATRIX <- file.path(RNA_DIR, "count_matrix_featureCounts.txt")

# Pre-computed DESeq2 CSVs (optional; pipeline will run DESeq2 if absent)
RNA_DESEQ2_DIR <- file.path(RNA_DIR, "deseq2")

# =============================================================================
# SAMPLE METADATA
# =============================================================================

# RNA: full BAM stem (from processed_files/) mapped to biological identity
# NOTE: The mapping below is based on the sample sheet you provided.
# Stems marked PLACEHOLDER need to be confirmed with Natalia.
###THIS IS A WRONG ANNOTATION NEEDS FIXING ACCORDING TO THE EXCEL https://swanseauniversity-my.sharepoint.com/:x:/r/personal/l_francis_swansea_ac_uk/_layouts/15/Doc.aspx?sourcedoc=%7B40F65BAE-8115-415F-A09F-24D37E766C66%7D&file=RNA-seq_and_ATAC_seq_directories_for_KDMi_and_PLK1.xlsx&action=default&mobileredirect=true&DefaultItemOpen=1&wdOrigin=APPHOME-WEB.DIRECT%2CAPPHOME-WEB.UNAUTH%2CAPPHOME-WEB.SHELL.SIGNIN%2CAPPHOME-WEB.FILEBROWSER.RECENT&wdPreviousSession=41755b44-2218-4cf2-836b-34a4e52fc40d&wdPreviousSessionSrc=AppHomeWeb&ct=1777028613271

RNA_SAMPLE_MAP <- tribble(
  ~file_stem,                                          ~sample_id,           ~cell_line, ~condition, ~replicate,
  # CAOV3 BI6727 (RNA_CA_T*)
  "RNA_CA_T1_EKRN230042750-1A_HKMHTDSX7_L2",  "CAOV3_BI6727_REP1", "CAOV3",  "BI6727", "REP1",
  "RNA_CA_T2_EKRN230042751-1A_HKMCFDSX7_L1",  "CAOV3_BI6727_REP2", "CAOV3",  "BI6727", "REP2",
  "RNA_CA_T3_EKRN230042752-1A_HKMHTDSX7_L2",  "CAOV3_BI6727_REP3", "CAOV3",  "BI6727", "REP3",
  # OVCAR3 BI6727 (RNA_OV_T*)
  "RNA_OV_T1_EKRN230042754-1A_HKMHTDSX7_L2",  "OVCAR3_BI6727_REP1","OVCAR3", "BI6727", "REP1",
  "RNA_OV_T2_EKRN230042755-1A_HKMHTDSX7_L2",  "OVCAR3_BI6727_REP2","OVCAR3", "BI6727", "REP2",
  "RNA_OV_T3_EKRN230042756-1A_HKMGHDSX7_L1",  "OVCAR3_BI6727_REP3","OVCAR3", "BI6727", "REP3",
  # CAOV3 BI2536 (CA_J*)
  "CA_J5_3_EKRN230042760-1A_HKMHTDSX7_L2",    "CAOV3_BI2536_REP1", "CAOV3",  "BI2536", "REP1",
  "CA_J4_3_EKRN230042761-1A_HKMHTDSX7_L2",    "CAOV3_BI2536_REP2", "CAOV3",  "BI2536", "REP2",
  # PLACEHOLDER - confirm 3rd CAOV3 BI2536 BAM with Natalia:
  "AS_206_2D_EKRN230042769-1A_HKMHTDSX7_L2",  "CAOV3_BI2536_REP3", "CAOV3",  "BI2536", "REP3",
  # CAOV3 DMSO (SK_J* based on sample sheet CAD2/3/4 ordering)
  "SK_J5_3_EKRN230042758-1A_HKMHTDSX7_L2",    "CAOV3_DMSO_REP1",   "CAOV3",  "DMSO",   "REP1",
  "SK_J4_3_EKRN230042759-1A_HKMHTDSX7_L2",    "CAOV3_DMSO_REP2",   "CAOV3",  "DMSO",   "REP2",
  "OV_J5_3_EKRN230042762-1A_HKMHTDSX7_L2",    "CAOV3_DMSO_REP3",   "CAOV3",  "DMSO",   "REP3",
  # OVCAR3 BI2536 (OV_J4 + AS_186 + AS_200_2D)
  "OV_J4_3_EKRN230042763-1A_HKMHTDSX7_L2",    "OVCAR3_BI2536_REP1","OVCAR3", "BI2536", "REP1",
  "AS_186_2D_EKRN230042764-1A_HKMHTDSX7_L2",  "OVCAR3_BI2536_REP2","OVCAR3", "BI2536", "REP2",
  "AS_200_2D_EKRN230042765-1A_HKMHTDSX7_L2",  "OVCAR3_BI2536_REP3","OVCAR3", "BI2536", "REP3",
  # OVCAR3 DMSO (AS_200_3D, AS_197_2D, AS_197_3D)
  "AS_200_3D_EKRN230042766-1A_HKMHTDSX7_L2",  "OVCAR3_DMSO_REP1",  "OVCAR3", "DMSO",   "REP1",
  "AS_197_2D_EKRN230042767-1A_HKMHTDSX7_L2",  "OVCAR3_DMSO_REP2",  "OVCAR3", "DMSO",   "REP2",
  "AS_197_3D_EKRN230042768-1A_HKMHTDSX7_L2",  "OVCAR3_DMSO_REP3",  "OVCAR3", "DMSO",   "REP3"
) |>
  mutate(
    bam_path  = file.path(RNA_BAM_DIR, paste0(file_stem, ".bam")),
    condition = factor(condition, levels = c("DMSO","BI2536","BI6727")),
    cell_line = factor(cell_line,  levels = c("CAOV3","OVCAR3")),
    replicate = factor(replicate,  levels = c("REP1","REP2","REP3"))
  )

# ATAC: nf-core uses clean CL_COND_REP naming
ATAC_SAMPLE_MAP <- expand_grid(
  cell_line = c("CAOV3","OVCAR3"),
  condition = c("DMSO","BI2536","BI6727"),
  replicate = c("REP1","REP2","REP3")
) |>
  mutate(
    sample_id = paste(cell_line, condition, replicate, sep = "_"),
    condition = factor(condition, levels = c("DMSO","BI2536","BI6727")),
    cell_line = factor(cell_line,  levels = c("CAOV3","OVCAR3")),
    replicate = factor(replicate,  levels = c("REP1","REP2","REP3")),
    bam_mlib  = file.path(ATAC_MLIB_BASE,
                           paste0(sample_id, ".mLb.clN.sorted.bam")),
    insert_size_metrics = file.path(
      ATAC_PICARD_DIR,
      paste0(sample_id, ".mLb.clN.CollectMultipleMetrics.insert_size_metrics")
    )
  )

# =============================================================================
# CONSTANTS
# =============================================================================

CONDITIONS <- c("DMSO","BI2536","BI6727")
CELL_LINES <- c("CAOV3","OVCAR3")
REPLICATES <- c("REP1","REP2","REP3")

PALETTE <- c(
  DMSO="#4C72B0", BI2536="#DD8452", BI6727="#55A868",
  CAOV3="#C44E52", OVCAR3="#8172B2",
  REP1="#e8c63e", REP2="#6fbfcc", REP3="#b36fb8"
)

COMPARISONS <- list(
  list(cl="CAOV3",  ctrl="DMSO", trt="BI2536"),
  list(cl="CAOV3",  ctrl="DMSO", trt="BI6727"),
  list(cl="OVCAR3", ctrl="DMSO", trt="BI2536"),
  list(cl="OVCAR3", ctrl="DMSO", trt="BI6727")
)

GENE_SETS <- list(
  "Cell-cycle regulators" = c(
    "PLK1","CDC25C","AURKA","AURKB","CCNB1","CCNB2",
    "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","CDK1","CDK2","CCNA2","BUB1","BUB1B"
  ),
  "DNA damage / checkpoint" = c(
    "ATM","ATR","CHEK1","CHEK2","GADD45A","GADD45B","GADD45G",
    "RRM2","BRCA1","BRCA2","RAD51","FANCD2","H2AFX"
  ),
  "Apoptosis / autophagy" = c(
    "BCL2","BCL2L1","MCL1","BAX","BAD","BID","PUMA",
    "CASP3","CASP7","CASP9","CYCS","BECN1","ATG5","ATG7","BIRC5"
  )
)

# =============================================================================
# HELPERS
# =============================================================================

save_plot <- function(p, path, width=8, height=6) {
  dir.create(dirname(path), showWarnings=FALSE, recursive=TRUE)
  ggsave(path, plot=p, width=width, height=height, dpi=150, device="pdf")
  message("Saved -> ", path)
}

save_pheatmap <- function(ph, path, width=10, height=8) {
  dir.create(dirname(path), showWarnings=FALSE, recursive=TRUE)
  pdf(path, width=width, height=height); print(ph); dev.off()
  message("Saved -> ", path)
}

log2cpm <- function(counts, pseudo=1) {
  lib <- colSums(counts)
  log2(sweep(counts, 2, lib/1e6, "/") + pseudo)
}

# =============================================================================
# HELPER: featureCounts command to run on the server
# =============================================================================

generate_rna_featurecounts_cmd <- function(
    gtf_path,
    out_file = RNA_COUNT_MATRIX,
    threads  = 8,
    strandedness = 2   # 0=unstranded, 1=forward, 2=reverse
) {
  bams    <- RNA_SAMPLE_MAP$bam_path
  bam_str <- paste(bams, collapse=" \\\n  ")
  cmd <- paste0(
    "featureCounts \\\n",
    "  -T ", threads, " \\\n",
    "  -p --countReadPairs \\\n",
    "  -s ", strandedness, " \\\n",
    "  -a ", gtf_path, " \\\n",
    "  -o ", out_file, " \\\n",
    "  ", bam_str
  )
  message("Run this on rbgo-server2:\n\n", cmd, "\n")
  invisible(cmd)
}

# =============================================================================
# DATA LOADERS
# =============================================================================

load_rna_counts <- function(path=RNA_COUNT_MATRIX) {
  if (!file.exists(path)) {
    message("! RNA count matrix not found: ", path)
    message("  Call generate_rna_featurecounts_cmd('/path/to/genes.gtf') for the shell command.")
    return(NULL)
  }
  message("Loading RNA counts: ", path)
  raw <- read.delim(path, comment.char="#", check.names=FALSE)
  meta_cols <- c("Geneid","Chr","Start","End","Strand","Length")
  mat <- raw |> select(-any_of(setdiff(meta_cols,"Geneid"))) |>
    column_to_rownames("Geneid")
  stem_to_id <- setNames(RNA_SAMPLE_MAP$sample_id, RNA_SAMPLE_MAP$file_stem)
  names(mat) <- sapply(names(mat), function(col) {
    stem <- sub("\\.bam$","", basename(col))
    if (stem %in% names(stem_to_id)) stem_to_id[[stem]] else col
  })
  keep <- intersect(names(mat), RNA_SAMPLE_MAP$sample_id)
  message("  ", nrow(mat), " genes x ", length(keep), " samples")
  mat[, keep, drop=FALSE]
}

load_atac_counts <- function(path=ATAC_COUNTS_FILE) {
  if (!file.exists(path)) stop("ATAC featureCounts not found:\n  ", path)
  message("Loading ATAC featureCounts: ", path)
  raw <- read.delim(path, comment.char="#", check.names=FALSE)
  mat <- raw |>
    select(-any_of(c("Chr","Start","End","Strand","Length"))) |>
    column_to_rownames("Geneid")
  names(mat) <- sapply(names(mat), function(col) {
    stem <- sub("\\.(mRp|mLb).*","", basename(col))
    if (stem %in% ATAC_SAMPLE_MAP$sample_id) stem else col
  })
  keep <- intersect(names(mat), ATAC_SAMPLE_MAP$sample_id)
  message("  ", nrow(mat), " peaks x ", length(keep), " samples")
  mat[, keep, drop=FALSE]
}

load_atac_dds <- function(path=ATAC_DESEQ2_RDS) {
  if (!file.exists(path)) { message("ATAC DDS not found: ", path); return(NULL) }
  message("Loading ATAC DESeq2 RDS: ", path)
  readRDS(path)
}

load_homer_annot <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  df <- read.delim(path, check.names=FALSE, quote="")
  names(df)[1] <- "peak"
  annot_col <- grep("Annotation",    names(df), value=TRUE, ignore.case=TRUE)[1]
  gene_col  <- grep("Gene Name|Nearest.*Gene", names(df), value=TRUE, ignore.case=TRUE)[1]
  dist_col  <- grep("Distance.*TSS", names(df), value=TRUE, ignore.case=TRUE)[1]
  out <- tibble(peak=df$peak)
  if (!is.na(annot_col)) out$annotation_raw <- df[[annot_col]]
  if (!is.na(gene_col))  out$nearest_gene   <- df[[gene_col]]
  if (!is.na(dist_col))  out$distance_tss   <- as.numeric(df[[dist_col]])
  out |> mutate(annotation = case_when(
    str_detect(annotation_raw, regex("promoter",  ignore_case=TRUE)) ~ "Promoter",
    str_detect(annotation_raw, regex("^exon",     ignore_case=TRUE)) ~ "Exon",
    str_detect(annotation_raw, regex("intron",    ignore_case=TRUE)) ~ "Intron",
    str_detect(annotation_raw, regex("3'|TTS",    ignore_case=TRUE)) ~ "3'UTR",
    str_detect(annotation_raw, regex("5'",        ignore_case=TRUE)) ~ "5'UTR",
    TRUE ~ "Intergenic"
  ))
}

load_picard_insert_size <- function(f) {
  if (!file.exists(f)) return(NULL)
  lines <- readLines(f)
  hist_start <- grep("^insert_size\t", lines)
  if (!length(hist_start)) return(NULL)
  read.delim(f, skip=hist_start-1, header=TRUE, nrows=1000, check.names=FALSE) |>
    as_tibble() |> select(insert_size=1, count=2) |>
    filter(!is.na(insert_size), count>0)
}

# =============================================================================
# DESeq2: RNA
# =============================================================================

run_rna_deseq2 <- function(counts, sample_map, out_dir) {
  csv_dir <- file.path(out_dir,"rnaseq","deseq2")
  dir.create(csv_dir, recursive=TRUE, showWarnings=FALSE)
  results_list <- list()
  for (comp in COMPARISONS) {
    key <- paste0(comp$cl,"_",comp$trt,"_vs_",comp$ctrl)
    csv <- file.path(csv_dir, paste0(key,"_deseq2.csv"))
    if (file.exists(csv)) {
      message("  Loading RNA DESeq2: ", basename(csv))
      df <- read.csv(csv); names(df) <- tolower(names(df))
      if (!"gene" %in% names(df)) df$gene <- rownames(df)
      results_list[[key]] <- df; next
    }
    message("  Running RNA DESeq2: ", key)
    samps <- sample_map |> filter(cell_line==comp$cl, condition %in% c(comp$ctrl,comp$trt))
    keep  <- intersect(samps$sample_id, colnames(counts))
    if (length(keep)<4) { warning("Too few samples for ",key); next }
    sub   <- round(counts[,keep]); sub <- sub[rowSums(sub)>10,]
    cd    <- samps |> filter(sample_id %in% keep) |>
      select(sample_id,condition) |> column_to_rownames("sample_id") |>
      mutate(condition=factor(condition, levels=c(comp$ctrl,comp$trt)))
    dds   <- DESeqDataSetFromMatrix(sub[,rownames(cd)], cd, ~condition)
    dds   <- DESeq(dds, quiet=TRUE)
    res   <- results(dds, contrast=c("condition",comp$trt,comp$ctrl), alpha=0.05) |>
      as.data.frame() |> rownames_to_column("gene")
    write.csv(res, csv, row.names=FALSE)
    results_list[[key]] <- res
    message("    Saved: ", basename(csv))
  }
  results_list
}

# =============================================================================
# DESeq2: ATAC (pairwise from nf-core DDS)
# =============================================================================

run_atac_deseq2 <- function(dds, out_dir) {
  csv_dir <- file.path(out_dir,"atacseq","deseq2")
  dir.create(csv_dir, recursive=TRUE, showWarnings=FALSE)
  results_list <- list()
  if (is.null(dds)) return(results_list)
  message("ATAC DDS colData: ", paste(names(colData(dds)), collapse=", "))

  for (comp in COMPARISONS) {
    key <- paste0(comp$cl,"_",comp$trt,"_vs_",comp$ctrl)
    csv <- file.path(csv_dir, paste0(key,"_diff_atac.csv"))
    if (file.exists(csv)) {
      message("  Loading ATAC DESeq2: ", basename(csv))
      results_list[[key]] <- read.csv(csv); next
    }
    message("  Running ATAC DESeq2: ", key)
    tryCatch({
      samp_ids <- ATAC_SAMPLE_MAP$sample_id[ATAC_SAMPLE_MAP$cell_line==comp$cl]
      samp_ids <- intersect(samp_ids, colnames(dds))
      dds_sub  <- dds[, samp_ids]
      cd <- as.data.frame(colData(dds_sub))
      if ("condition" %in% names(cd)) {
        colData(dds_sub)$condition <- factor(cd$condition,
                                              levels=c(comp$ctrl,comp$trt))
      } else {
        inferred <- str_extract(colnames(dds_sub),"DMSO|BI2536|BI6727")
        colData(dds_sub)$condition <- factor(inferred,
                                              levels=c(comp$ctrl,comp$trt))
      }
      dds_sub <- dds_sub[, !is.na(colData(dds_sub)$condition)]
      design(dds_sub) <- ~condition
      dds_sub <- DESeq(dds_sub, quiet=TRUE)
      res <- results(dds_sub, contrast=c("condition",comp$trt,comp$ctrl),
                     alpha=0.05) |> as.data.frame() |> rownames_to_column("peak")
      annot <- load_homer_annot(ATAC_CONSENSUS_ANNOT)
      if (!is.null(annot)) res <- left_join(res, annot, by="peak")
      write.csv(res, csv, row.names=FALSE)
      results_list[[key]] <- res
      message("    Saved: ", basename(csv))
    }, error=function(e) warning("ATAC DESeq2 failed ",key,": ",conditionMessage(e)))
  }
  results_list
}

# =============================================================================
# FIGURES
# =============================================================================

# ── PCA ───────────────────────────────────────────────────────────────────────
plot_pca <- function(log_mat, sample_map, out_dir, assay="rna") {
  message(toupper(assay),": PCA")
  subdir <- if(assay=="rna") "rnaseq" else "atacseq"
  shapes <- c(REP1=16,REP2=15,REP3=17)
  for (cl in CELL_LINES) {
    samps <- intersect(sample_map$sample_id[sample_map$cell_line==cl], colnames(log_mat))
    if (length(samps)<3) next
    rv  <- apply(log_mat[,samps],1,var)
    top <- names(sort(rv,decreasing=TRUE))[1:min(5000,nrow(log_mat))]
    pca <- prcomp(t(log_mat[top,samps]), scale.=TRUE)
    var <- summary(pca)$importance["Proportion of Variance",]
    df  <- as_tibble(pca$x[,1:2]) |>
      mutate(sample_id=rownames(pca$x),
             condition=sample_map$condition[match(sample_id,sample_map$sample_id)],
             replicate=sample_map$replicate[match(sample_id,sample_map$sample_id)])
    p <- ggplot(df,aes(PC1,PC2,colour=condition,shape=replicate,label=sample_id))+
      geom_point(size=4,stroke=0.8)+
      geom_text_repel(size=2.5,show.legend=FALSE)+
      scale_colour_manual(values=PALETTE[CONDITIONS])+
      scale_shape_manual(values=shapes)+
      labs(title=paste0(if(assay=="rna")"RNA-seq" else "ATAC-seq"," PCA - ",cl),
           x=paste0("PC1 (",percent(var[1],0.1),")"),
           y=paste0("PC2 (",percent(var[2],0.1),")"),
           colour="Condition",shape="Replicate")+
      theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
    save_plot(p, file.path(out_dir,subdir,paste0("pca_",cl,".pdf")))
  }
}

# ── DEG heatmap ────────────────────────────────────────────────────────────────
plot_deg_heatmap <- function(log_counts, sample_map, deseq_results, out_dir, n_top=50) {
  message("RNA-seq: DEG heatmap")
  for (cl in CELL_LINES) {
    sig_genes <- character(0)
    for (comp in COMPARISONS) {
      if (comp$cl!=cl) next
      key <- paste0(comp$cl,"_",comp$trt,"_vs_",comp$ctrl)
      df  <- deseq_results[[key]]; if (is.null(df)) next
      names(df) <- tolower(names(df))
      gene_col  <- intersect(c("gene","geneid"),names(df))[1]
      top <- df |> filter(!is.na(padj),!is.na(log2foldchange),
                          padj<0.05,abs(log2foldchange)>1) |>
        slice_min(padj,n=n_top) |> pull(.data[[gene_col]])
      sig_genes <- union(sig_genes,top)
    }
    if (length(sig_genes)<2)
      sig_genes <- names(sort(apply(log_counts,1,var),decreasing=TRUE))[1:n_top]
    samps <- intersect(sample_map$sample_id[sample_map$cell_line==cl],colnames(log_counts))
    mat   <- log_counts[intersect(sig_genes,rownames(log_counts)),samps]
    z_mat <- t(scale(t(mat))); z_mat[is.nan(z_mat)] <- 0
    ann_col <- sample_map |> filter(sample_id %in% samps) |>
      select(sample_id,Condition=condition,Replicate=replicate) |>
      column_to_rownames("sample_id")
    ph <- pheatmap(z_mat,
      annotation_col=ann_col[colnames(z_mat),,drop=FALSE],
      annotation_colors=list(Condition=PALETTE[CONDITIONS],Replicate=PALETTE[REPLICATES]),
      color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
      breaks=seq(-3,3,length.out=101),cluster_rows=TRUE,cluster_cols=TRUE,
      show_rownames=nrow(z_mat)<=60,fontsize_row=7,fontsize_col=8,
      main=paste("Top DEGs Heatmap -",cl),silent=TRUE)
    save_pheatmap(ph,file.path(out_dir,"rnaseq",paste0("heatmap_top_degs_",cl,".pdf")),
                  width=11,height=max(8,nrow(z_mat)*0.22))
  }
}

# ── Volcano ────────────────────────────────────────────────────────────────────
plot_volcano <- function(results_list, out_dir, lfc_cut=1, padj_cut=0.05,
                          n_label=15, subdir="rnaseq") {
  message(toupper(subdir),": Volcano plots")
  for (comp in COMPARISONS) {
    key <- paste0(comp$cl,"_",comp$trt,"_vs_",comp$ctrl)
    df  <- results_list[[key]]; if (is.null(df)) next
    names(df) <- tolower(names(df))
    lfc_col  <- intersect(c("log2foldchange","lfc"),names(df))[1]
    gene_col <- intersect(c("gene","geneid","peak"),names(df))[1]
    df <- df |> rename(lfc=!!lfc_col) |>
      filter(!is.na(padj),!is.na(lfc)) |>
      mutate(neg_log10p=pmin(-log10(padj),300),
             dir=case_when(padj<padj_cut&lfc> lfc_cut~"Up",
                           padj<padj_cut&lfc< -lfc_cut~"Down",TRUE~"NS"))
    top <- df |> filter(dir!="NS") |> slice_min(padj,n=n_label)
    p <- ggplot(df,aes(lfc,neg_log10p,colour=dir))+
      geom_point(data=filter(df,dir=="NS"),alpha=0.3,size=0.7)+
      geom_point(data=filter(df,dir!="NS"),alpha=0.7,size=1.2)+
      geom_text_repel(data=top,aes(label=.data[[gene_col]]),
                      size=2.5,max.overlaps=20,show.legend=FALSE)+
      geom_hline(yintercept=-log10(padj_cut),linetype="dashed",linewidth=0.5)+
      geom_vline(xintercept=c(-lfc_cut,lfc_cut),linetype="dashed",linewidth=0.5)+
      scale_colour_manual(values=c(Up="#d62728",Down="#1f77b4",NS="lightgrey"))+
      annotate("text",x=-Inf,y=Inf,
               label=paste0("Up:",sum(df$dir=="Up"),"  Down:",sum(df$dir=="Down")),
               hjust=-0.1,vjust=1.5,size=3.5)+
      labs(title=paste0("Volcano: ",comp$cl," ",comp$trt," vs ",comp$ctrl),
           x="log2 FC",y="-log10 adj.p",colour=NULL)+
      theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
    save_plot(p,file.path(out_dir,subdir,paste0("volcano_",key,".pdf")))
  }
}

# ── GSEA ────────────────────────────────────────────────────────────────────────
run_gsea_analysis <- function(deseq_results, out_dir) {
  message("RNA-seq: GSEA")
  csv_dir <- file.path(out_dir,"rnaseq","gsea_csv")
  dir.create(csv_dir,showWarnings=FALSE,recursive=TRUE)
  h_sets  <- msigdbr(species="Homo sapiens",category="H")
  c2_sets <- msigdbr(species="Homo sapiens",category="C2",subcategory="CP:REACTOME")
  gene_sets <- c(split(h_sets$gene_symbol,h_sets$gs_name),
                 split(c2_sets$gene_symbol,c2_sets$gs_name))
  for (comp in COMPARISONS) {
    key <- paste0(comp$cl,"_",comp$trt,"_vs_",comp$ctrl)
    df  <- deseq_results[[key]]; if (is.null(df)) next
    names(df) <- tolower(names(df))
    lfc_col  <- intersect(c("log2foldchange","lfc"),names(df))[1]
    gene_col <- intersect(c("gene","geneid"),names(df))[1]
    ranks <- df |> filter(!is.na(.data[[lfc_col]]),!is.na(pvalue),pvalue>0) |>
      mutate(m=sign(.data[[lfc_col]])*-log10(pvalue)) |>
      arrange(desc(m)) |> { setNames(.$m, toupper(.[[gene_col]])) }()
    set.seed(42)
    res <- tryCatch(fgsea(gene_sets,ranks,minSize=10,maxSize=500,nPermSimple=1000),
                    error=function(e){warning(conditionMessage(e));NULL})
    if (is.null(res)) next
    res_df <- as.data.frame(res)|>select(-leadingEdge)|>arrange(padj)
    write.csv(res_df,file.path(csv_dir,paste0(key,"_gsea.csv")),row.names=FALSE)
    # dotplot
    PRIORITY <- c("HALLMARK_G2M_CHECKPOINT","HALLMARK_MITOTIC_SPINDLE",
                  "HALLMARK_E2F_TARGETS","HALLMARK_DNA_REPAIR",
                  "HALLMARK_APOPTOSIS","HALLMARK_P53_PATHWAY",
                  "REACTOME_CELL_CYCLE","REACTOME_DNA_REPLICATION","REACTOME_APOPTOSIS")
    df2 <- res_df |> filter(padj<0.25) |>
      mutate(neg_log10fdr=pmin(-log10(padj),10),priority=pathway %in% PRIORITY) |>
      arrange(desc(priority),desc(abs(NES))) |> slice_head(n=25) |>
      mutate(pathway=str_replace_all(pathway,"_"," "),pathway=fct_reorder(pathway,NES))
    if (nrow(df2)>0) {
      p <- ggplot(df2,aes(NES,pathway))+geom_vline(xintercept=0,linewidth=0.5)+
        geom_point(aes(size=neg_log10fdr,colour=NES))+
        scale_colour_gradient2(low="#1f77b4",mid="white",high="#d62728",midpoint=0,name="NES")+
        scale_size_continuous(name="-log10 FDR",range=c(2,8))+
        labs(title=paste("GSEA -",key),x="NES",y=NULL)+
        theme_bw(base_size=11)+theme(plot.title=element_text(face="bold"),axis.text.y=element_text(size=8))
      save_plot(p,file.path(out_dir,"rnaseq",paste0("gsea_dotplot_",key,".pdf")),
                width=9,height=max(5,nrow(df2)*0.35))
    }
    message("  Saved GSEA: ",key)
  }
}

# ── Gene-set heatmaps ────────────────────────────────────────────────────────────
plot_gene_set_heatmaps <- function(log_counts, sample_map, out_dir) {
  message("RNA-seq: Gene-set heatmaps")
  ann_col <- sample_map |>
    select(sample_id,`Cell line`=cell_line,Condition=condition,Replicate=replicate) |>
    column_to_rownames("sample_id")
  ann_colours <- list(`Cell line`=setNames(PALETTE[CELL_LINES],CELL_LINES),
                      Condition=setNames(PALETTE[CONDITIONS],CONDITIONS),
                      Replicate=setNames(PALETTE[REPLICATES],REPLICATES))
  for (set_name in names(GENE_SETS)) {
    present <- intersect(GENE_SETS[[set_name]],rownames(log_counts))
    if (length(present)<2) next
    mat <- log_counts[present,]; z_mat <- t(scale(t(mat))); z_mat[is.nan(z_mat)] <- 0
    ph <- pheatmap(z_mat,annotation_col=ann_col[colnames(z_mat),,drop=FALSE],
                   annotation_colors=ann_colours,
                   color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
                   breaks=seq(-2,2,length.out=101),cluster_rows=TRUE,cluster_cols=TRUE,
                   show_rownames=TRUE,fontsize_row=8,
                   main=paste("Gene-set Heatmap:",set_name),silent=TRUE)
    save_pheatmap(ph,file.path(out_dir,"rnaseq",
                               paste0("geneset_heatmap_",gsub("[/ ]","_",set_name),".pdf")),
                  width=14,height=max(5,length(present)*0.38))
  }
}

# ── UpSet overlap ─────────────────────────────────────────────────────────────
plot_deg_overlap <- function(deseq_results, out_dir, lfc_cut=1, padj_cut=0.05) {
  message("RNA-seq: UpSet DEG overlap")
  for (cl in CELL_LINES) {
    deg_sets <- list()
    for (comp in COMPARISONS) {
      if (comp$cl!=cl) next
      key <- paste0(comp$cl,"_",comp$trt,"_vs_",comp$ctrl)
      df  <- deseq_results[[key]]; if (is.null(df)) next
      names(df) <- tolower(names(df))
      gene_col  <- intersect(c("gene","geneid"),names(df))[1]
      deg_sets[[comp$trt]] <- df |>
        filter(!is.na(padj),padj<padj_cut,abs(log2foldchange)>lfc_cut) |>
        pull(.data[[gene_col]])
    }
    if (length(deg_sets)<2) next
    all_g <- unique(unlist(deg_sets))
    mat   <- sapply(deg_sets,function(s) as.integer(all_g %in% s))
    rownames(mat) <- all_g
    path <- file.path(out_dir,"rnaseq",paste0("upset_deg_overlap_",cl,".pdf"))
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    pdf(path,width=8,height=5)
    upset(as.data.frame(mat),sets=names(deg_sets),
          sets.bar.color=unname(PALETTE[names(deg_sets)]),
          order.by="freq",main.bar.color="grey30",text.scale=1.3,mb.ratio=c(0.6,0.4))
    title(paste("DEG Overlap -",cl),line=3); dev.off()
    message("Saved -> ",path)
  }
}

# ── Insert size (real Picard/MultiQC data) ────────────────────────────────────
plot_insert_size <- function(out_dir) {
  message("ATAC-seq: Insert size")
  if (file.exists(ATAC_INSERT_SIZE_MQC)) {
    mqc  <- read.delim(ATAC_INSERT_SIZE_MQC,check.names=FALSE)
    long <- mqc |> rename(insert_size=1) |>
      pivot_longer(-insert_size,names_to="samp",values_to="count") |>
      mutate(samp=str_replace(samp," .*",""),
             cell_line=str_extract(samp,"CAOV3|OVCAR3"),
             condition=str_extract(samp,"DMSO|BI2536|BI6727"),
             replicate=str_extract(samp,"REP[1-3]")) |>
      filter(!is.na(cell_line),!is.na(condition),count>0) |>
      group_by(samp) |> mutate(freq=count/sum(count)) |> ungroup()
  } else {
    long <- ATAC_SAMPLE_MAP |>
      mutate(data=map(insert_size_metrics,load_picard_insert_size)) |>
      filter(!map_lgl(data,is.null)) |> unnest(data) |>
      group_by(sample_id) |> mutate(freq=count/sum(count)) |> ungroup()
  }
  plots <- lapply(CELL_LINES, function(cl) {
    sub <- filter(long,cell_line==cl)
    ggplot(sub,aes(insert_size,freq,colour=condition,
                   group=interaction(condition,replicate)))+
      geom_line(linewidth=0.7,alpha=0.8)+
      scale_colour_manual(values=PALETTE[CONDITIONS])+
      scale_x_continuous(limits=c(0,800))+
      labs(title=paste(cl,"Insert Size"),x="Insert size (bp)",
           y="Relative frequency",colour="Condition")+
      theme_bw(base_size=11)+theme(plot.title=element_text(face="bold"))
  })
  p <- wrap_plots(plots,ncol=1)+
    plot_annotation(title="ATAC-seq Insert Size Distributions (Supplementary)",
                    theme=theme(plot.title=element_text(face="bold",size=13)))
  save_plot(p,file.path(out_dir,"atacseq","supp_insert_size.pdf"),width=9,height=9)
}

# ── FRiP ────────────────────────────────────────────────────────────────────────
plot_frip <- function(out_dir) {
  message("ATAC-seq: FRiP scores")
  files <- list.files(ATAC_FRIP_DIR,pattern="FRiP_mqc\\.tsv$",full.names=TRUE)
  if (!length(files)) { message("  No FRiP files found"); return(invisible(NULL)) }
  df <- map_dfr(files,~read.delim(.x,header=FALSE,col.names=c("sample","frip"))) |>
    mutate(cell_line=str_extract(sample,"CAOV3|OVCAR3"),
           condition=str_extract(sample,"DMSO|BI2536|BI6727")) |>
    filter(!is.na(cell_line))
  p <- ggplot(df,aes(condition,frip,fill=condition))+
    geom_col(colour="black",linewidth=0.4,width=0.6,position=position_dodge(0.7))+
    geom_hline(yintercept=0.2,linetype="dashed",colour="red",linewidth=0.7)+
    facet_wrap(~cell_line)+
    scale_fill_manual(values=PALETTE[CONDITIONS],guide="none")+
    labs(title="ATAC-seq FRiP Scores (Supplementary)",x=NULL,y="FRiP")+
    theme_bw(base_size=11)+theme(plot.title=element_text(face="bold"))
  save_plot(p,file.path(out_dir,"atacseq","supp_frip_scores.pdf"),width=9,height=5)
}

# ── ATAC MA + Volcano ────────────────────────────────────────────────────────────
plot_atac_diff <- function(diff_atac, out_dir, lfc_cut=1, padj_cut=0.05) {
  message("ATAC-seq: MA + Volcano")
  key_genes <- unique(unlist(GENE_SETS))
  for (comp in COMPARISONS) {
    key <- paste0(comp$cl,"_",comp$trt,"_vs_",comp$ctrl)
    df  <- diff_atac[[key]]; if (is.null(df)) next
    names(df) <- tolower(names(df))
    lfc_col <- intersect(c("log2foldchange","lfc"),names(df))[1]
    bm_col  <- intersect(c("basemean","mean_acc"),names(df))[1]
    df <- df |> rename(lfc=!!lfc_col) |> filter(!is.na(padj),!is.na(lfc)) |>
      mutate(neg_log10p=pmin(-log10(padj),300),
             sig=padj<padj_cut&abs(lfc)>lfc_cut,
             key_gene=if("nearest_gene" %in% names(.)) nearest_gene %in% key_genes else FALSE,
             cat=case_when(sig&key_gene~"Key gene",sig&lfc>0~"Gained",
                           sig&lfc<0~"Lost",TRUE~"NS"))
    cols <- scale_colour_manual(
      values=c("Key gene"="#ff7f0e","Gained"="#d62728","Lost"="#1f77b4","NS"="lightgrey"),
      name=NULL)
    ma <- ggplot(df,aes(log2(.data[[bm_col]]+1),lfc,colour=cat))+
      geom_point(data=filter(df,cat=="NS"),alpha=0.2,size=0.5)+
      geom_point(data=filter(df,cat!="NS"),alpha=0.6,size=0.9)+
      geom_hline(yintercept=c(-lfc_cut,lfc_cut),linetype="dashed",linewidth=0.5)+
      cols+labs(title=paste("MA -",key),x="log2 mean acc+1",y="log2 FC")+
      theme_bw(base_size=11)+theme(plot.title=element_text(face="bold"),legend.position="bottom")
    vol <- ggplot(df,aes(lfc,neg_log10p,colour=cat))+
      geom_point(data=filter(df,cat=="NS"),alpha=0.2,size=0.5)+
      geom_point(data=filter(df,cat!="NS"),alpha=0.6,size=0.9)+
      geom_hline(yintercept=-log10(padj_cut),linetype="dashed",linewidth=0.5)+
      geom_vline(xintercept=c(-lfc_cut,lfc_cut),linetype="dashed",linewidth=0.5)+
      cols+labs(title=paste("Volcano -",key),x="log2 FC",y="-log10 adj.p")+
      theme_bw(base_size=11)+theme(plot.title=element_text(face="bold"),legend.position="bottom")
    combined <- ma+vol+plot_layout(guides="collect")&theme(legend.position="bottom")
    save_plot(combined,file.path(out_dir,"atacseq",paste0("diff_atac_",key,".pdf")),width=13,height=6)
  }
}

# ── Peak annotation bar ──────────────────────────────────────────────────────────
plot_peak_annotation <- function(diff_atac, out_dir, lfc_cut=1, padj_cut=0.05) {
  message("ATAC-seq: Peak annotation bar chart")
  regions <- c("Promoter","Exon","Intron","Intergenic","3'UTR","5'UTR")
  rcols   <- c(Promoter="#e41a1c",Exon="#377eb8",Intron="#4daf4a",
               Intergenic="#984ea3","3'UTR"="#ff7f00","5'UTR"="#a65628")
  frames <- lapply(COMPARISONS,function(comp){
    key <- paste0(comp$cl,"_",comp$trt,"_vs_",comp$ctrl)
    df  <- diff_atac[[key]]; if(is.null(df)) return(NULL)
    names(df) <- tolower(names(df))
    if(!"annotation" %in% names(df)) return(NULL)
    df |> filter(!is.na(padj),padj<padj_cut,abs(log2foldchange)>lfc_cut) |>
      mutate(direction=if_else(log2foldchange>0,"Gained","Lost"),
             comparison=paste0(comp$cl,"\n",comp$trt," vs ",comp$ctrl),
             annotation=factor(annotation,levels=regions))
  })
  all_df <- bind_rows(frames)
  if (nrow(all_df)==0) return(invisible(NULL))
  plot_df <- all_df |> count(comparison,annotation,direction) |>
    mutate(n=if_else(direction=="Lost",-n,n))
  p <- ggplot(plot_df,aes(annotation,n,fill=annotation,alpha=direction))+
    geom_col(colour="grey20",linewidth=0.3,width=0.7)+
    geom_hline(yintercept=0,linewidth=0.5)+
    facet_wrap(~comparison,ncol=2)+
    scale_fill_manual(values=rcols,guide="none")+
    scale_alpha_manual(values=c(Gained=0.9,Lost=0.5),name=NULL)+
    labs(title="Differential Peaks by Genomic Region",x=NULL,
         y="Peak count (+gained / -lost)")+
    theme_bw(base_size=11)+
    theme(plot.title=element_text(face="bold"),
          axis.text.x=element_text(angle=35,hjust=1),strip.text=element_text(size=8))
  save_plot(p,file.path(out_dir,"atacseq","peak_region_annotation.pdf"),width=12,height=10)
}

# ── RNA-ATAC correlation ─────────────────────────────────────────────────────────
plot_rna_atac_correlation <- function(rna_res, diff_atac, out_dir) {
  message("Integrated: RNA-ATAC correlation")
  pw_lookup <- imap(GENE_SETS,~tibble(gene=.x,pathway=.y)) |> bind_rows()
  pw_pal <- c("Cell-cycle regulators"="#e41a1c","DNA damage / checkpoint"="#377eb8",
              "Apoptosis / autophagy"="#4daf4a","Other"="lightgrey")
  for (comp in COMPARISONS) {
    key  <- paste0(comp$cl,"_",comp$trt,"_vs_",comp$ctrl)
    rna  <- rna_res[[key]]; atac <- diff_atac[[key]]
    if (is.null(rna)||is.null(atac)) next
    names(rna) <- tolower(names(rna)); names(atac) <- tolower(names(atac))
    if (!"nearest_gene" %in% names(atac)) next
    prom <- if("annotation" %in% names(atac))
              filter(atac,str_detect(tolower(annotation),"promoter")) else atac
    ag <- prom |> group_by(nearest_gene) |>
      summarise(atac_lfc=mean(log2foldchange,na.rm=TRUE)) |> rename(gene=nearest_gene)
    gc <- intersect(c("gene","geneid"),names(rna))[1]
    mg <- rna |> select(gene=!!gc,rna_lfc=log2foldchange) |>
      inner_join(ag,by="gene") |> left_join(pw_lookup,by="gene") |>
      replace_na(list(pathway="Other")) |> filter(!is.na(rna_lfc),!is.na(atac_lfc))
    if (nrow(mg)<10) next
    rho <- cor(mg$atac_lfc,mg$rna_lfc,method="spearman",use="complete.obs")
    p <- ggplot(mg,aes(atac_lfc,rna_lfc,colour=pathway))+
      geom_point(data=filter(mg,pathway=="Other"),alpha=0.15,size=0.8)+
      geom_point(data=filter(mg,pathway!="Other"),alpha=0.8,size=2.5)+
      geom_text_repel(data=filter(mg,pathway!="Other"),aes(label=gene),
                      size=2.3,max.overlaps=15,show.legend=FALSE)+
      geom_smooth(data=mg,aes(atac_lfc,rna_lfc),method="lm",colour="black",
                  linewidth=0.7,se=TRUE,inherit.aes=FALSE)+
      geom_hline(yintercept=0,linewidth=0.4)+geom_vline(xintercept=0,linewidth=0.4)+
      scale_colour_manual(values=pw_pal,name="Gene set")+
      labs(title=paste("RNA-ATAC Correlation -",key),
           subtitle=paste0("Spearman rho = ",round(rho,2)),
           x="Promoter accessibility log2FC (ATAC)",y="Expression log2FC (RNA)")+
      theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
    save_plot(p,file.path(out_dir,"integrated",paste0("rna_atac_correlation_",key,".pdf")))
  }
}

# ── Four-quadrant ─────────────────────────────────────────────────────────────────
plot_four_quadrant <- function(rna_res, diff_atac, out_dir, lfc_cut=1, padj_cut=0.05) {
  message("Integrated: Four-quadrant")
  q_pal <- c("+ATAC/+RNA (activated)"="#2ca02c","-ATAC/-RNA (repressed)"="#d62728",
             "-ATAC/+RNA (discordant)"="#9467bd","+ATAC/-RNA (discordant)"="#ff7f0e",
             "Not significant"="lightgrey")
  for (comp in COMPARISONS) {
    key <- paste0(comp$cl,"_",comp$trt,"_vs_",comp$ctrl)
    rna <- rna_res[[key]]; atac <- diff_atac[[key]]
    if (is.null(rna)||is.null(atac)||!"nearest_gene" %in% tolower(names(atac))) next
    names(rna) <- tolower(names(rna)); names(atac) <- tolower(names(atac))
    prom <- if("annotation" %in% names(atac))
              filter(atac,str_detect(tolower(annotation),"promoter")) else atac
    ag <- prom |> group_by(nearest_gene) |>
      summarise(atac_lfc=mean(log2foldchange,na.rm=TRUE),atac_padj=min(padj,na.rm=TRUE)) |>
      rename(gene=nearest_gene)
    gc <- intersect(c("gene","geneid"),names(rna))[1]
    mg <- rna |> select(gene=!!gc,rna_lfc=log2foldchange,rna_padj=padj) |>
      inner_join(ag,by="gene") |> filter(!is.na(rna_padj),!is.na(atac_padj)) |>
      mutate(rna_sig=rna_padj<padj_cut&abs(rna_lfc)>lfc_cut,
             atac_sig=atac_padj<padj_cut&abs(atac_lfc)>lfc_cut,
             quad=case_when(
               rna_sig&atac_sig&rna_lfc>0&atac_lfc>0~"+ATAC/+RNA (activated)",
               rna_sig&atac_sig&rna_lfc<0&atac_lfc<0~"-ATAC/-RNA (repressed)",
               rna_sig&atac_sig&rna_lfc>0&atac_lfc<0~"-ATAC/+RNA (discordant)",
               rna_sig&atac_sig&rna_lfc<0&atac_lfc>0~"+ATAC/-RNA (discordant)",
               TRUE~"Not significant"))
    cnts <- count(mg,quad)
    p <- ggplot(mg,aes(atac_lfc,rna_lfc,colour=quad))+
      geom_point(data=filter(mg,quad=="Not significant"),alpha=0.15,size=0.7)+
      geom_point(data=filter(mg,quad!="Not significant"),alpha=0.7,size=1.5)+
      geom_hline(yintercept=c(-lfc_cut,0,lfc_cut),
                 linetype=c("dashed","solid","dashed"),linewidth=0.5)+
      geom_vline(xintercept=c(-lfc_cut,0,lfc_cut),
                 linetype=c("dashed","solid","dashed"),linewidth=0.5)+
      scale_colour_manual(values=q_pal,
        labels=~paste0(.x," (",cnts$n[match(.x,cnts$quad)],")"),name=NULL)+
      labs(title=paste("Four-Quadrant RNA-ATAC -",key),
           x="Promoter accessibility log2FC (ATAC)",y="Expression log2FC (RNA)")+
      guides(colour=guide_legend(override.aes=list(size=3,alpha=1),ncol=1))+
      theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
    save_plot(p,file.path(out_dir,"integrated",paste0("four_quadrant_",key,".pdf")),width=9,height=7)
  }
}

# ── Drug-specific signatures ──────────────────────────────────────────────────────
plot_drug_specific <- function(rna_res, diff_atac, out_dir, lfc_cut=1, padj_cut=0.05) {
  message("Integrated: Drug-specific signatures")
  for (cl in CELL_LINES) {
    make_set <- function(drug, results) {
      key <- paste0(cl,"_",drug,"_vs_DMSO")
      df  <- results[[key]]; if(is.null(df)) return(character(0))
      names(df) <- tolower(names(df))
      gc <- intersect(c("gene","geneid","nearest_gene"),names(df))[1]
      df |> filter(!is.na(padj),padj<padj_cut,abs(log2foldchange)>lfc_cut) |>
        pull(.data[[gc]]) |> unique()
    }
    for (dtype in list(list(res=rna_res,label="DEGs",outpfx="drug_specific_degs"),
                       list(res=diff_atac,label="ATAC peaks",outpfx="drug_specific_atac"))) {
      s36 <- make_set("BI2536",dtype$res); s67 <- make_set("BI6727",dtype$res)
      if (!length(c(s36,s67))) next
      summary_df <- tibble(
        Group=c("BI2536 only","Shared","BI6727 only"),
        n=c(length(setdiff(s36,s67)),length(intersect(s36,s67)),length(setdiff(s67,s36))),
        fill=c(PALETTE["BI2536"],"#888888",PALETTE["BI6727"])
      ) |> mutate(Group=factor(Group,levels=Group))
      p <- ggplot(summary_df,aes(Group,n,fill=Group))+
        geom_col(colour="black",linewidth=0.5,width=0.6)+
        geom_text(aes(label=n),vjust=-0.4,size=4)+
        scale_fill_manual(values=setNames(summary_df$fill,summary_df$Group),guide="none")+
        labs(title=paste(dtype$label,"- Drug-Specific -",cl),x=NULL,y=paste("Number of",dtype$label))+
        ylim(0,max(summary_df$n)*1.12)+
        theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"))
      save_plot(p,file.path(out_dir,"integrated",paste0(dtype$outpfx,"_",cl,".pdf")),width=6,height=5)
    }
  }
}

# =============================================================================
# MAIN
# =============================================================================

run_pipeline <- function(out_dir   = file.path(BASE_DIR, "results"),
                          skip_gsea = FALSE) {
  message(strrep("=",60))
  message("PLK1 inhibitor RNA-seq + ATAC-seq Pipeline")
  message("Output: ", out_dir)
  message(strrep("=",60))
  for (sub in c("rnaseq","atacseq","integrated"))
    dir.create(file.path(out_dir,sub),recursive=TRUE,showWarnings=FALSE)

  # RNA -----------------------------------------------------------------------
  rna_raw <- load_rna_counts()
  if (!is.null(rna_raw)) {
    smap     <- RNA_SAMPLE_MAP |> filter(sample_id %in% colnames(rna_raw))
    log_rna  <- log2cpm(rna_raw)
    rna_res  <- run_rna_deseq2(rna_raw, smap, out_dir)
    plot_pca(log_rna, smap, out_dir, assay="rna")
    plot_deg_heatmap(log_rna, smap, rna_res, out_dir)
    plot_volcano(rna_res, out_dir, subdir="rnaseq")
    plot_gene_set_heatmaps(log_rna, smap, out_dir)
    plot_deg_overlap(rna_res, out_dir)
    if (!skip_gsea)
      tryCatch(run_gsea_analysis(rna_res, out_dir),
               error=function(e) message("GSEA skipped: ",conditionMessage(e)))
  } else {
    rna_res <- list()
  }

  # ATAC -----------------------------------------------------------------------
  atac_raw <- tryCatch(load_atac_counts(), error=function(e){message(conditionMessage(e));NULL})
  if (!is.null(atac_raw)) {
    log_atac <- log2cpm(atac_raw)
    plot_pca(log_atac, ATAC_SAMPLE_MAP, out_dir, assay="atac")
  }
  plot_insert_size(out_dir)
  plot_frip(out_dir)
  atac_dds  <- load_atac_dds()
  diff_atac <- run_atac_deseq2(atac_dds, out_dir)
  plot_atac_diff(diff_atac, out_dir)
  plot_volcano(diff_atac, out_dir, subdir="atacseq")
  plot_peak_annotation(diff_atac, out_dir)

  # Integrated -----------------------------------------------------------------
  if (length(rna_res)>0 && length(diff_atac)>0) {
    plot_rna_atac_correlation(rna_res, diff_atac, out_dir)
    plot_four_quadrant(rna_res, diff_atac, out_dir)
    plot_drug_specific(rna_res, diff_atac, out_dir)
  }

  message(strrep("=",60))
  message("Done. Results -> ", normalizePath(out_dir))
  message(strrep("=",60))
}

# =============================================================================
# STATUS CHECK on source()
# =============================================================================
message("\n--- PLK1 pipeline loaded ---")
key_files <- c(ATAC_COUNTS_FILE, ATAC_DESEQ2_RDS, ATAC_CONSENSUS_ANNOT,
               ATAC_INSERT_SIZE_MQC, RNA_COUNT_MATRIX)
for (f in key_files)
  message(if(file.exists(f)) "[OK]     " else "[MISSING]", " ", basename(f))
message("")
message("RNA count matrix missing? Run:")
message("  generate_rna_featurecounts_cmd('/path/to/genes.gtf')")
message("\nTo run full pipeline:")
message("  run_pipeline()")
message("----------------------------\n")
