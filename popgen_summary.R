# Local sweep flags based on nucleotide diversity, Tajima's D and FST.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

ploidy <- "ta"  # da or ta; statistic files are selected explicitly below
candidate_file <- sprintf("results/candidate_coordinates/%s_candidate_coordinates.tsv", ploidy)
ploidy_x <- if (ploidy == "da") "2x" else "4x"

pi_file <- sprintf("data/aa_bt_1kb_mask_%s_pi.txt", ploidy_x)
taj_file <- sprintf("data/aa_bt_1kb_mask_%s_tajima_d.txt", ploidy_x)
fst_file <- sprintf("data/aa_bt_1kb_mask_%s_fst.txt", ploidy_x)
out_dir <- file.path("results/local_sweeps", ploidy_x)

pad_bp <- 50000
gene_buffer_bp <- 2000
min_windows_gene <- 3
p_low <- 0.05
p_high <- 0.95
same_ploidy_only <- TRUE
exclude_labels <- character(0)
require_fst <- TRUE
taj_abs_cutoff <- 2
use_any_bin_for_taj <- TRUE

stopifnot(ploidy %in% c("da", "ta"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
genes <- read_tsv(candidate_file, show_col_types = FALSE)
if (!nrow(genes)) stop("No candidate genes found.")
if (anyNA(genes) || any(genes$gene_start < 1 | genes$gene_end < genes$gene_start)) {
  stop("Invalid candidate coordinates.")
}

pi_df <- read.table(pi_file, header = TRUE) %>%
  filter(no_sites >= 10) %>%
  rename(bin_start = window_pos_1, pi = avg_pi, lineage = pop) %>%
  mutate(bin_end = bin_start + 999L)
taj_df <- read.table(taj_file, header = TRUE) %>%
  filter(no_sites >= 10) %>%
  rename(bin_start = window_pos_1, taj = tajima_d, lineage = pop) %>%
  mutate(bin_end = bin_start + 999L)
fst_df <- read.table(fst_file, header = TRUE) %>%
  filter(no_snps >= 10) %>%
  rename(bin_start = window_pos_1, fst = avg_hudson_fst) %>%
  mutate(bin_end = bin_start + 999L)

if (same_ploidy_only) {
  fst_df <- fst_df %>%
    filter(!is.na(str_extract(pop1, "_[24]x$")),
           str_extract(pop1, "_[24]x$") == str_extract(pop2, "_[24]x$"))
}
if (length(exclude_labels)) {
  fst_df <- fst_df %>% filter(!paste(pop1, "vs", pop2) %in% exclude_labels)
}

result_rows <- vector("list", nrow(genes))

for (i in seq_len(nrow(genes))) {
  chrom <- genes$chromosome[i]
  gstart <- genes$gene_start[i]
  gend <- genes$gene_end[i]
  win_start <- max(0, gstart - pad_bp)
  win_end <- gend + pad_bp
  excl_start <- gstart - gene_buffer_bp
  excl_end <- gend + gene_buffer_bp
  
  pi_local <- pi_df %>%
    filter(chromosome == chrom, bin_start <= win_end, bin_end >= win_start)
  taj_local <- taj_df %>%
    filter(chromosome == chrom, bin_start <= win_end, bin_end >= win_start)
  fst_local <- fst_df %>%
    filter(chromosome == chrom, bin_start <= win_end, bin_end >= win_start)
  
  pi_gene <- pi_local %>% filter(bin_start <= gend, bin_end >= gstart)
  taj_gene <- taj_local %>% filter(bin_start <= gend, bin_end >= gstart)
  fst_gene <- fst_local %>% filter(bin_start <= gend, bin_end >= gstart)
  pi_bg <- pi_local %>% filter(!(bin_start <= excl_end & bin_end >= excl_start))
  taj_bg <- taj_local %>% filter(!(bin_start <= excl_end & bin_end >= excl_start))
  fst_bg <- fst_local %>% filter(!(bin_start <= excl_end & bin_end >= excl_start))
  
  pi_summary <- pi_gene %>%
    group_by(lineage) %>%
    summarise(n_bins = n(), mean_pi = mean(pi, na.rm = TRUE), .groups = "drop") %>%
    filter(n_bins >= min_windows_gene)
  taj_summary <- taj_gene %>%
    group_by(lineage) %>%
    summarise(n_bins = n(), mean_taj = mean(taj, na.rm = TRUE), .groups = "drop") %>%
    filter(n_bins >= min_windows_gene)
  fst_summary <- fst_gene %>%
    group_by(pop1, pop2) %>%
    summarise(n_bins = n(), mean_fst = mean(fst, na.rm = TRUE), .groups = "drop") %>%
    filter(n_bins >= min_windows_gene)
  
  pi_thresholds <- pi_bg %>%
    group_by(lineage) %>%
    summarise(pi_p05 = quantile(pi, p_low, na.rm = TRUE), .groups = "drop")
  taj_thresholds <- taj_bg %>%
    group_by(lineage) %>%
    summarise(taj_p05 = quantile(taj, p_low, na.rm = TRUE), .groups = "drop")
  fst_thresholds <- fst_bg %>%
    group_by(pop1, pop2) %>%
    summarise(fst_p95 = quantile(fst, p_high, na.rm = TRUE), .groups = "drop")
  
  low_pi <- pi_summary %>%
    inner_join(pi_thresholds, by = "lineage") %>%
    filter(mean_pi <= pi_p05) %>%
    pull(lineage)
  low_taj <- taj_summary %>%
    inner_join(taj_thresholds, by = "lineage") %>%
    filter(mean_taj <= taj_p05) %>%
    pull(lineage)
  high_fst <- fst_summary %>%
    inner_join(fst_thresholds, by = c("pop1", "pop2")) %>%
    filter(mean_fst >= fst_p95)
  
  # The absolute-D rule is independent of the minimum-bin requirement.
  taj_abs <- taj_gene %>%
    group_by(lineage) %>%
    summarise(hit = if (use_any_bin_for_taj) {
      any(abs(taj) >= taj_abs_cutoff, na.rm = TRUE)
    } else {
      abs(mean(taj, na.rm = TRUE)) >= taj_abs_cutoff
    }, .groups = "drop") %>%
    filter(hit) %>%
    pull(lineage)
  
  low_div <- unique(c(low_pi, low_taj, taj_abs))
  matched_fst <- any(high_fst$pop1 %in% low_div | high_fst$pop2 %in% low_div)
  
  result_rows[[i]] <- tibble(
    gene = genes$gene_name[i], chrom = chrom,
    n_pi_bins_gene = sum(pi_summary$n_bins),
    n_taj_bins_gene = sum(taj_summary$n_bins),
    n_fst_bins_gene = sum(fst_summary$n_bins),
    mean_pi = mean(pi_gene$pi, na.rm = TRUE),
    mean_tajimaD = mean(taj_gene$taj, na.rm = TRUE),
    mean_fst = mean(fst_gene$fst, na.rm = TRUE),
    any_low_div_local = length(low_div) > 0,
    any_high_fst_local = matched_fst,
    sweep_local = length(low_div) > 0 && (!require_fst || matched_fst),
    low_div_lineages = paste(low_div, collapse = ";"),
    taj_abs2_lineages = paste(taj_abs, collapse = ";"),
    pairs_high_fst = if (nrow(high_fst)) {
      paste(paste(high_fst$pop1, "vs", high_fst$pop2), collapse = ";")
    } else ""
  )
}

results <- bind_rows(result_rows)

summary_table <- results %>%
  summarise(total_genes = n(), sweep_genes = sum(sweep_local),
            perc_sweep = round(100 * sweep_genes / total_genes, 1))
signal_counts <- results %>%
  summarise(n_low_div_any = sum(any_low_div_local),
            n_sweep_with_fst = sum(any_low_div_local & any_high_fst_local),
            n_taj_abs2_lineage_hits = sum(nzchar(taj_abs2_lineages)))
lineage_counts <- results %>%
  select(low_div_lineages) %>%
  tidyr::separate_rows(low_div_lineages, sep = ";") %>%
  filter(nzchar(low_div_lineages)) %>%
  count(lineage = low_div_lineages, sort = TRUE)

write_tsv(results, file.path(out_dir, paste0(ploidy, "_local_sweeps.tsv")))
write_tsv(summary_table, file.path(out_dir, paste0(ploidy, "_sweep_summary.tsv")))
write_tsv(signal_counts, file.path(out_dir, paste0(ploidy, "_signal_counts.tsv")))
write_tsv(lineage_counts, file.path(out_dir, paste0(ploidy, "_lineage_counts.tsv")))
print(summary_table)