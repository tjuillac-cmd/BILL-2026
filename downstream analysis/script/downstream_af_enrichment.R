#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  required_pkgs <- c("data.table", "dplyr", "tidyr", "stringr", "ggplot2", "readr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(sprintf("Missing required packages: %s", paste(missing_pkgs, collapse = ", ")))
  }
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    perm = 10000L,
    fdr = 0.05,
    outdir = "d:/vcf_visualiser/results2.1/tables",
    seed = 123L
  )
  if (length(args) == 0) return(out)

  for (a in args) {
    if (startsWith(a, "--perm=")) {
      out$perm <- as.integer(sub("^--perm=", "", a))
    } else if (startsWith(a, "--fdr=")) {
      out$fdr <- as.numeric(sub("^--fdr=", "", a))
    } else if (startsWith(a, "--outdir=")) {
      out$outdir <- sub("^--outdir=", "", a)
    } else if (startsWith(a, "--seed=")) {
      out$seed <- as.integer(sub("^--seed=", "", a))
    }
  }
  out
}

opts <- parse_args()
if (is.na(opts$perm) || opts$perm < 1) stop("--perm must be a positive integer")
if (is.na(opts$fdr) || opts$fdr <= 0 || opts$fdr >= 1) stop("--fdr must be in (0,1)")
if (is.na(opts$seed)) stop("--seed must be an integer")

input_files <- c(
  "d:/vcf_visualiser/results2.1/tables/P25-4-snp_details.tsv",
  "d:/vcf_visualiser/results2.1/tables/P25-8-snp_details.tsv",
  "d:/vcf_visualiser/results2.1/tables/P27-4-snp_details.tsv",
  "d:/vcf_visualiser/results2.1/tables/P27-8-snp_details.tsv"
)

required_cols <- c("Chrom", "Pos", "Type", "Gene", "BAM_AF", "Filter")

if (!dir.exists(opts$outdir)) {
  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)
}

load_one <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing input file: %s", path))
  dt <- fread(path, sep = "\t", na.strings = c("N/A", "NA", ""), showProgress = FALSE)

  missing <- setdiff(required_cols, colnames(dt))
  if (length(missing) > 0) {
    stop(sprintf("File %s is missing required columns: %s", path, paste(missing, collapse = ", ")))
  }

  base <- basename(path)
  m <- str_match(base, "^(P\\d+)-(4|8)-snp_details\\.tsv$")
  if (is.na(m[1, 1])) stop(sprintf("Unexpected filename format: %s", base))

  sample_id <- m[1, 2]
  timepoint <- m[1, 3]

  dt[, file := normalizePath(path, winslash = "/", mustWork = TRUE)]
  dt[, sample_id := sample_id]
  dt[, timepoint := timepoint]
  dt[, condition := factor(timepoint, levels = c("4", "8"))]
  dt[, variant_id := paste(sample_id, timepoint, Chrom, Pos, Type, sep = "|")]
  dt[, BAM_AF := suppressWarnings(as.numeric(BAM_AF))]
  dt[, Pos := suppressWarnings(as.integer(Pos))]
  dt
}

all_dt <- rbindlist(lapply(input_files, load_one), use.names = TRUE, fill = TRUE)

# Keep SNVs with valid AF and compute paired AF differences: triangle AF = AF(P27) - AF(P25).
af_dt <- all_dt %>%
  filter(Type == "SNV") %>%
  filter(!is.na(BAM_AF)) %>%
  filter(BAM_AF >= 0, BAM_AF <= 1)

paired_af <- af_dt %>%
  filter(sample_id %in% c("P25", "P27")) %>%
  group_by(Chrom, Pos, Ref, Alt, Type, condition, sample_id) %>%
  summarise(BAM_AF = mean(BAM_AF), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = sample_id, values_from = BAM_AF) %>%
  filter(!is.na(P25), !is.na(P27)) %>%
  mutate(af_diff = P27 - P25)

n4 <- sum(paired_af$condition == "4")
n8 <- sum(paired_af$condition == "8")
if (n4 == 0 || n8 == 0) stop("AF-difference dataset has zero paired rows in one or both conditions")

af_diff_4 <- paired_af %>% filter(condition == "4") %>% pull(af_diff)
af_diff_8 <- paired_af %>% filter(condition == "8") %>% pull(af_diff)

wilcox_res <- wilcox.test(af_diff ~ condition, data = paired_af, alternative = "two.sided", exact = FALSE)
ks_res <- suppressWarnings(ks.test(af_diff_4, af_diff_8, alternative = "two.sided"))

median_4 <- median(af_diff_4)
median_8 <- median(af_diff_8)
median_diff <- median_4 - median_8

# Effect sizes for two independent groups (condition 4 vs 8).
pair_cmp <- outer(af_diff_4, af_diff_8, "-")
wins <- sum(pair_cmp > 0)
losses <- sum(pair_cmp < 0)
ties <- sum(pair_cmp == 0)
u_stat <- wins + 0.5 * ties
rank_biserial <- (wins - losses) / (n4 * n8)
cliffs_delta <- rank_biserial

af_results <- data.frame(
  metric = c(
    "n_paired_diff_condition_4",
    "n_paired_diff_condition_8",
    "median_diff_condition_4",
    "median_diff_condition_8",
    "median_diff_4_minus_8",
    "wilcox_W",
    "wilcox_p",
    "ks_D",
    "ks_p",
    "rank_biserial",
    "cliffs_delta"
  ),
  value = c(
    n4,
    n8,
    median_4,
    median_8,
    median_diff,
    unname(wilcox_res$statistic),
    wilcox_res$p.value,
    unname(ks_res$statistic),
    ks_res$p.value,
    rank_biserial,
    cliffs_delta
  )
)

write_tsv(af_results, file.path(opts$outdir, "af_test_results.tsv"))

subtitle_txt <- sprintf(
  "Mann-Whitney p = %.3g | KS p = %.3g | froid -4 (n = %d), chaud -8 (n = %d)",
  wilcox_res$p.value, ks_res$p.value, n4, n8
)

p <- ggplot(paired_af, aes(x = af_diff, color = condition, fill = condition)) +
  geom_density(alpha = 0.25, adjust = 1.0) +
  geom_vline(
    data = data.frame(condition = factor(c("4", "8"), levels = c("4", "8")), med = c(median_4, median_8)),
    aes(xintercept = med, color = condition),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("4" = "#1b9e77", "8" = "#d95f02"), labels = c("4" = "froid traitement", "8" = "chaud traitement")) +
  scale_fill_manual(values = c("4" = "#1b9e77", "8" = "#d95f02"), labels = c("4" = "froid traitement", "8" = "chaud traitement")) +
  labs(
    title = "distribution de la différence d'AF (\u25b3AF) par condition",
    subtitle = subtitle_txt,
    x = "\u25b3AF = AF(P27) - AF(P25)",
    y = "Density",
    color = "Condition",
    fill = "Condition"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(opts$outdir, "af_density_plot.png"),
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)

p_violin <- ggplot(paired_af, aes(x = condition, y = af_diff, fill = condition, color = condition)) +
  geom_violin(alpha = 0.35, trim = FALSE, linewidth = 0.4) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.7, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.5, color = "gray40") +
  scale_color_manual(values = c("4" = "#1b9e77", "8" = "#d95f02"), labels = c("4" = "froid traitement", "8" = "chaud traitement")) +
  scale_fill_manual(values = c("4" = "#1b9e77", "8" = "#d95f02"), labels = c("4" = "froid traitement", "8" = "chaud traitement")) +
  scale_x_discrete(labels = c("4" = "froid traitement", "8" = "chaud traitement")) +
  labs(
    title = "violin + boxplot de la difference d'AF (\u25b3AF) par condition",
    subtitle = subtitle_txt,
    x = "Condition",
    y = "\u25b3AF = AF(P27) - AF(P25)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(
  filename = file.path(opts$outdir, "af_violin_boxplot.png"),
  plot = p_violin,
  width = 7,
  height = 5,
  dpi = 300
)

# Gene-level parsing from SNVs with valid AF.
gene_base <- af_dt %>%
  select(variant_id, condition, Pos, Gene) %>%
  mutate(Gene = ifelse(is.na(Gene), "", Gene))

gene_tokens <- gene_base %>%
  separate_rows(Gene, sep = ",") %>%
  mutate(gene_token = str_trim(Gene)) %>%
  filter(gene_token != "") %>%
  filter(str_detect(gene_token, "^CyHV3_ORF")) %>%
  distinct(variant_id, condition, Pos, gene_token)

if (nrow(gene_tokens) == 0) stop("No ORF tokens found after parsing Gene field")

obs_counts <- gene_tokens %>%
  count(gene_token, condition, name = "count") %>%
  tidyr::pivot_wider(names_from = condition, values_from = count, values_fill = 0, names_prefix = "count_") %>%
  mutate(
    count_4 = ifelse(is.na(count_4), 0L, as.integer(count_4)),
    count_8 = ifelse(is.na(count_8), 0L, as.integer(count_8)),
    count_total = count_4 + count_8,
    diff_4_minus_8 = count_4 - count_8
  )

gene_universe <- sort(unique(gene_tokens$gene_token))
obs_map <- setNames(obs_counts$diff_4_minus_8, obs_counts$gene_token)
obs_abs_map <- setNames(abs(obs_counts$diff_4_minus_8), obs_counts$gene_token)

# Permutation test on condition labels at variant level.
set.seed(opts$seed)
variant_to_cond <- af_dt %>% distinct(variant_id, condition)
variants <- variant_to_cond$variant_id
cond_labels <- as.character(variant_to_cond$condition)

variant_gene <- gene_tokens %>% distinct(variant_id, gene_token)

perm_extreme <- setNames(rep(0L, length(gene_universe)), gene_universe)

for (b in seq_len(opts$perm)) {
  shuffled <- sample(cond_labels, replace = FALSE)
  perm_df <- data.frame(variant_id = variants, perm_condition = shuffled, stringsAsFactors = FALSE)

  perm_counts <- variant_gene %>%
    left_join(perm_df, by = "variant_id") %>%
    count(gene_token, perm_condition, name = "count") %>%
    tidyr::pivot_wider(names_from = perm_condition, values_from = count, values_fill = 0, names_prefix = "count_") %>%
    mutate(
      count_4 = ifelse(is.na(count_4), 0L, as.integer(count_4)),
      count_8 = ifelse(is.na(count_8), 0L, as.integer(count_8)),
      diff_4_minus_8 = count_4 - count_8
    )

  perm_map <- setNames(perm_counts$diff_4_minus_8, perm_counts$gene_token)
  perm_all <- rep(0L, length(gene_universe)); names(perm_all) <- gene_universe
  perm_all[names(perm_map)] <- perm_map

  perm_extreme <- perm_extreme + as.integer(abs(perm_all) >= obs_abs_map)
}

perm_p <- (perm_extreme + 1) / (opts$perm + 1)
perm_tbl <- data.frame(
  gene_token = names(perm_p),
  perm_p = as.numeric(perm_p),
  stringsAsFactors = FALSE
)

# Poisson enrichment using observed per-gene covered positions as length proxy.
length_tbl <- gene_tokens %>%
  group_by(gene_token) %>%
  summarise(gene_span = n_distinct(Pos), .groups = "drop")

total_mut <- nrow(gene_tokens)
length_sum <- sum(length_tbl$gene_span)
if (length_sum <= 0) stop("Unable to compute gene span for Poisson model")

poisson_tbl <- obs_counts %>%
  left_join(length_tbl, by = "gene_token") %>%
  mutate(
    gene_span = ifelse(is.na(gene_span), 0, gene_span),
    lambda = total_mut * (gene_span / length_sum),
    poisson_p = ppois(q = count_total - 1, lambda = lambda, lower.tail = FALSE)
  ) %>%
  select(gene_token, gene_span, lambda, poisson_p)

full_tbl <- obs_counts %>%
  full_join(perm_tbl, by = "gene_token") %>%
  full_join(poisson_tbl, by = "gene_token") %>%
  mutate(
    count_4 = ifelse(is.na(count_4), 0L, as.integer(count_4)),
    count_8 = ifelse(is.na(count_8), 0L, as.integer(count_8)),
    count_total = ifelse(is.na(count_total), count_4 + count_8, as.integer(count_total)),
    diff_4_minus_8 = ifelse(is.na(diff_4_minus_8), count_4 - count_8, as.integer(diff_4_minus_8)),
    perm_p = ifelse(is.na(perm_p), 1, perm_p),
    poisson_p = ifelse(is.na(poisson_p), 1, poisson_p),
    perm_fdr = p.adjust(perm_p, method = "BH"),
    poisson_fdr = p.adjust(poisson_p, method = "BH"),
    significant_perm = perm_fdr <= opts$fdr,
    significant_poisson = poisson_fdr <= opts$fdr,
    significant_any = significant_perm | significant_poisson
  ) %>%
  arrange(perm_fdr, poisson_fdr, desc(count_total), gene_token)

ann_vcf_files <- c(
  "d:/vcf_visualiser/results2.1/tables/ann.vcf/P25-4.ann.vcf",
  "d:/vcf_visualiser/results2.1/tables/ann.vcf/P25-8.ann.vcf",
  "d:/vcf_visualiser/results2.1/tables/ann.vcf/P27-4.ann.vcf",
  "d:/vcf_visualiser/results2.1/tables/ann.vcf/P27-8.ann.vcf"
)

sig_genes <- full_tbl %>%
  filter(significant_any) %>%
  pull(gene_token) %>%
  unique()

# Include specific ORFs in addition to enriched genes for mutation-type and dN/dS analyses.
additional_orfs <- c("CyHV3_ORF64", "CyHV3_ORF89")
listed_genes <- unique(c(sig_genes, additional_orfs))

parse_ann_vcf <- function(path, sig_gene_set) {
  if (!file.exists(path)) return(data.frame())
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  if (length(lines) == 0) return(data.frame())

  sample_label <- sub("\\.ann\\.vcf$", "", basename(path))
  parsed <- vector("list", length(lines))
  out_idx <- 0L

  for (line in lines) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 8) next

    info <- fields[8]
    ann_m <- str_match(info, "(?:^|;)ANN=([^;]+)")
    if (is.na(ann_m[1, 2])) next

    ann_entries <- strsplit(ann_m[1, 2], ",", fixed = TRUE)[[1]]
    for (ann in ann_entries) {
      ann_fields <- strsplit(ann, "\\|")[[1]]
      if (length(ann_fields) < 4) next

      mutation_type_raw <- trimws(ann_fields[2])
      gene_name <- trimws(ann_fields[4])
      if (gene_name == "" || mutation_type_raw == "") next
      if (!(gene_name %in% sig_gene_set)) next

      mutation_types <- strsplit(mutation_type_raw, "&", fixed = TRUE)[[1]]
      mutation_types <- trimws(mutation_types)
      mutation_types <- mutation_types[mutation_types != ""]
      if (length(mutation_types) == 0) next

      for (mt in mutation_types) {
        out_idx <- out_idx + 1L
        parsed[[out_idx]] <- data.frame(
          sample_id = sample_label,
          chrom = fields[1],
          pos = suppressWarnings(as.integer(fields[2])),
          ref = fields[4],
          alt = fields[5],
          gene_token = gene_name,
          mutation_type = mt,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (out_idx == 0L) return(data.frame())
  bind_rows(parsed[seq_len(out_idx)])
}

mutation_type_raw <- bind_rows(lapply(ann_vcf_files, parse_ann_vcf, sig_gene_set = listed_genes))

if (nrow(mutation_type_raw) > 0) {
  mutation_type_tbl <- mutation_type_raw %>%
    distinct(sample_id, chrom, pos, ref, alt, gene_token, mutation_type) %>%
    count(gene_token, mutation_type, name = "mutation_count") %>%
    group_by(gene_token) %>%
    mutate(gene_total = sum(mutation_count)) %>%
    ungroup() %>%
    arrange(desc(gene_total), gene_token, mutation_type) %>%
    mutate(gene_token = factor(gene_token, levels = unique(gene_token)))

  p_mut_type <- ggplot(mutation_type_tbl, aes(x = gene_token, y = mutation_count, fill = mutation_type)) +
    geom_col(width = 0.82) +
    labs(
      title = "Mutation Type Distribution per Enriched Gene",
      subtitle = "FDR-significant genes + CyHV3_ORF64 + CyHV3_ORF89",
      x = "Gene",
      y = "Number of mutations",
      fill = "Mutation type"
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  mutation_type_tbl <- data.frame(
    gene_token = character(),
    mutation_type = character(),
    mutation_count = integer(),
    stringsAsFactors = FALSE
  )

  p_mut_type <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No mutation annotations found for significant genes") +
    theme_void() +
    labs(title = "Mutation Type Distribution per Enriched Gene")
}

write_tsv(mutation_type_tbl, file.path(opts$outdir, "mutation_type_distribution_enriched_genes.tsv"))

# Compare missense/synonymous ratios between cold (-4) and hot (-8) for enriched-gene variants.
ratio_source <- mutation_type_raw %>%
  distinct(sample_id, chrom, pos, ref, alt, gene_token, mutation_type) %>%
  mutate(
    condition_group = case_when(
      str_detect(sample_id, "-4$") ~ "froid",
      str_detect(sample_id, "-8$") ~ "chaud",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(condition_group)) %>%
  filter(mutation_type %in% c("missense_variant", "synonymous_variant"))

ratio_counts <- ratio_source %>%
  count(condition_group, mutation_type, name = "count") %>%
  tidyr::complete(
    condition_group = c("froid", "chaud"),
    mutation_type = c("missense_variant", "synonymous_variant"),
    fill = list(count = 0)
  )

ratio_matrix <- matrix(
  c(
    ratio_counts$count[ratio_counts$condition_group == "froid" & ratio_counts$mutation_type == "missense_variant"],
    ratio_counts$count[ratio_counts$condition_group == "froid" & ratio_counts$mutation_type == "synonymous_variant"],
    ratio_counts$count[ratio_counts$condition_group == "chaud" & ratio_counts$mutation_type == "missense_variant"],
    ratio_counts$count[ratio_counts$condition_group == "chaud" & ratio_counts$mutation_type == "synonymous_variant"]
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    c("froid", "chaud"),
    c("missense_variant", "synonymous_variant")
  )
)

ratio_froid <- if (sum(ratio_matrix["froid", ]) > 0) {
  ratio_matrix["froid", "missense_variant"] / sum(ratio_matrix["froid", ])
} else {
  NA_real_
}
ratio_chaud <- if (sum(ratio_matrix["chaud", ]) > 0) {
  ratio_matrix["chaud", "missense_variant"] / sum(ratio_matrix["chaud", ])
} else {
  NA_real_
}

if (sum(ratio_matrix) > 0) {
  chisq_res <- suppressWarnings(chisq.test(ratio_matrix, correct = FALSE))
  fisher_res <- fisher.test(ratio_matrix)

  mutation_ratio_test_tbl <- data.frame(
    metric = c(
      "froid_missense_count",
      "froid_synonymous_count",
      "chaud_missense_count",
      "chaud_synonymous_count",
      "froid_missense_ratio",
      "chaud_missense_ratio",
      "chisq_statistic",
      "chisq_p",
      "fisher_odds_ratio",
      "fisher_p"
    ),
    value = c(
      ratio_matrix["froid", "missense_variant"],
      ratio_matrix["froid", "synonymous_variant"],
      ratio_matrix["chaud", "missense_variant"],
      ratio_matrix["chaud", "synonymous_variant"],
      ratio_froid,
      ratio_chaud,
      unname(chisq_res$statistic),
      chisq_res$p.value,
      unname(fisher_res$estimate),
      fisher_res$p.value
    )
  )
} else {
  mutation_ratio_test_tbl <- data.frame(
    metric = c(
      "froid_missense_count",
      "froid_synonymous_count",
      "chaud_missense_count",
      "chaud_synonymous_count",
      "froid_missense_ratio",
      "chaud_missense_ratio",
      "chisq_statistic",
      "chisq_p",
      "fisher_odds_ratio",
      "fisher_p"
    ),
    value = c(0, 0, 0, 0, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
  )
}

write_tsv(mutation_ratio_test_tbl, file.path(opts$outdir, "mutation_ratio_missense_synonymous_tests.tsv"))

fisher_p_heat <- mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "fisher_p"]
fisher_or_heat <- mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "fisher_odds_ratio"]
fisher_heatmap_caption <- if (!is.na(fisher_p_heat) && !is.na(fisher_or_heat)) {
  sprintf(
    "2x2 Fisher exact test (missense vs synonymous by froid/chaud): p=%.3g, OR=%.3f",
    fisher_p_heat, fisher_or_heat
  )
} else {
  "2x2 Fisher exact test unavailable"
}

# dN/dS proxy per listed gene from missense/synonymous counts (with 0.5 pseudocount).
dnds_source <- mutation_type_raw %>%
  distinct(sample_id, chrom, pos, ref, alt, gene_token, mutation_type) %>%
  mutate(
    condition_group = case_when(
      str_detect(sample_id, "-4$") ~ "froid",
      str_detect(sample_id, "-8$") ~ "chaud",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(condition_group)) %>%
  filter(mutation_type %in% c("missense_variant", "synonymous_variant"))

if (nrow(dnds_source) > 0) {
  dnds_by_cond <- dnds_source %>%
    count(gene_token, condition_group, mutation_type, name = "n") %>%
    tidyr::complete(
      gene_token = listed_genes,
      condition_group = c("froid", "chaud"),
      mutation_type = c("missense_variant", "synonymous_variant"),
      fill = list(n = 0)
    )

  dnds_overall <- dnds_source %>%
    count(gene_token, mutation_type, name = "n") %>%
    mutate(condition_group = "overall") %>%
    tidyr::complete(
      gene_token = listed_genes,
      condition_group = "overall",
      mutation_type = c("missense_variant", "synonymous_variant"),
      fill = list(n = 0)
    )

  dnds_long <- bind_rows(dnds_by_cond, dnds_overall)

  dnds_tbl <- dnds_long %>%
    tidyr::pivot_wider(names_from = mutation_type, values_from = n, values_fill = 0) %>%
    mutate(
      dN = missense_variant,
      dS = synonymous_variant,
      dnds = (dN + 0.5) / (dS + 0.5),
      log2_dnds = log2(dnds)
    ) %>%
    select(gene_token, condition_group, dN, dS, dnds, log2_dnds)

  dnds_tbl <- dnds_tbl %>%
    mutate(
      condition_group = factor(condition_group, levels = c("froid", "chaud", "overall")),
      gene_token = factor(gene_token, levels = rev(listed_genes))
    ) %>%
    arrange(gene_token, condition_group)

  p_dnds <- ggplot(dnds_tbl, aes(x = condition_group, y = gene_token, fill = log2_dnds)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.2f", dnds)), size = 3.1) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "#f7f7f7",
      high = "#b2182b",
      midpoint = 0,
      name = "log2(dN/dS)"
    ) +
    labs(
      title = "dN/dS per Gene (listed genes)",
      subtitle = "Listed genes = FDR-significant genes + CyHV3_ORF64 + CyHV3_ORF89",
      x = "Condition",
      y = "Gene",
      caption = fisher_heatmap_caption
    ) +
    theme_bw(base_size = 12)
} else {
  dnds_tbl <- data.frame(
    gene_token = character(),
    condition_group = character(),
    dN = integer(),
    dS = integer(),
    dnds = numeric(),
    log2_dnds = numeric(),
    stringsAsFactors = FALSE
  )

  p_dnds <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No missense/synonymous annotations for listed genes") +
    theme_void() +
    labs(
      title = "dN/dS per Gene (listed genes)",
      caption = fisher_heatmap_caption
    )
}

write_tsv(dnds_tbl, file.path(opts$outdir, "dnds_per_gene_listed_genes.tsv"))
ggsave(
  filename = file.path(opts$outdir, "dnds_heatmap_listed_genes.png"),
  plot = p_dnds,
  width = 8,
  height = 5.5,
  dpi = 300
)

if (nrow(mutation_type_raw) > 0) {
  froid_ratio <- mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "froid_missense_ratio"]
  chaud_ratio <- mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "chaud_missense_ratio"]
  chi_p <- mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "chisq_p"]
  fisher_p <- mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "fisher_p"]
  fisher_or <- mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "fisher_odds_ratio"]

  metric_caption <- sprintf(
    "Missense ratio: froid=%.3f, chaud=%.3f | Chi-square p=%.3g | Fisher p=%.3g | Fisher OR=%.3f",
    froid_ratio, chaud_ratio, chi_p, fisher_p, fisher_or
  )

  p_mut_type <- p_mut_type +
    geom_text(
      data = mutation_type_tbl,
      aes(label = mutation_count),
      position = position_stack(vjust = 0.5),
      size = 3,
      color = "white",
      show.legend = FALSE
    ) +
    labs(caption = metric_caption)
}

ggsave(
  filename = file.path(opts$outdir, "mutation_type_stacked_barplot_enriched_genes.png"),
  plot = p_mut_type,
  width = 10,
  height = 6,
  dpi = 300
)


gene_bar_tbl <- full_tbl %>%
  mutate(
    enrichment_score = -log10(pmin(perm_fdr, poisson_fdr)),
    enrichment_score = ifelse(is.finite(enrichment_score), enrichment_score, 0),
    condition_label = case_when(
      diff_4_minus_8 > 0 ~ "traitement froid",
      diff_4_minus_8 < 0 ~ "traitement chaud",
      TRUE ~ "sans difference"
    )
  ) %>%
  arrange(desc(enrichment_score), desc(count_total)) %>%
  slice_head(n = 20) %>%
  mutate(gene_token = factor(gene_token, levels = rev(gene_token)))

p_gene <- ggplot(gene_bar_tbl, aes(x = gene_token, y = enrichment_score, fill = condition_label)) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = c("traitement froid" = "#1b9e77", "traitement chaud" = "#d95f02", "sans difference" = "#757575")
  ) +
  labs(
    title = "enrichissement g\u00e9nique selon la condition",
    subtitle = "Top 20 g\u00e8nes selon -log10(FDR min entre permutation et Poisson)",
    x = "G\u00e8ne (ORF)",
    y = "Score d'enrichissement (-log10 FDR)",
    fill = "Condition enrichie"
  ) +
  coord_flip() +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(opts$outdir, "gene_enrichment_barplot.png"),
  plot = p_gene,
  width = 9,
  height = 6,
  dpi = 300
)
write_tsv(full_tbl, file.path(opts$outdir, "gene_enrichment_full.tsv"))
write_tsv(full_tbl %>% filter(significant_any), file.path(opts$outdir, "significant_genes.tsv"))

cat("Paired AF-difference rows by condition:\n")
cat(sprintf("  cold (-4): %d\n", n4))
cat(sprintf("  hot (-8): %d\n", n8))
cat(sprintf("Median △AF cold (-4): %.6g\n", median_4))
cat(sprintf("Median △AF hot (-8): %.6g\n", median_8))
cat(sprintf("Median difference (cold-hot): %.6g\n", median_diff))
cat(sprintf("Wilcoxon p-value: %.6g\n", wilcox_res$p.value))
cat(sprintf("KS p-value: %.6g\n", ks_res$p.value))
cat(sprintf("Rank-biserial effect size: %.6g\n", rank_biserial))
cat(sprintf("Cliff's delta: %.6g\n", cliffs_delta))
cat(sprintf("Significant genes (perm FDR <= %.3f): %d\n", opts$fdr, sum(full_tbl$significant_perm)))
cat(sprintf("Significant genes (poisson FDR <= %.3f): %d\n", opts$fdr, sum(full_tbl$significant_poisson)))
cat(sprintf("Significant genes (any method): %d\n", sum(full_tbl$significant_any)))
cat(sprintf(
  "Missense/synonymous counts (froid): %d / %d\n",
  mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "froid_missense_count"],
  mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "froid_synonymous_count"]
))
cat(sprintf(
  "Missense/synonymous counts (chaud): %d / %d\n",
  mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "chaud_missense_count"],
  mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "chaud_synonymous_count"]
))
cat(sprintf("Chi-square p-value: %.6g\n", mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "chisq_p"]))
cat(sprintf("Fisher exact p-value: %.6g\n", mutation_ratio_test_tbl$value[mutation_ratio_test_tbl$metric == "fisher_p"]))
cat(sprintf("dN/dS rows (gene x condition): %d\n", nrow(dnds_tbl)))
cat(sprintf("Outputs written to: %s\n", normalizePath(opts$outdir, winslash = "/", mustWork = TRUE)))
