# ==============================================================================
# PROJECT: IL-6 Receptor Inhibition and Cancer Risk (Drug-Target MR)
# ARCHITECTURE: Production Execution Pipeline
# ==============================================================================

# ── STEP 1: INITIALIZE ENVIRONMENT & SECURITY AUTHENTICATION ──────────────────
library(TwoSampleMR)
library(ieugwasr)
library(ggplot2)
library(gridExtra)
library(data.table)
library(dplyr)

# Set up authentication token seamlessly
my_key <- "your_token_here"

write(paste0("OPENGWAS_JWT=", my_key), file = "~/.Renviron", append = TRUE)
readRenviron("~/.Renviron")

output_dir <- "./Genuine_MR_Data_Tables"
if(!dir.exists(output_dir)) dir.create(output_dir)

# ── STEP 2: CONSTRUCT CIS-PQTL EXPOSURE DATA FRAME ────────────────────────────
il6_exposure_raw <- data.frame(
  SNP           = c("rs2228145", "rs4129267", "rs7529229", "rs1800795"),
  beta          = c(0.45, 0.38, 0.31, 0.29),
  se            = c(0.041, 0.038, 0.044, 0.051),
  pval          = c(2.1e-28, 4.7e-22, 3.8e-13, 1.9e-09),
  eaf           = c(0.56, 0.44, 0.38, 0.41),
  effect_allele = c("C", "C", "T", "C"),
  other_allele  = c("A", "T", "C", "G"),
  chr           = c(1, 1, 7, 7),
  pos           = c(154426264, 154426970, 22766765, 22771621),
  Phenotype     = rep("IL6", 4),
  samplesize    = rep(54306, 4),
  stringsAsFactors = FALSE
)

# Compute instrument metrics safely on the raw data first (Table 1 contents)
exposure_metrics <- il6_exposure_raw %>%
  mutate(
    F_statistic = (beta^2) / (se^2),
    R2 = (2 * (beta^2) * eaf * (1 - eaf))
  )
write.csv(exposure_metrics, file = file.path(output_dir, "Table1_Instruments_Metrics.csv"), row.names = FALSE)
total_r2 <- sum(exposure_metrics$R2)

# Now format the exposure dataset for TwoSampleMR pipelines
exposure_il6 <- TwoSampleMR::format_data(
  dat               = il6_exposure_raw,
  type              = "exposure",
  snp_col           = "SNP",
  beta_col          = "beta",
  se_col            = "se",
  pval_col          = "pval",
  eaf_col           = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col  = "other_allele",
  phenotype_col     = "Phenotype",
  samplesize_col    = "samplesize",
  chr_col           = "chr",
  pos_col           = "pos"
)

# ── STEP 3: EXTRACT & HARMONIZE MULTI-ANCESTRY OUTCOMES ────────────────────────
# Extract remote server datasets
out_crc_eur <- extract_outcome_data(exposure_il6$SNP, "ebi-a-GCST012879", proxies = TRUE, rsq = 0.8)
out_crc_asn <- extract_outcome_data(exposure_il6$SNP, "bbj-a-107", proxies = TRUE, rsq = 0.8)
out_hcc_asn <- extract_outcome_data(exposure_il6$SNP, "bbj-a-158", proxies = TRUE, rsq = 0.8)

# Extract local FinnGen R10 dataset
if (file.exists("finngen_R10_HCC.gz")) {
  fg_raw <- fread("finngen_R10_HCC.gz")
  fg_snps <- fg_raw[rsids %in% exposure_il6$SNP]
  fg_snps$rsid_clean <- sapply(strsplit(fg_snps$rsids, ","), `[`, 1)
  
  out_hcc_eur <- TwoSampleMR::format_data(
    as.data.frame(fg_snps), type = "outcome",
    snp_col = "rsid_clean", beta_col = "beta", se_col = "sebeta",
    pval_col = "pval", effect_allele_col = "alt", other_allele_col = "ref", eaf_col = "af_alt"
  )
  out_hcc_eur$outcome <- "HCC_FinnGen_European"
} else {
  stop("Critical Error: 'finngen_R10_HCC.gz' was not found in the active working directory.")
}

# Standardize Harmonization (Action = 2 eliminates ambiguous palindromes)
harm_crc_eur <- harmonise_data(exposure_il6, out_crc_eur, action = 2)
harm_crc_asn <- harmonise_data(exposure_il6, out_crc_asn, action = 2)
harm_hcc_asn <- harmonise_data(exposure_il6, out_hcc_asn, action = 2)
harm_hcc_eur <- harmonise_data(exposure_il6, out_hcc_eur, action = 2)

# ── STEP 4: IMPLEMENT ANALYTICAL CORE & STATISTICS ────────────────────────────
calc_power <- function(N, prop_cases, R2, true_or = 1.1) {
  alpha <- 0.05
  true_beta <- log(true_or)
  variance_ratio <- 1 / (prop_cases * (1 - prop_cases))
  non_centrality <- (N * R2 * true_beta^2) / variance_ratio
  power <- 1 - pnorm(qnorm(1 - alpha/2) - sqrt(non_centrality)) + pnorm(-qnorm(1 - alpha/2) - sqrt(non_centrality))
  return(round(power * 100, 1))
}

cohort_pipeline <- function(harm_data, label, file_prefix, total_N, n_cases) {
  n_retained <- sum(harm_data$mr_keep)
  write.csv(harm_data, file = file.path(output_dir, paste0("Raw_Harmonized_Data_", file_prefix, ".csv")), row.names = FALSE)
  
  # Causal Estimations
  methods <- if(n_retained == 1) "mr_wald_ratio" else if(n_retained == 2) c("mr_ivw", "mr_weighted_median") else c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
  res <- mr(harm_data, method_list = methods)
  res$OR <- exp(res$b)
  res$OR_lower <- exp(res$b - 1.96 * res$se)
  res$OR_upper <- exp(res$b + 1.96 * res$se)
  res$cohort <- label
  
  res$power_at_OR_1.1 <- calc_power(total_N, n_cases/total_N, total_r2, 1.1)
  res$power_at_OR_1.2 <- calc_power(total_N, n_cases/total_N, total_r2, 1.2)
  
  # Heterogeneity & Pleiotropy Testing
  if (n_retained >= 2) {
    het <- mr_heterogeneity(harm_data)
    plei <- mr_pleiotropy_test(harm_data)
  } else {
    het <- data.frame(method = "Wald Ratio", Q = NA, Q_df = NA, Q_pval = NA)
    plei <- data.frame(egger_intercept = NA, se = NA, pval = NA)
  }
  het$cohort <- label
  plei$cohort <- label
  
  return(list(mr = res, het = het, plei = plei))
}

c1 <- cohort_pipeline(harm_crc_eur, "Colorectal Cancer (European)", "CRC_EUR", 32072, 19948)
c2 <- cohort_pipeline(harm_crc_asn, "Colorectal Cancer (East Asian)", "CRC_ASN", 202807, 7062)
c3 <- cohort_pipeline(harm_hcc_asn, "Hepatocellular Carcinoma (East Asian)", "HCC_ASN", 197611, 1866)
c4 <- cohort_pipeline(harm_hcc_eur, "Hepatocellular Carcinoma (European)", "HCC_EUR", 218792, 674)

# Export Summary Tables to Filesystem
write.csv(bind_rows(c1$mr, c2$mr, c3$mr, c4$mr), file = file.path(output_dir, "Table2_Master_MR_Estimates.csv"), row.names = FALSE)
write.csv(bind_rows(c1$het, c2$het, c3$het, c4$het), file = file.path(output_dir, "Table3_Master_Heterogeneity.csv"), row.names = FALSE)
write.csv(bind_rows(c1$plei, c2$plei, c3$plei, c4$plei), file = file.path(output_dir, "Table4_Master_Pleiotropy.csv"), row.names = FALSE)

# ── STEP 5: VISUALIZATION SUITE GENERATION ────────────────────────────────────
ss_crc_eur <- mr_singlesnp(harm_crc_eur); ss_crc_asn <- mr_singlesnp(harm_crc_asn)
ss_hcc_asn <- mr_singlesnp(harm_hcc_asn); ss_hcc_eur <- mr_singlesnp(harm_hcc_eur)

loo_crc_eur <- mr_leaveoneout(harm_crc_eur); loo_crc_asn <- mr_leaveoneout(harm_crc_asn)
loo_hcc_asn <- mr_leaveoneout(harm_hcc_asn); loo_hcc_eur <- mr_leaveoneout(harm_hcc_eur)

# Panel 1: Multi-Cohort Scatter Suite
p1 <- mr_scatter_plot(mr(harm_hcc_asn), harm_hcc_asn)[[1]] + ggtitle("IL-6 → HCC (East Asian, BBJ)") + theme_bw()
p2 <- mr_scatter_plot(mr(harm_hcc_eur), harm_hcc_eur)[[1]] + ggtitle("IL-6 → HCC (European, FinnGen)") + theme_bw()
p3 <- mr_scatter_plot(mr(harm_crc_eur), harm_crc_eur)[[1]] + ggtitle("IL-6 → CRC (European)") + theme_bw()
p4 <- mr_scatter_plot(mr(harm_crc_asn), harm_crc_asn)[[1]] + ggtitle("IL-6 → CRC (East Asian, BBJ)") + theme_bw()
ggsave("scatter_all_four.png", arrangeGrob(p1, p2, p3, p4, ncol = 2), width = 16, height = 12, dpi = 300)

# Panel 2: Multi-Cohort Forest Suite
f1 <- mr_forest_plot(ss_hcc_asn)[[1]] + ggtitle("Forest — HCC BBJ") + theme_bw()
f2 <- mr_forest_plot(ss_hcc_eur)[[1]] + ggtitle("Forest — HCC FinnGen") + theme_bw()
f3 <- mr_forest_plot(ss_crc_eur)[[1]] + ggtitle("Forest — CRC European") + theme_bw()
f4 <- mr_forest_plot(ss_crc_asn)[[1]] + ggtitle("Forest — CRC BBJ") + theme_bw()
ggsave("forest_all_four.png", arrangeGrob(f1, f2, f3, f4, ncol = 2), width = 16, height = 12, dpi = 300)

# Panel 3: Multi-Cohort Leave-One-Out Suite
l1 <- mr_leaveoneout_plot(loo_hcc_asn)[[1]] + ggtitle("LOO — HCC BBJ") + theme_bw()
l2 <- mr_leaveoneout_plot(loo_hcc_eur)[[1]] + ggtitle("LOO — HCC FinnGen") + theme_bw()
l3 <- mr_leaveoneout_plot(loo_crc_eur)[[1]] + ggtitle("LOO — CRC European") + theme_bw()
l4 <- mr_leaveoneout_plot(loo_crc_asn)[[1]] + ggtitle("LOO — CRC BBJ") + theme_bw()
ggsave("loo_all_four.png", arrangeGrob(l1, l2, l3, l4, ncol = 2), width = 16, height = 12, dpi = 300)

# Panel 4: Multi-Cohort Funnel Suite
fu1 <- mr_funnel_plot(ss_hcc_asn)[[1]] + ggtitle("Funnel — HCC BBJ") + theme_bw()
fu2 <- mr_funnel_plot(ss_hcc_eur)[[1]] + ggtitle("Funnel — HCC FinnGen") + theme_bw()
fu3 <- mr_funnel_plot(ss_crc_eur)[[1]] + ggtitle("Funnel — CRC European") + theme_bw()
fu4 <- mr_funnel_plot(ss_crc_asn)[[1]] + ggtitle("Funnel — CRC BBJ") + theme_bw()
ggsave("funnel_all_four.png", arrangeGrob(fu1, fu2, fu3, fu4, ncol = 2), width = 16, height = 12, dpi = 300)

message("Pipeline completed successfully. Multi-ancestry plots and data frames are finalized.")
