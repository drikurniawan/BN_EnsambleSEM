library(lavaan)
library(scales)
library(dplyr)
library(semTools)
library(bestNormalize)
library(scales)
library(tidyverse)

setwd("D:/2026/disertasi/paper1")
# Step 1: Load your survey data
data_2024_raw <- read.csv("data_survey_2024.csv")
data_2025_raw <- read.csv("data_survey_2025.csv")

data_2024_raw <- data_2024_raw[,-c(1,2,8)]
data_2025_raw <- data_2025_raw[,-c(1,2,3,6,7,11,12)]

summary(data_2024_raw)
summary(data_2025_raw)
#===============================================================================
# imputation
# Install mice package if not already installed
if (!require(mice)) {
  install.packages("mice")
  library(mice)
}

# Perform MICE imputation
# m = number of multiple imputations
# method = imputation method for each variable (can be left default)
# maxit = number of iterations

df24 = data_2024_raw
df25 = data_2025_raw
imputed24 <- mice(df24, m = 5, method = "pmm", maxit = 5, seed = 123)
imputed25 <- mice(df25, m = 5, method = "pmm", maxit = 5, seed = 123)
# View imputation summary
summary(imputed24)
summary(imputed25)
# Get one completed dataset (first imputation)
completed_data24 <- complete(imputed24, 1)
completed_data25 <- complete(imputed25, 1)
cat("\nCompleted dataset after imputation:\n")

#===============================================================================
survey_data_2024 <- completed_data24
survey_data_2025 <- completed_data25
# Step 2: Detect straightliners of all 5's
library(dplyr)
library(tidyr)

# 1. Define your 35 target indicators
indicator_cols <- c(
  paste0("CTA", 1:5),
  paste0("KBB", 1:5),
  paste0("SPI", 1:5),
  paste0("RBB", 1:5),
  paste0("KAB", 1:5),
  paste0("NIAT", 1:10)
)

ITEMS <- list(
  CTA = paste0("CTA",  1:5),
  KBB = paste0("KBB",  1:5),
  SPI = paste0("SPI",  1:5),
  RBB = paste0("RBB",  1:5),
  KAB = paste0("KAB",  1:5),
  NIAT = paste0("NIAT", 1:10)
)

# 2. remove streightliners

remove_straightliners <- function(data, indicator_cols) {
  cleaned_data <- data %>%
    rowwise() %>%
    # Menghitung standar deviasi di seluruh indikator yang ditentukan
    mutate(respondent_sd = sd(c_across(all_of(indicator_cols)), na.rm = TRUE)) %>%
    ungroup() %>%
    # Hanya mempertahankan baris dengan varians (SD > 0)
    filter(respondent_sd > 0) %>%
    # Menghapus kolom SD sementara
    select(-respondent_sd)
  
  return(cleaned_data)
}
cleaned_data_2024 <- remove_straightliners(survey_data_2024,indicator_cols)
cleaned_data_2025 <- remove_straightliners(survey_data_2025,indicator_cols)

# 3. Verify the final sample size
nrow(cleaned_data_2024)
nrow(cleaned_data_2025)
# install.packages(c("lavaan", "semTools"))
library(lavaan)
library(semTools)
model_definition_v0 <- "
  NIAT =~ NIAT1 + NIAT2 + NIAT3 + NIAT4 + NIAT5 + NIAT6 + NIAT7 + NIAT8 + NIAT9 + NIAT10
  CTA =~ CTA1 + CTA2 + CTA3 + CTA4 + CTA5
  KBB =~ KBB1 + KBB2 + KBB3 + KBB4 + KBB5
  SPI =~ SPI1 + SPI2 + SPI3 + SPI4 + SPI5
  RBB =~ RBB1 + RBB2 + RBB3 + RBB4 + RBB5
  KAB =~ KAB1 + KAB2 + KAB3 + KAB4 + KAB5
"

# 1. Fit a basic confirmatory factor analysis (CFA) model

cfa_fit <- cfa(model_definition_v0, data = cleaned_data_2024, estimator = "MLR", std.lv    = TRUE)

fit_measures <- fitMeasures(cfa_fit, 
                            c("chisq", "df", "pvalue", "cfi", "tli", 
                              "rmsea", "srmr"))
sink("BN_SEM_141026/fit_model_v0_full.txt")
summary(cfa_fit, fit.measures = TRUE, standardized = TRUE)
print(fit_measures)
# ── A6. Reliability (Cronbach alpha) ──────────────────────────────────────────
cat("--- Reliability Estimates ---\n")
reliability_df <- do.call(rbind, lapply(names(ITEMS), function(cn) {
  items_present <- intersect(ITEMS[[cn]], names(cleaned_data_2024))
  if (length(items_present) < 2) return(NULL)
  d_sub <- data_2024_raw[, items_present]
  d_sub <- d_sub[complete.cases(d_sub), ]
  alpha_val <- tryCatch(psych::alpha(d_sub)$total$raw_alpha, error = function(e) NA)
  data.frame(Construct = cn,
             n_items   = length(items_present),
             Alpha     = round(alpha_val, 3),
             Threshold = "≥ .70")
}))
print(reliability_df)

# ── AVE dan CR ────────────────────────────────────────────────
CONSTRUCTS  <- names(ITEMS)
compute_AVE_CR <- function(fit, constructs) {
  params <- parameterEstimates(fit, standardized = TRUE) %>%
    filter(op == "=~")
  do.call(rbind, lapply(constructs, function(cn) {
    lam <- params[params$lhs == cn, "std.all"]
    AVE <- mean(lam^2)
    CR  <- sum(lam)^2 / (sum(lam)^2 + sum(1 - lam^2))
    data.frame(
      Construct   = cn,
      N_items     = length(lam),
      Min_loading = round(min(lam), 3),
      Avg_loading = round(mean(lam), 3),
      AVE         = round(AVE, 3),
      CR          = round(CR, 3),
      AVE_verdict = ifelse(AVE >= .50, "✓", "⚠ < .50"),
      CR_verdict  = ifelse(CR  >= .70, "✓", "⚠ < .70")
    )
  }))
}

ave_cr_v3 <- compute_AVE_CR(cfa_fit, CONSTRUCTS)
cat("\n--- Convergent Validity (AVE & CR) ---\n")
print(ave_cr_v3)

# Compute HTMT
htmt(model_definition_v0, data = cleaned_data_2024)
sink()
# EFA
# install.packages("psych")
library(psych)
library(dplyr)

# 1. Select only your items
indicator_cols <- c(
  paste0("CTA", c(1:5)), paste0("KBB", 1:5), paste0("SPI", 1:5),
  paste0("RBB", 1:5), paste0("KAB", 1:5), paste0("NIAT", 1:10)
)
efa_data <- cleaned_data_2024 %>% select(all_of(indicator_cols))

# 2. Run EFA specifying 6 factors (or let it determine freely)
# We use 'oblimin' rotation because your factors are expected to correlate
efa_result <- fa(efa_data, nfactors = 6, rotate = "oblimin", fm = "pa")
sink("BN_SEM_141026/EFA_model_v0.txt")
# 3. Print loadings, suppressing weak ones (< 0.3) for scannability
print(efa_result$loadings, cut = 0.3)
sink()

#-----------
library(lavaan)
library(semTools)

# Ensure this is the EXACT same model string used for the 2025 wave
unified_tracking_model <- "
  # 1. Unified Brand Evaluation (CTA + KBB + SPI items)
  CIVIC_ID =~ CTA1 + CTA4 + CTA5 + 
               KBB1 + KBB2 + KBB3 + KBB4 + KBB5 +
               SPI1 + SPI2 + SPI3 + SPI4 + SPI5
               
  # 2. Independent Initial Capability(KAB items)
  CAPABILITY =~ KAB3 + KAB4 + KAB5

  # 3. Willingness to Sacrifice (RBB items)
  SACRIFICE  =~ RBB2 + RBB3 + RBB5
               
  # 4. Target Intentions (All NIAT items)
  INTENTION   =~ NIAT1 + NIAT2 + NIAT3 + NIAT4 + NIAT5 + 
               NIAT6 + NIAT7 + NIAT8 + NIAT9 + NIAT10
"
UNIFIED_ITEMS <- list(
  CIVIC_ID = c("CTA1","CTA4","CTA5","KBB1","KBB2","KBB3","KBB4","SPI1","SPI2","SPI3","SPI5"),
  SACRIFICE = paste0("RBB", c(2,3,5)),
  CAPABILITY = paste0("KAB",  3:5),
  INTENTION = paste0("NIAT", 1:10)
)

setwd("D:/2026/disertasi/paper1/BN_SEM_141026")
# 1. Calculate the HTMT matrix for the 2024 Wave
htmt_2024 <- htmt(unified_tracking_model, data = cleaned_data_2024)
htmt_2025 <- htmt(unified_tracking_model, data = cleaned_data_2025)
sink("htmt_unified_tracking_model_2024_2025.txt")
print("--- 2024 HTMT Results Revised Model---")
print(htmt_2024)
print(htmt_2025)
sink()

cfa_fit24 <- cfa(unified_tracking_model, data = cleaned_data_2024, estimator = "MLR", std.lv    = TRUE)
cfa_fit25 <- cfa(unified_tracking_model, data = cleaned_data_2025, estimator = "MLR", std.lv    = TRUE)

fit_measures24 <- fitMeasures(cfa_fit24, 
                            c("chisq", "df", "pvalue", "cfi", "tli", 
                              "rmsea", "srmr"))

fit_measures25 <- fitMeasures(cfa_fit25, 
                              c("chisq", "df", "pvalue", "cfi", "tli", 
                                "rmsea", "srmr"))

sink("fit_unified_tracking_model_full.txt")
cat("FIT MEASURE 2024")
cat("\n")
print(fit_measures24)
# ── Reliability (Cronbach alpha) ──────────────────────────────────────────
cat("--- Reliability Estimates ---\n")
reliability_df <- do.call(rbind, lapply(names(UNIFIED_ITEMS), function(cn) {
  items_present <- intersect(UNIFIED_ITEMS[[cn]], names(cleaned_data_2024))
  if (length(items_present) < 2) return(NULL)
  d_sub <- data_2024_raw[, items_present]
  d_sub <- d_sub[complete.cases(d_sub), ]
  alpha_val <- tryCatch(psych::alpha(d_sub)$total$raw_alpha, error = function(e) NA)
  data.frame(Construct = cn,
             n_items   = length(items_present),
             Alpha     = round(alpha_val, 3),
             Threshold = "≥ .70")
}))
print(reliability_df)
# ── AVE dan CR ────────────────────────────────────────────────
CONSTRUCTS  <- names(UNIFIED_ITEMS)
compute_AVE_CR <- function(fit, constructs) {
  params <- parameterEstimates(fit, standardized = TRUE) %>%
    filter(op == "=~")
  do.call(rbind, lapply(constructs, function(cn) {
    lam <- params[params$lhs == cn, "std.all"]
    AVE <- mean(lam^2)
    CR  <- sum(lam)^2 / (sum(lam)^2 + sum(1 - lam^2))
    data.frame(
      Construct   = cn,
      N_items     = length(lam),
      Min_loading = round(min(lam), 3),
      Avg_loading = round(mean(lam), 3),
      AVE         = round(AVE, 3),
      CR          = round(CR, 3),
      AVE_verdict = ifelse(AVE >= .50, "✓", "⚠ < .50"),
      CR_verdict  = ifelse(CR  >= .70, "✓", "⚠ < .70")
    )
  }))
}

ave_cr_2024 <- compute_AVE_CR(cfa_fit24, CONSTRUCTS)
cat("\n--- Convergent Validity (AVE & CR) ---\n")
print(ave_cr_2024)
cat("\n--- Discriminant Validity (HTMT) ---\n")
print(htmt_2024)
cat("FIT MEASURE 2025")
cat("\n")
print(fit_measures25)
# ── A6. Reliability (Cronbach alpha) ──────────────────────────────────────────
cat("--- Reliability Estimates ---\n")
reliability_df <- do.call(rbind, lapply(names(UNIFIED_ITEMS), function(cn) {
  items_present <- intersect(UNIFIED_ITEMS[[cn]], names(cleaned_data_2024))
  if (length(items_present) < 2) return(NULL)
  d_sub <- data_2024_raw[, items_present]
  d_sub <- d_sub[complete.cases(d_sub), ]
  alpha_val <- tryCatch(psych::alpha(d_sub)$total$raw_alpha, error = function(e) NA)
  data.frame(Construct = cn,
             n_items   = length(items_present),
             Alpha     = round(alpha_val, 3),
             Threshold = "≥ .70")
}))
print(reliability_df)
# ── AVE dan CR ────────────────────────────────────────────────
CONSTRUCTS  <- names(UNIFIED_ITEMS)
compute_AVE_CR <- function(fit, constructs) {
  params <- parameterEstimates(fit, standardized = TRUE) %>%
    filter(op == "=~")
  do.call(rbind, lapply(constructs, function(cn) {
    lam <- params[params$lhs == cn, "std.all"]
    AVE <- mean(lam^2)
    CR  <- sum(lam)^2 / (sum(lam)^2 + sum(1 - lam^2))
    data.frame(
      Construct   = cn,
      N_items     = length(lam),
      Min_loading = round(min(lam), 3),
      Avg_loading = round(mean(lam), 3),
      AVE         = round(AVE, 3),
      CR          = round(CR, 3),
      AVE_verdict = ifelse(AVE >= .50, "✓", "⚠ < .50"),
      CR_verdict  = ifelse(CR  >= .70, "✓", "⚠ < .70")
    )
  }))
}

ave_cr_2025 <- compute_AVE_CR(cfa_fit25, CONSTRUCTS)
cat("\n--- Convergent Validity (AVE & CR) ---\n")
print(ave_cr_2025)
cat("\n--- Discriminant Validity (HTMT) ---\n")
print(htmt_2025)
sink()

# 1. Add a wave identifier variable to each clean dataset
cleaned_data_2024 <- cleaned_data_2024 %>% mutate(wave = "2024")
cleaned_data_2025 <- cleaned_data_2025 %>% mutate(wave = "2025")

write.csv(cleaned_data_2024, "cleaned_data_2024.csv", row.names = FALSE)
write.csv(cleaned_data_2025, "cleaned_data_2025.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
cleaned_data_2024 <- read.csv("cleaned_data_2024.csv")
cleaned_data_2025 <- read.csv("cleaned_data_2025.csv")
# 2. Combine both waves into a single longitudinal tracking dataset
combined_tracker_data <- bind_rows(cleaned_data_2024, cleaned_data_2025)

# 4. Run multi-group CFA to test for Metric Invariance
# (This ensures factor loadings are equivalent across 2024 and 2025)
metric_fit <- cfa(
  model = unified_tracking_model, 
  data = combined_tracker_data, 
  group = "wave", 
  group.equal = "loadings"
)

fit_measures <- fitMeasures(metric_fit, 
                              c("chisq", "df", "pvalue", "cfi", "tli", 
                                "rmsea", "srmr"))

sink("metric_invariance_2024_2025.txt")
# 5. Check if your model fit holds up under tracking constraints
print(fit_measures)
summary(metric_fit, fit.measures = TRUE, standardized = TRUE)
sink()

# 1. Fit the Scalar Invariance Model (Fixes both loadings and intercepts across waves)
scalar_fit <- cfa(
  model = unified_tracking_model, 
  data = combined_tracker_data, 
  group = "wave", 
  group.equal = c("loadings", "intercepts")
)
sink("scalar_fit_latent_mean_diff_24_25.txt")
# 2. Extract the true latent mean differences 
# (Group 2024 will automatically be set as the baseline 0.00 reference point)
summary(scalar_fit, fit.measures = TRUE, standardized = TRUE)
sink()

# Extract Factor Scores
#Step 1: Extract Factor Scores and Prepare the Datasets
library(lavaan)
library(dplyr)

fit_combine_track <- sem(
  model = unified_tracking_model, 
  data = combined_tracker_data, 
  group = "wave", 
  group.equal = c("loadings"), # Memaksakan factor loading sama di setiap wave
  estimator = "MLR"
)


# 3. Print out the calculated parameters
sink("fit_combine_track.txt")
summary(fit_combine_track, standardized = TRUE)
sink()



# 1. Predict factor scores (this outputs a LIST of 4 dataframes)
predicted_list <- lavPredict(fit_combine_track)

# 2. Correctly stack the groups on top of each other using bind_rows
# (lavaan outputs them in the exact order of your groups)
latent_scores <- bind_rows(as.data.frame(predicted_list[[1]]),
                           as.data.frame(predicted_list[[2]]))

# 3. Double-check row counts match up perfectly
# Your combined_tracker_data has 11099 + 2147 + 1450 + 16464 = 31,160 rows.
# Let's ensure latent_scores has exactly 31,160 rows as well:
cat("Rows in survey data:", nrow(combined_tracker_data), "\n")
cat("Rows in factor scores:", nrow(latent_scores), "\n")

# 4. Safely bind them together now that row counts match perfectly
base_variables <- combined_tracker_data %>% 
  select(wave, province_id, provinsi, k_lingkungan, k_umur, k_gender, k_pendidikan)

full_bn_dataset <- bind_cols(base_variables, latent_scores) %>%
  mutate(
    k_gender       = as.factor(k_gender),        # Categorical
    wave      = as.factor(wave),       # Categorical
    k_umur = as.factor(k_umur)   # Added: Converted to factor (Categorical)
  )

write.csv(full_bn_dataset,"full_bn_dataset.csv", row.names = FALSE)
