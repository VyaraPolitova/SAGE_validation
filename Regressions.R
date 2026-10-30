library(tidyverse)    
library(janitor)      
library(gtsummary)    
library(gt)           
library(psych)  
library(moments)
library(Hmisc)
library(lavaan)
library(semTools)
library(semPlot)
library(patchwork)
library(corrplot)
library(car)
library(mvnormtest)
library(lmtest)
library(sandwich)
library(biotools)
library(boot)
library(AER)
library(MVN)
library(conflicted)
library(MASS)
library(effectsize)

data <- read.csv("data/SAGE_MMSE_MoCA-data.csv")
head(data, 10)
names(data)

data <- data %>%
  mutate(Gender = dplyr::recode(Gender, "ж" = "Female", "м" = "Male")) %>%
  rename(MMSE_total = 'MMSE.total',
         MoCA_total = 'MoCA.total',
         SAGE_total = 'SAGEtotal') #%>%


data_healthy <- data %>%
  dplyr::filter(Status == "H")
data_clinical <- data %>%
  dplyr::filter(Status == "C")

colSums(is.na(data))
colSums(is.na(data_healthy))

data$Age_c <- data$Age - mean(data$Age, na.rm = TRUE)

# --- Primary OLS models ---
model_mmse <- lm(MMSE_total ~ Age_c + Years_of_education + Status, data = data)
model_moca <- lm(MoCA_total ~ Age_c + Years_of_education + Status, data = data)
model_sage <- lm(SAGE_total ~ Age_c + Years_of_education + Status, data = data)

# --- Primary OLS models age squared ---
model_mmse_quad <- lm(MMSE_total ~ Age_c + I(Age_c^2) + Years_of_education + Status, data = data)
model_moca_quad <- lm(MoCA_total ~ Age_c + I(Age_c^2) + Years_of_education + Status, data = data)
model_sage_quad <- lm(SAGE_total ~ Age_c + I(Age_c^2) + Years_of_education + Status, data = data)

#Exploratory
#model <- lm(cbind(MMSE_total, MoCA_total, SAGE_total) ~ Age + Years_of_education + Status, data = data)
#summary(model)


# --- Core Assumptions checks ---

# --- Independence of errors --- generally not needed
# Formal Durbin-Watson Test (Checks for lag-1 correlation)
lmtest::dwtest(model_mmse_quad)
lmtest::dwtest(model_moca_quad)
lmtest::dwtest(model_sage_quad)

# --- Linearity --- FIXED WITH AGE^2
#Residuals vs. Fitted plot
par(mfrow = c(2,2), mar = c(2,2,2,2))
plot(model_mmse_quad)
plot(model_moca_quad)
plot(model_sage_quad) 
#per predictor
car::crPlots(model_mmse_quad) 
car::crPlots(model_moca_quad)
car::crPlots(model_sage_quad)
#formal linearity check - Ramsey p > 0.05
lmtest::resettest(model_mmse_quad, power = 2:3, type = "fitted")
lmtest::resettest(model_moca_quad, power = 2:3, type = "fitted")
lmtest::resettest(model_sage_quad, power = 2:3, type = "fitted")

# --- Multicollinearity --- IT IS OK
# Years_of_education VIF should be < 5 (ideally < 2.5)
car::vif(model_mmse_quad)
car::vif(model_moca_quad)
car::vif(model_sage_quad)

# --- Homoscedasticity --- NEEDS FIXING
# p > 0.05 → fail to reject H₀ 
lmtest::bptest(model_mmse_quad) 
lmtest::bptest(model_moca_quad) 
lmtest::bptest(model_sage_quad) 

# --- Normality of residuals --- 
# p > 0.05
shapiro.test(residuals(model_mmse_quad)) 
shapiro.test(residuals(model_moca_quad)) 
shapiro.test(residuals(model_sage_quad)) 

par(mfrow = c(1, 3))
car::qqPlot(model_mmse_quad)
car::qqPlot(model_moca_quad)
car::qqPlot(model_sage_quad)

# --- Influence diagnostics --- Outliers and Influential Observations
# Cook's D < 0.5
influencePlot(model_mmse_quad)
influencePlot(model_moca_quad)
influencePlot(model_sage_quad)

data[c(21, 26, 75, 88), c("MMSE_total", "MoCA_total", "SAGE_total", "Age", "Years_of_education", "Status")]

# --- Robust regression (MM-estimation) ---
model_mmse_rlm <- rlm(MMSE_total ~ Age_c + I(Age_c^2) + Years_of_education + Status, data = data, method = "MM")
model_moca_rlm <- rlm(MoCA_total ~ Age_c + I(Age_c^2) + Years_of_education + Status, data = data, method = "MM")
model_sage_rlm <- rlm(SAGE_total ~ Age_c + I(Age_c^2) + Years_of_education + Status, data = data, method = "MM")

get_rlm_p <- function(model) {
  cf <- summary(model)$coefficients
  cbind(cf, "p value" = round(2 * (1 - pnorm(abs(cf[, "t value"]))), 5))
}

get_rlm_p(model_mmse_rlm)
get_rlm_p(model_moca_rlm)
get_rlm_p(model_sage_rlm)

# --- Tobit for ceiling effects --- 
tobit_mmse <- tobit(MMSE_total ~ Age_c + I(Age_c^2) + Years_of_education + Status, right = 30, data = data, robust = TRUE)
tobit_moca <- tobit(MoCA_total ~ Age_c + I(Age_c^2) + Years_of_education + Status, right = 30, data = data, robust = TRUE)
tobit_sage <- tobit(SAGE_total ~ Age_c + I(Age_c^2) + Years_of_education + Status, right = 22, data = data, robust = TRUE)



strip_zero <- function(x) sub("^(-?)0\\.", "\\1.", x)
p_fmt <- function(p) ifelse(p < .001, "< .001", strip_zero(sprintf("%.3f", p)))

partial_f2 <- function(full_model, data, dv, drop_var) {
  full_terms <- attr(terms(full_model), "term.labels")
  reduced_terms <- setdiff(full_terms, drop_var)
  reduced_formula <- as.formula(paste(dv, "~", paste(reduced_terms, collapse = " + ")))
  reduced_model <- lm(reduced_formula, data = data)
  r2_full <- summary(full_model)$r.squared
  r2_reduced <- summary(reduced_model)$r.squared
  round((r2_full - r2_reduced) / (1 - r2_full), 3)
}

partial_f2(model_mmse_quad, data, "MMSE_total", "Years_of_education")
partial_f2(model_moca_quad, data, "MoCA_total", "Years_of_education")
partial_f2(model_sage_quad, data, "SAGE_total", "Years_of_education")

#Helper function formatting
clean_predictor_names <- function(x) {
  x <- gsub("^Age_c$", "Age (mean centered)", x)
  x <- gsub("^I\\(Age_c\\^2\\)$", "Age² (mean centered)", x)
  x <- gsub("^Years_of_education$", "Years of education", x)
  x <- gsub("^\\(Intercept\\)$", "Intercept", x)
  # Handles dummy-coded levels of StatusN if they appear in the table
  x <- gsub("^StatusN.*", "Diagnostic status (Normative vs. Clinical)", x)
  return(x)
}

# --- Table generated functions ---
#OLS table function
generate_table1_ols <- function(model, data, dv) {
  hc3 <- coeftest(model, vcov = vcovHC(model, type = "HC3"))
  betas <- standardize_parameters(model, method = "basic")
  coef_names <- rownames(hc3)
  
  results <- data.frame(
    Predictor = coef_names,
    b = hc3[, "Estimate"],
    SE_HC3 = hc3[, "Std. Error"],
    t_stat = hc3[, "t value"],
    p_val = hc3[, "Pr(>|t|)"],
    CI_Lower = hc3[, "Estimate"] - 1.96 * hc3[, "Std. Error"],
    CI_Upper = hc3[, "Estimate"] + 1.96 * hc3[, "Std. Error"],
    stringsAsFactors = FALSE
  )
  
  results$Predictor <- clean_predictor_names(results$Predictor)
  
  results$Beta <- sapply(coef_names, function(x) {
    val <- betas$Std_Coefficient[betas$Parameter == x]
    if (length(val) == 0) NA_real_ else val
  })
  results$Beta[results$Predictor == "(Intercept)"] <- NA_real_
  
  results$b_CI <- sprintf("%.2f [%.2f, %.2f]", results$b, results$CI_Lower, results$CI_Upper)
  results$Beta <- sapply(results$Beta, function(x) if (is.na(x)) NA_character_ else strip_zero(sprintf("%.2f", x)))
  results$SE_HC3 <- sprintf("%.2f", results$SE_HC3)
  results$t_stat <- sprintf("%.2f", results$t_stat)
  results$p_val  <- p_fmt(results$p_val)
  
  list(Coefficients = results[, c("Predictor", "b_CI", "Beta", "SE_HC3", "t_stat", "p_val")])
}


#Robust r table function
generate_table2_rlm <- function(model) {
  s <- summary(model)$coefficients
  results <- data.frame(
    Predictor = rownames(s),
    b = s[, "Value"],
    SE = s[, "Std. Error"],
    t_stat = s[, "t value"],
    p_val = 2 * (1 - pnorm(abs(s[, "t value"]))),
    stringsAsFactors = FALSE
  )
  
  results$Predictor <- clean_predictor_names(results$Predictor)
  
  results$b <- sprintf("%.2f", results$b)
  results$SE <- sprintf("%.2f", results$SE)
  results$t_stat <- sprintf("%.2f", results$t_stat)
  results$p_val <- p_fmt(results$p_val)
  results
}


#Tobit table function
generate_table3_tobit <- function(model) {
  s <- summary(model)$coefficients
  results <- data.frame(
    Predictor = rownames(s),
    b = s[, "Estimate"],
    SE = s[, "Std. Error"],
    z_stat = s[, "z value"],
    p_val = s[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  results <- results[results$Predictor != "Log(scale)", ]
  results$Predictor <- clean_predictor_names(results$Predictor)
  
  results$b <- sprintf("%.2f", results$b)
  results$SE <- sprintf("%.2f", results$SE)
  results$z_stat <- sprintf("%.2f", results$z_stat)
  results$p_val <- p_fmt(results$p_val)
  results
}


#Summary table education coeffs 
extract_education_row <- function(ols_mod, rlm_mod, tob_mod) {
  hc3 <- coeftest(ols_mod, vcov = vcovHC(ols_mod, type = "HC3"))
  b_ols <- sprintf("%.2f", hc3["Years_of_education", "Estimate"])
  se_ols <- sprintf("%.2f", hc3["Years_of_education", "Std. Error"])
  p_raw <- hc3["Years_of_education", "Pr(>|t|)"]
  p_txt <- ifelse(p_raw < .001, "< .001", paste0("= ", strip_zero(sprintf("%.3f", p_raw))))
  
  beta_results <- standardize_parameters(ols_mod, method = "basic")
  beta_val <- beta_results[beta_results$Parameter == "Years_of_education", "Std_Coefficient"]
  beta <- strip_zero(sprintf("%.2f", beta_val))
  
  rlm_s <- summary(rlm_mod)$coefficients
  b_rlm <- sprintf("%.2f", rlm_s["Years_of_education", "Value"])
  se_rlm <- sprintf("%.2f", rlm_s["Years_of_education", "Std. Error"])
  
  tob_s <- summary(tob_mod)$coefficients
  b_tob <- sprintf("%.2f", tob_s["Years_of_education", "Estimate"])
  se_tob <- sprintf("%.2f", tob_s["Years_of_education", "Std. Error"])
  
  data.frame(
    OLS_HC3 = sprintf("%s (%s), *p*  %s", b_ols, se_ols, p_txt),
    Std_Beta = beta,
    RLM_MM = sprintf("%s (%s)", b_rlm, se_rlm),
    Tobit = sprintf("%s (%s)", b_tob, se_tob),
    stringsAsFactors = FALSE
  )
}

#Model fit effect size summary, with Holm
build_model_fit_summary <- function(m_mmse, m_moca, m_sage, dat) {
  raw_p <- c(
    coeftest(m_mmse, vcov = vcovHC(m_mmse, type = "HC3"))["Years_of_education", "Pr(>|t|)"],
    coeftest(m_moca, vcov = vcovHC(m_moca, type = "HC3"))["Years_of_education", "Pr(>|t|)"],
    coeftest(m_sage, vcov = vcovHC(m_sage, type = "HC3"))["Years_of_education", "Pr(>|t|)"]
  )
  p_holm <- p.adjust(raw_p, method = "holm")
  
  data.frame(
    Instrument = c("MMSE", "MoCA", "SAGE"),
    R2 = c(
      strip_zero(sprintf("%.2f", summary(m_mmse)$r.squared)),
      strip_zero(sprintf("%.2f", summary(m_moca)$r.squared)),
      strip_zero(sprintf("%.2f", summary(m_sage)$r.squared))
    ),
    Adj_R2 = c(
      strip_zero(sprintf("%.2f", summary(m_mmse)$adj.r.squared)),
      strip_zero(sprintf("%.2f", summary(m_moca)$adj.r.squared)),
      strip_zero(sprintf("%.2f", summary(m_sage)$adj.r.squared))
    ),
    Education_f2 = c(
      partial_f2(m_mmse, dat, "MMSE_total", "Years_of_education"),
      partial_f2(m_moca, dat, "MoCA_total", "Years_of_education"),
      partial_f2(m_sage, dat, "SAGE_total", "Years_of_education")
    ),
    p_holm = p_fmt(p_holm),
    stringsAsFactors = FALSE
  )
}


# --- Tables function ---
theme_apa_gt <- function(gt_tbl, source_note_text = NULL) {
  tbl <- gt_tbl %>%
    tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) %>%
    tab_style(style = cell_borders(sides = "top", color = "black", weight = px(1)),
              locations = cells_source_notes()) %>%
    tab_options(
      table.border.top.style = "solid", table.border.top.width = px(2), table.border.top.color = "black",
      table.border.bottom.style = "solid", table.border.bottom.width = px(2), table.border.bottom.color = "black",
      column_labels.border.bottom.style = "solid", column_labels.border.bottom.width = px(1), column_labels.border.bottom.color = "black",
      column_labels.border.top.style = "none",
      table_body.border.bottom.style = "none", table_body.hlines.style = "none", table_body.vlines.style = "none",
      stub.border.style = "none", heading.border.bottom.style = "none",
      data_row.padding = px(4), table.font.names = "Times New Roman", table.font.size = px(12), table.width = pct(80)
    )
  if (!is.null(source_note_text)) tbl <- tbl %>% tab_source_note(source_note = md(source_note_text))
  tbl
}

build_table1_gt <- function(m_mmse, m_moca, m_sage, dat) {
  df_mmse <- cbind(Screener = "MMSE", generate_table1_ols(m_mmse, dat, "MMSE_total")$Coefficients)
  df_moca <- cbind(Screener = "MoCA", generate_table1_ols(m_moca, dat, "MoCA_total")$Coefficients)
  df_sage <- cbind(Screener = "SAGE", generate_table1_ols(m_sage, dat, "SAGE_total")$Coefficients)
  rbind(df_mmse, df_moca, df_sage) %>%
    gt(groupname_col = "Screener") %>%
    cols_label(Predictor = "Predictor", b_CI = md("*b* [95% CI]"), Beta = md("β"),
               SE_HC3 = md("*SE*<sub>HC3</sub>"), t_stat = md("*t*"), p_val = md("*p*")) %>%
    theme_apa_gt(source_note_text = "*Note.* Unstandardized coefficients (*b*) with 95% CIs, robust standard errors (*SE*<sub>HC3</sub>), and standardized coefficients (β) are presented. Leading zeros omitted for bounded parameters (β, *p*).")
}

build_table2_gt <- function(m_mmse, m_moca, m_sage) {
  df_mmse <- cbind(Screener = "MMSE", generate_table2_rlm(m_mmse))
  df_moca <- cbind(Screener = "MoCA", generate_table2_rlm(m_moca))
  df_sage <- cbind(Screener = "SAGE", generate_table2_rlm(m_sage))
  rbind(df_mmse, df_moca, df_sage) %>%
    gt(groupname_col = "Screener") %>%
    cols_label(Predictor = "Predictor", b = md("*b*"), SE = md("*SE*"), t_stat = md("*t*"), p_val = md("*p*")) %>%
    theme_apa_gt(source_note_text = "*Note.* Robust MM-estimation (`rlm`) models with bisquare weighting. Two-tailed *p*-values via asymptotic normal approximation.")
}

build_table3_gt <- function(m_mmse, m_moca, m_sage) {
  df_mmse <- cbind(Screener = "MMSE", generate_table3_tobit(m_mmse))
  df_moca <- cbind(Screener = "MoCA", generate_table3_tobit(m_moca))
  df_sage <- cbind(Screener = "SAGE", generate_table3_tobit(m_sage))
  rbind(df_mmse, df_moca, df_sage) %>%
    gt(groupname_col = "Screener") %>%
    cols_label(Predictor = "Predictor", b = md("*b*"), SE = md("*SE*"), z_stat = md("*z*"), p_val = md("*p*")) %>%
    theme_apa_gt(source_note_text = "*Note.* Tobit models account for right-censoring ceiling limits (MMSE = 30, MoCA = 30, SAGE = 22).")
}

build_table4_gt <- function(m_mmse_ols, m_mmse_rlm, m_mmse_tob,
                            m_moca_ols, m_moca_rlm, m_moca_tob,
                            m_sage_ols, m_sage_rlm, m_sage_tob) {
  df_comp <- rbind(
    MMSE = extract_education_row(m_mmse_ols, m_mmse_rlm, m_mmse_tob),
    MoCA = extract_education_row(m_moca_ols, m_moca_rlm, m_moca_tob),
    SAGE = extract_education_row(m_sage_ols, m_sage_rlm, m_sage_tob)
  )
  df_comp <- cbind(Instrument = rownames(df_comp), df_comp)
  df_comp %>%
    gt() %>%
    cols_label(Instrument = "Cognitive Screener", OLS_HC3 = md("OLS (*SE*<sub>HC3</sub>), *p*"),
               Std_Beta = md("Std. β"), RLM_MM = md("Robust MM (*SE*)"), Tobit = md("Tobit (*SE*)")) %>%
    fmt_markdown(columns = OLS_HC3) %>%
    theme_apa_gt(source_note_text = "*Note.* Comparison of education coefficients across specification models. Standardized β calculated from primary OLS models. *p*-values shown are unadjusted; see Table 5 for Holm-adjusted *p*-values.")
}

build_table5_gt <- function(m_mmse, m_moca, m_sage, dat) {
  build_model_fit_summary(m_mmse, m_moca, m_sage, dat) %>%
    gt() %>%
    cols_label(Instrument = "Cognitive Screener", R2 = md("*R*²"), Adj_R2 = md("Adj. *R*²"),
               Education_f2 = md("Education *f*²"), p_holm = md("*p* (Holm-adjusted)")) %>%
    theme_apa_gt(source_note_text = "*Note.* Education *f*² represents the unique variance explained by education, controlling for age and diagnostic group. *p*-values Holm-adjusted for multiple comparisons across the three instruments.")
}


gt_table1 <- build_table1_gt(model_mmse_quad, model_moca_quad, model_sage_quad, data)
gt_table2 <- build_table2_gt(model_mmse_rlm, model_moca_rlm, model_sage_rlm)
gt_table3 <- build_table3_gt(tobit_mmse, tobit_moca, tobit_sage)
gt_table4 <- build_table4_gt(
  model_mmse_quad, model_mmse_rlm, tobit_mmse,
  model_moca_quad, model_moca_rlm, tobit_moca,
  model_sage_quad, model_sage_rlm, tobit_sage
)
gt_table5 <- build_table5_gt(model_mmse_quad, model_moca_quad, model_sage_quad, data)

gt_table1
gt_table2
gt_table3
gt_table4
gt_table5

#F-stat and dfs
summary(model_mmse_quad)$fstatistic
summary(model_moca_quad)$fstatistic
summary(model_sage_quad)$fstatistic



# --- MANOVA --- 

# --- Correlation among outcomes --- good - reference the convergent validity in the text
cor(data %>% dplyr::select(MMSE_total, MoCA_total, SAGE_total), use = "complete.obs")

model <- lm(cbind(MMSE_total, MoCA_total, SAGE_total) ~ Age_c + I(Age_c^2) + Years_of_education + Status, data = data)

#using Pillai due to normality violations
mv_manova <- manova(cbind(MMSE_total, MoCA_total, SAGE_total) ~ Age_c + I(Age_c^2) + Years_of_education + Status, data = data)
summary(mv_manova, test = "Pillai")

# Sensitivity
data_sens <- data[-c(21, 26, 88), ]

model_mmse_sens <- lm(MMSE_total ~ Age_c + I(Age_c^2) + Years_of_education + Status, data = data_sens)
model_moca_sens <- lm(MoCA_total ~ Age_c + I(Age_c^2) + Years_of_education + Status, data = data_sens)
model_sage_sens <- lm(SAGE_total ~ Age_c + I(Age_c^2) + Years_of_education + Status, data = data_sens)

mv_manova_sens <- manova(cbind(MMSE_total, MoCA_total, SAGE_total) ~ Age_c + I(Age_c^2) + Years_of_education + Status, data = data_sens)
summary(mv_manova_sens, test = "Pillai")

# Multivariate lm cross-check (should match MANOVA)
car::Anova(model, test.statistic = "Pillai")

#multivariate normality
resid_matrix <- cbind(residuals(model_mmse_quad), residuals(model_moca_quad), residuals(model_sage))
MVN::mvn(data = resid_matrix, mvn = "mardia")   

resid_matrix_sens <- cbind(residuals(model_mmse_sens), residuals(model_moca_sens), residuals(model_sage_sens))
MVN::mvn(data = resid_matrix_sens, mvn = "mardia")

# --- Multivariate outliers (Mahalanobis distance) ---
D2 <- mahalanobis(resid_matrix, colMeans(resid_matrix), cov(resid_matrix))
sort(pchisq(D2, df = 3, lower.tail = FALSE))[1:10] # smallest p-values = most extreme cases

# --- Homogeneity of covariance matrices across Status groups ---
# it is overly sensitive
boxM(cbind(data$MMSE_total, data$MoCA_total, data$SAGE_total), data$Status)

