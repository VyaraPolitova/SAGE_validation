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
library(effectsize)
library(conflicted)
library(pwr)

data <- read.csv("data/SAGE_MMSE_MoCA-data.csv")
head(data, 10)

data <- data %>%
  mutate(Gender = case_match(Gender, "ж" ~ "Female", "м" ~ "Male")) %>%
  rename(MMSE_total = 'MMSE.total',
         MoCA_total = 'MoCA.total',
         SAGE_total = 'SAGEtotal')

data_healthy <- data %>%
  dplyr::filter(Status == "N")

data_clinical <- data %>%
  dplyr::filter(Status == "C")

colSums(is.na(data))
colSums(is.na(data_healthy))

strip_zero <- function(x) sub("^(-?)0\\.", "\\1.", x)

#Pipeline
known_groups_pipeline <- function(var_name, healthy_df, clinical_df, full_df) {
  formula_obj <- as.formula(paste(var_name, "~ Status"))
  
  #Assumption checks
  #Normality - already estimated in the descriptive stats
  sh_h <- shapiro.test(healthy_df[[var_name]])
  sh_c <- shapiro.test(clinical_df[[var_name]])
  #Homogeneity of Variance (Brown-Forsythe variant) - not needed
  lev  <- leveneTest(formula_obj, data = full_df, center = median)
  #Wilcox test/Mann-Whitney
  mwu <- wilcox.test(formula_obj, data = full_df, exact = FALSE, correct = TRUE)
  #Effect Size - in the second table
  rb_eff <- rank_biserial(formula_obj, data = full_df)
  #Post-hoc power
  d_approx <- 2 * rb_eff$r_rank_biserial / sqrt(1 - rb_eff$r_rank_biserial^2)
  pw <- pwr.t2n.test(n1 = nrow(healthy_df), n2 = nrow(clinical_df), d = d_approx, sig.level = 0.05)
  
  list(
    assumptions = tibble(
      Instrument = var_name,
      Shapiro_W_Healthy  = sprintf("%.3f", sh_h$statistic),
      Shapiro_p_Healthy  = ifelse(sh_h$p.value < .001, "< .001", sprintf("%.3f", sh_h$p.value)),
      Shapiro_W_Clinical = sprintf("%.3f", sh_c$statistic),
      Shapiro_p_Clinical = ifelse(sh_c$p.value < .001, "< .001", sprintf("%.3f", sh_c$p.value)),
      Levene_F = sprintf("%.3f", lev$`F value`[1]),
      Levene_p = ifelse(lev$`Pr(>F)`[1] < .001, "< .001", sprintf("%.3f", lev$`Pr(>F)`[1]))
    ),
    results = tibble(
      Instrument = var_name,
      W = sprintf("%.0f", mwu$statistic),
      p_value = mwu$p.value,
      rank_biserial_r = strip_zero(sprintf("%.2f", rb_eff$r_rank_biserial)),
      CI = sprintf(
        "[%s, %s]",
        strip_zero(sprintf("%.2f", rb_eff$CI_low)),
        strip_zero(sprintf("%.2f", rb_eff$CI_high))
      ),
      power = strip_zero(sprintf("%.2f", pw$power))
    )
  )
}


kg_mmse <- known_groups_pipeline("MMSE_total", data_healthy, data_clinical, data)
kg_moca <- known_groups_pipeline("MoCA_total", data_healthy, data_clinical, data)
kg_sage <- known_groups_pipeline("SAGE_total", data_healthy, data_clinical, data)

assumptions_table <- bind_rows(kg_mmse$assumptions, kg_moca$assumptions, kg_sage$assumptions)  
results_table <- bind_rows(kg_mmse$results, kg_moca$results, kg_sage$results) %>%
  mutate(
    p_holm = p.adjust(p_value, method = "holm"),
    p_show = ifelse(p_holm < .001, "< .001", sprintf("%.3f", p_holm))
  )


# --- Table 1 Mann-Whitney results + effect size + power ---
results_table %>%
  gtsummary::select(
    Instrument,
    W,
    p_show,
    rank_biserial_r,
    CI,
    power
  ) %>%
  mutate(
    Instrument = dplyr::recode(
      Instrument,
      "MMSE_total" = "MMSE",
      "MoCA_total" = "MoCA",
      "SAGE_total" = "SAGE"
    )
  ) %>%
  gt() %>%
  cols_label(
    Instrument = "Instrument",
    W = md("Mann-Whitney $U$"),
    p_show = md("$p$-value"),
    rank_biserial_r = md("Rank-biserial $r$"),
    CI = "95% CI",
    power = "Post-hoc power"
  ) %>%
  tab_source_note(source_note = md("*Note.* $p$-values were calculated using Mann-Whitney $U$ tests with continuity correction and adjusted using the Holm method for multiple comparisons across the three instruments.")) %>%
  gt::tab_style(style = gt::cell_text(weight = "bold"), locations = gt::cells_column_labels()) %>%
  gt::tab_options(
    table.border.top.style = "solid",
    table.border.top.width = gt::px(2),
    table.border.top.color = "black",
    
    table.border.bottom.style = "solid",
    table.border.bottom.width = gt::px(2),
    table.border.bottom.color = "black",
    
    column_labels.border.bottom.style = "solid",
    column_labels.border.bottom.width = gt::px(1),
    column_labels.border.bottom.color = "black",
    column_labels.border.top.style = "none",
    
    table_body.border.bottom.style = "solid",
    table_body.border.bottom.color = "black",
    table_body.border.bottom.width = gt::px(0.5),
    table_body.hlines.style = "none",
    table_body.vlines.style = "none",
    
    stub.border.style = "none",
    heading.border.bottom.style = "none",
    row_group.border.top.style = "none",
    row_group.border.bottom.style = "none",
    data_row.padding = gt::px(4),
    table.font.names = "Times New Roman",
    table.font.size = gt::px(12),
    table.width = gt::pct(50)
  )

# --- Table 2 Assumption checks --- NOT NEEDED
assumptions_table %>%
  gt() %>%
  cols_label(
    Instrument = "Instrument",
    Shapiro_W_Healthy = "W (Normative)",
    Shapiro_p_Healthy = "p (Normative)",
    Shapiro_W_Clinical = "W (Clinical)",
    Shapiro_p_Clinical = "p (Clinical)",
    Levene_F = "Levene's F",
    Levene_p = "Levene's p"
  ) %>%
  tab_spanner(label = "Shapiro-Wilk Normality",
              columns = c(Shapiro_W_Healthy, Shapiro_p_Healthy, Shapiro_W_Clinical, Shapiro_p_Clinical)) %>%
  tab_spanner(label = "Homogeneity of Variance", columns = c(Levene_F, Levene_p)) %>%
  tab_header(title = "Known-Groups Validity: Assumption Checks") %>%
  gt::tab_style(style = gt::cell_text(weight = "bold"), locations = gt::cells_column_labels()) %>%
  gt::tab_style(style = gt::cell_text(weight = "bold"), locations = gt::cells_column_spanners()) %>%
  gt::tab_options(
    table.border.top.style = "solid",
    table.border.top.width = gt::px(2),
    table.border.top.color = "black",
    table.border.bottom.style = "solid",
    table.border.bottom.width = gt::px(2),
    table.border.bottom.color = "black",
    column_labels.border.bottom.style = "solid",
    column_labels.border.bottom.width = gt::px(1),
    column_labels.border.bottom.color = "black",
    column_labels.border.top.style = "none",
    table_body.border.bottom.style = "none",
    table_body.hlines.style = "none",
    table_body.vlines.style = "none",
    stub.border.style = "none",
    heading.border.bottom.style = "none",
    row_group.border.top.style = "none",
    row_group.border.bottom.style = "none",
    data_row.padding = gt::px(4),
    table.font.names = "Times New Roman",
    table.font.size = gt::px(12),
    table.width = gt::pct(50)
  )

