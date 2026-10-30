library(tidyverse)    
library(janitor)      
library(gtsummary)    
library(gt)           
library(psych)  
library(moments)
library(e1071)
library(conflicted)


data <- read.csv("data/SAGE_MMSE_MoCA-data.csv")
head(data, 10)

data <- data %>%
  mutate(Gender = dplyr::recode(Gender, "ж" = "Female", "м" = "Male")) %>%
  rename(MMSE_total = 'MMSE-total',
         MoCA_total = 'MoCA-total',
         SAGE_total = 'SAGEtotal')

data_healthy <- data %>%
  dplyr::filter(Status == "N")
data_clinical <- data %>%
  dplyr::filter(Status == "C")

colSums(is.na(data))
colSums(is.na(data_healthy))

# N
print(nrow(data))
number_healthy <- nrow(data_healthy)
number_clinical <- nrow(data_clinical)
print(number_healthy)
print(number_clinical)


#Descriptive stat for sample (demographics)
tbl_descriptives_sample <- data %>%
  dplyr::mutate(Status = dplyr::recode(Status, "N" = "Normative", "C" = "Clinical"),
                Status = factor(Status, levels = c("Normative", "Clinical"))) %>%
  dplyr::select(Status, Gender, Age, Years_of_education) %>%
  tbl_summary(
    by = Status,
    type = list(
      Age ~ "continuous2",
      Years_of_education ~ "continuous2"
    ),
    statistic = list(
      Age ~ c("{mean} ({sd})", "{median} [{p25}, {p75}]", "{min} - {max}"),
      Years_of_education ~ c("{mean} ({sd})", "{median} [{p25}, {p75}]", "{min} - {max}"),
      Gender ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1),
    label = list(
      Gender ~ "Gender",
      Age ~ "Age",
      Years_of_education ~ "Years of Education"
    )
  ) %>%
  add_overall(col_label = "**Total**") %>%
  add_p(pvalue_fun = ~ gtsummary::style_pvalue(.x, digits = 3)) %>%
  modify_table_body(
    ~ .x %>%
      dplyr::mutate(
        label = dplyr::case_when(
          label == "Mean (SD)" ~ "M (SD)",
          label == "Median [Q1, Q3]" ~ "Mdn [IQR]",
          label == "Min - Max" ~ "Range",
          TRUE ~ label
        )
      )
  ) %>%
  modify_header(
    label = "**Variable**",
    p.value = gt::html("<i>p</i>-value"),
    stat_0 = gt::html("<b>Total</b><br>N = {N}"),
    stat_1 = gt::html("<b>Normative</b><br>n = {n}"),
    stat_2 = gt::html("<b>Clinical</b><br>n = {n}")
  ) %>%
  bold_labels() %>%
  modify_footnote(
    all_stat_cols() ~ NA,
    p.value ~ NA
  ) %>%
  as_gt() %>%
  gt::text_transform(
    locations = gt::cells_body(
      columns = p.value,
      rows = variable == "Gender" & row_type == "label"
    ),
    fn = function(x) paste0(x, "<sup>a</sup>")
  ) %>%
  gt::text_transform(
    locations = gt::cells_body(
      columns = p.value,
      rows = variable %in% c("Age", "Years_of_education") &
        row_type == "label"
    ),
    fn = function(x) paste0(x, "<sup>b</sup>")
  ) %>%
  gt::tab_style(
    style = gt::cell_text(style = "italic"),
    locations = gt::cells_body(
      columns = label,
      rows = label %in% c("M (SD)", "Mdn [IQR]")
    )
  ) %>%
  gt::tab_style(
    style = gt::cell_borders(
      sides = "top",
      color = "black",
      weight = gt::px(1)
    ),
    locations = gt::cells_source_notes()
  ) %>%
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
    
    #table_body.border.top.style = "solid",   
    #table_body.border.top.width = gt::px(2),  
    #table_body.border.top.color = "black",  
    table_body.border.bottom.style = "none",
    
    table_body.hlines.style = "none",
    table_body.vlines.style = "none",
    
    heading.border.bottom.style = "none",
    
    row_group.border.top.style = "none",
    row_group.border.bottom.style = "none",
    
    data_row.padding = gt::px(4),
    
    table.font.names = "Times New Roman",
    table.font.size = gt::px(12),
    table.width = gt::pct(50)
  ) %>%
  gt::tab_source_note(
    source_note = gt::md(
      "*Note.* Data are presented as *n* (%) for categorical variables and *M* (*SD*), *Mdn* *[IQR]*, and range for continuous variables. *p*-values represent comparisons between the Normative and Clinical groups. <br>
    <sup>a</sup> Pearson's chi-square test.<br>
    <sup>b</sup> Wilcoxon rank-sum test."
    )
  )

tbl_descriptives_sample


#Descriptive stats for tests
# Stats overall
stats_total <- data %>%
  dplyr::select(MMSE_total, MoCA_total, SAGE_total) %>%
  pivot_longer(
    everything(),
    names_to = "Instrument",
    values_to = "Score"
  ) %>%
  mutate(
    Instrument = dplyr::recode(
      Instrument,
      "MMSE_total" = "MMSE",
      "MoCA_total" = "MoCA",
      "SAGE_total" = "SAGE"
    ),
    Status = "Total"
  ) %>%
  group_by(Instrument, Status) %>%
  summarise(
    `Mean (SD)` = sprintf(
      "%.1f (%.1f)",
      mean(Score, na.rm = TRUE),
      sd(Score, na.rm = TRUE)
    ),
    `Median [IQR]` = sprintf(
      "%.1f [%.1f, %.1f]",
      median(Score, na.rm = TRUE),
      quantile(Score, .25, na.rm = TRUE),
      quantile(Score, .75, na.rm = TRUE)
    ),
    `Range (Min-Max)` = sprintf(
      "%.1f - %.1f",
      min(Score, na.rm = TRUE),
      max(Score, na.rm = TRUE)
    ),
    `Skewness` = sprintf(
      "%.2f",
      e1071::skewness(Score, na.rm = TRUE, type = 2)
    ),
    `Excess Kurtosis` = sprintf(
      "%.2f",
      e1071::kurtosis(Score, na.rm = TRUE, type = 2)  # already excess kurtosis — no -3 needed
    ),
    .groups = "drop"
  )

# Stats by group (Clinical vs Healthy)
stats_group <- data %>%
  dplyr::select(Status, MMSE_total, MoCA_total, SAGE_total) %>%
  pivot_longer(
    -Status,
    names_to = "Instrument",
    values_to = "Score"
  ) %>%
  mutate(
    Instrument = dplyr::recode(
      Instrument,
      "MMSE_total" = "MMSE",
      "MoCA_total" = "MoCA",
      "SAGE_total" = "SAGE"
    )
  ) %>%
  group_by(Instrument, Status) %>%
  summarise(
    `Mean (SD)` = sprintf(
      "%.1f (%.1f)",
      mean(Score, na.rm = TRUE),
      sd(Score, na.rm = TRUE)
    ),
    `Median [IQR]` = sprintf(
      "%.1f [%.1f, %.1f]",
      median(Score, na.rm = TRUE),
      quantile(Score, .25, na.rm = TRUE),
      quantile(Score, .75, na.rm = TRUE)
    ),
    `Range (Min-Max)` = sprintf(
      "%.1f - %.1f",
      min(Score, na.rm = TRUE),
      max(Score, na.rm = TRUE)
    ),
    Skewness = sprintf(
      "%.2f",
      e1071::skewness(Score, na.rm = TRUE, type = 2)
    ),
    `Excess Kurtosis` = sprintf(
      "%.2f",
      e1071::kurtosis(Score, na.rm = TRUE, type = 2)  # already excess kurtosis — no -3 needed
    ),
    .groups = "drop"
  )

# Combine, pivot to long (one row per statistic), then wide by Status
tools_final <- bind_rows(stats_total, stats_group) %>%
  pivot_longer(
    cols = c(
      `Mean (SD)`,
      `Median [IQR]`,
      `Range (Min-Max)`,
      `Skewness`,
      `Excess Kurtosis`
    ),
    names_to = "Statistic",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Status,
    values_from = Value
  ) %>%
  dplyr::select(Instrument, Statistic, any_of(c("Total", "N", "C"))) %>%
  mutate(
    Statistic = dplyr::recode(
      Statistic,
      "Mean (SD)" = "*M* (*SD*)",
      "Median [IQR]" = "*Mdn* [*IQR*]",
      "Range (Min-Max)" = "Range",
      "Skewness" = "Skewness",
      "Excess Kurtosis" = "Excess kurtosis"
    ),
    Statistic = factor(
      Statistic,
      levels = c(
        "*M* (*SD*)",
        "*Mdn* [*IQR*]",
        "Range",
        "Skewness",
        "Excess kurtosis"
      )
    )
  ) %>%
  arrange(Instrument, Statistic)

# Table output
tbl_tools_final <- tools_final %>%
  gt::gt(rowname_col = "Statistic", groupname_col = "Instrument") %>%
  gt::fmt_markdown(
    columns = gt::stub()
  ) %>%
  gt::cols_label(Total = "Total", N = "Normative", C = "Clinical") %>% 
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_row_groups()
  ) %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_column_labels()
  ) %>%
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

tbl_tools_final

#Additional

#Z-scores for skewness and kurtosis
max_score <- c(MMSE = 30, MoCA = 30, SAGE = 22)

long <- data %>%
  gtsummary::select(Status, MMSE_total, MoCA_total, SAGE_total) %>%
  pivot_longer(-Status, names_to = "Instrument", values_to = "Score") %>%
  mutate(Instrument = sub("_total", "", Instrument))

diag <- long %>%
  group_by(Instrument, Status) %>%
  summarise(
    n        = n(),
    skew     = e1071::skewness(Score, type = 2),
    exkurt   = e1071::kurtosis(Score, type = 2),
    se_skew  = sqrt(6*n*(n-1)/((n-2)*(n+1)*(n+3))),
    se_kurt  = 2*se_skew*sqrt((n^2-1)/((n-3)*(n+5))),
    z_skew   = skew/se_skew,
    z_kurt   = exkurt/se_kurt,
    SW_p     = shapiro.test(Score)$p.value,
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))
diag

#for the total sample
tot_samp <- data %>%
  gtsummary::select(MMSE_total, MoCA_total, SAGE_total) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Instrument",
    values_to = "Score"
  ) %>%
  mutate(Instrument = sub("_total", "", Instrument)) %>%
  group_by(Instrument) %>%
  summarise(
    n        = sum(!is.na(Score)),
    skew     = e1071::skewness(Score, na.rm = TRUE, type = 2),
    exkurt   = e1071::kurtosis(Score, na.rm = TRUE, type = 2),
    se_skew  = sqrt(6*n*(n-1)/((n-2)*(n+1)*(n+3))),
    se_kurt  = 2*se_skew*sqrt((n^2-1)/((n-3)*(n+5))),
    z_skew   = skew/se_skew,
    z_kurt   = exkurt/se_kurt,
    SW_p     = shapiro.test(Score[!is.na(Score)])$p.value,
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

tot_samp


#Numbers at max for Ceiling effects
ceiling_check <- function(df, var, max_val) {
  df %>%
    summarise(
      max_score = max({{ var }}, na.rm = TRUE),
      at_max = sum({{ var }} == max_val, na.rm = TRUE),
      pct_at_max = sprintf("%.2f", mean({{ var }} == max_val, na.rm = TRUE) * 100),
      within_1 = sum({{ var }} >= (max_val - 1), na.rm = TRUE),
      pct_within_1 = sprintf("%.2f", mean({{ var }} >= (max_val - 1), na.rm = TRUE) * 100)
    )
}

ceiling_check(data, MMSE_total, 30)

ceiling_check(data, SAGE_total, 22)

ceiling_check(data, MoCA_total, 30)

