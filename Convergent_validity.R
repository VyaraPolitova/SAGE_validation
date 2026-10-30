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
library(conflicted)
library(cowplot) 


data <- read.csv("data/SAGE_MMSE_MoCA-data.csv")
head(data, 10)

data <- data %>%
  mutate(Gender = dplyr::recode(Gender, "ж" = "Female", "м" = "Male")) %>%
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

#Spearman
Sp_corr <- data %>%
  gtsummary::select(MMSE_total, MoCA_total, SAGE_total)

sp_corr_result <- corr.test(Sp_corr, method = "spearman", adjust="holm", minlength = 25)


#Spearman correlations stratified by group
sp_corr_result_healthy <- corr.test(data_healthy %>% gtsummary::select(MMSE_total, MoCA_total, SAGE_total), method = "spearman", adjust="holm", minlength = 25)
sp_corr_result_clinical <- corr.test(data_clinical %>% gtsummary::select(MMSE_total, MoCA_total, SAGE_total), method = "spearman", adjust="holm", minlength = 25)

print(sp_corr_result_healthy$ci)
print(sp_corr_result_clinical$ci)

ci_df_strat <- function(corr_result, Status) {
  corr_result$ci %>%
    rownames_to_column("Pair") %>%
    mutate(
      Pair = dplyr::recode(Pair,
                    "MMSE_total-MoCA_total" = "MMSE - MoCA",
                    "MMSE_total-SAGE_total" = "MMSE - SAGE",
                    "MoCA_total-SAGE_total" = "MoCA - SAGE"),
      Group = Status
    ) %>%
    transmute(
      Group,
      Pair,
      rho = strip_zero(sprintf("%.2f", r)),
      CI = sprintf("[%s, %s]", strip_zero(sprintf("%.2f", lower)), strip_zero(sprintf("%.2f", upper))),
      p_value = ifelse(p < .001, "< .001", strip_zero(sprintf("%.3f", p)))
    )
}
  
ci_full_labeled <- ci_df_strat(sp_corr_result, "Full sample")
  

ci_healthy_labeled  <- ci_df_strat(sp_corr_result_healthy, "Normative")
ci_clinical_labeled <- ci_df_strat(sp_corr_result_clinical, "Clinical")

cor_table_by_group <- bind_rows(ci_full_labeled, ci_clinical_labeled, ci_healthy_labeled)

cor_table_by_group %>%
  gt(groupname_col = "Group") %>%
  cols_label(Pair = "Instrument Pair", rho = md("Spearman's ρ"), CI = "95% CI", p_value = md("$p$-value")) %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_column_labels()
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
    
    table_body.border.bottom.style = "none",
    table_body.hlines.style = "none",
    table_body.vlines.style = "none",
    stub.border.style = "none",
    heading.border.bottom.style = "none",
    
    data_row.padding = gt::px(4),
    table.font.names = "Times New Roman",
    table.font.size = gt::px(12),
    table.width = gt::pct(60)
    ) %>%
  gt::tab_source_note(
      source_note = gt::md(
        "*Note.* $p$-values were adjusted using the Holm method for multiple comparisons within each sample (full, normative, and clinical)."
      )
  )


#Plots
#Scatter plots together 
p1 <- ggplot(data, aes(MMSE_total, MoCA_total, color = Status)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Status)) +   
  geom_smooth(method = "lm", color = "gray40", linetype = "dashed", se = FALSE) +
  theme_minimal() +
  labs(x = "MMSE total scores", y = "MoCA total scores", color = "Group") +
  scale_color_discrete(name = "Group", labels = c("C" = "Clinical", "N" = "Normative")) 

p2 <- ggplot(data, aes(MMSE_total, SAGE_total, color = Status)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Status)) +   
  geom_smooth(method = "lm", color = "gray40", linetype = "dashed", se = FALSE) +
  theme_minimal() +
  labs(x = "MMSE total scores", y = "SAGE total scores", color = "Group") +
  scale_color_discrete(name = "Group", labels = c("C" = "Clinical", "N" = "Normative")) 

p3 <- ggplot(data, aes(MoCA_total, SAGE_total, color = Status)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Status)) +   
  geom_smooth(method = "lm", color = "gray40", linetype = "dashed", se = FALSE) +
  theme_minimal() +
  labs(x = "MoCA total scores", y = "SAGE total scores", color = "Group") +
  scale_color_discrete(name = "Group", labels = c("C" = "Clinical", "N" = "Normative")) 

plot_grid(p1, p2, p3)