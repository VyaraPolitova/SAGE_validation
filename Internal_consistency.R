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

data <- data %>%
  mutate(Gender = dplyr::recode(Gender, "ж" = "Female", "м" = "Male")) %>%
  rename(MMSE_total = 'MMSE-total',
         MoCA_total = 'MoCA-total',
         SAGE_total = 'SAGEtotal',
         `MoCA-Repetition` = `MoCA- Repetition`) 


data_healthy <- data %>%
  dplyr::filter(Status == "N")
data_clinical <- data %>%
  dplyr::filter(Status == "C")

colSums(is.na(data))
colSums(is.na(data_healthy))

#Dataframes with items only
#DOESN'T INCLUDE MoCA-learning and MoCA-recognition
SAGE_items_total <- data %>%
  gtsummary::select(starts_with("SAGE-") & !contains("total"))

MoCA_items_total <- data %>%
  gtsummary::select(starts_with("MoCA-") & !contains("total") & !contains("MoCA-recognition") & !contains("MoCA-learning1") & !contains("MoCA-learning2"))

names(SAGE_items_total)
names(MoCA_items_total)

 
# Item response distributions and variability check
sapply(SAGE_items_total, function(x) length(unique(x[!is.na(x)])))
sapply(MoCA_items_total, function(x) length(unique(x[!is.na(x)])))

sage_perc <- sapply(SAGE_items_total, function(x) {
  tbl <- table(x)
  prop.table(tbl) * 100
})

sage_perc

moca_perc <- sapply(MoCA_items_total, function(x) {
  tbl <- table(x)
  prop.table(tbl) * 100
})

moca_perc

#--- Pearson correlation matrix
cor(SAGE_items_total, use = "pairwise.complete.obs")
cor(MoCA_items_total, use = "pairwise.complete.obs")

#Collapsing SAGE-8 - only 1 person in it
data$`SAGE-8_c` <- ifelse(data$`SAGE-8` >= 2, 2, data$`SAGE-8`)
table(data$`SAGE-8_c`)

SAGE_items_total_v2 <- SAGE_items_total %>%
  mutate(`SAGE-8` = data$`SAGE-8_c`)   # replace original SAGE-8 with collapsed version

# Unidimensionality check
# Parallel analysis
#SAGE
set.seed(123)
sage_pa <- fa.parallel(SAGE_items_total_v2, fa = "fa", cor = "poly", plot = FALSE)
sage_pa$nfact
# MoCA
set.seed(123)
moca_pa <- fa.parallel(MoCA_items_total, fa = "fa", cor = "poly", plot = FALSE)
moca_pa$nfact

#--- Cronbach's alpha ---
sage_alpha_full <- psych::alpha(as.data.frame(SAGE_items_total_v2), check.keys=FALSE)

moca_alpha_full <- psych::alpha(as.data.frame(MoCA_items_total), check.keys=FALSE)

#Item-total correlations, alpha if item dropped
sage_alpha_full$item.stats %>% gtsummary::select(raw.r, r.drop)
sage_alpha_full$alpha.drop
moca_alpha_full$item.stats %>% gtsummary::select(raw.r, r.drop)
moca_alpha_full$alpha.drop

# Recodes each item to contiguous integer levels based on observed values,
# preserving rank order (assumes original coding is already ordinal with higher = better performance)
SAGE_items_recoded <- SAGE_items_total_v2 %>%
  mutate(across(everything(), ~ as.integer(factor(.x, levels = sort(unique(.x))))))
MoCA_items_recoded <- MoCA_items_total %>%
  mutate(across(everything(), ~ as.integer(factor(.x, levels = sort(unique(.x))))))
# Verify recoding worked (contiguous integers, no gaps)
sapply(SAGE_items_recoded, table, useNA = "ifany")
sapply(MoCA_items_recoded, table, useNA = "ifany")

# --- Diagnostic check for polychoric estimation ---
check_polychoric_pairs <- function(items) {
  
  item_names <- names(items)
  failed_pairs <- data.frame(
    item1 = character(),
    item2 = character()
  )
  
  for (i in 1:(length(item_names) - 1)) {
    for (j in (i + 1):length(item_names)) {
      
      result <- tryCatch({
        polychoric(items[, c(item_names[i], item_names[j])])
        NULL
      }, error = function(e) {
        c(item_names[i], item_names[j])
      })
      
      if (!is.null(result)) {
        failed_pairs <- rbind(
          failed_pairs,
          data.frame(
            item1 = result[1],
            item2 = result[2]
          )
        )
      }
    }
  }
  
  failed_pairs
}

sage_failed_pairs <- check_polychoric_pairs(SAGE_items_recoded)
moca_failed_pairs <- check_polychoric_pairs(MoCA_items_recoded)

sage_failed_pairs
moca_failed_pairs

# Polychoric correlations matrix
sage_polychoric_cor <- polychoric(SAGE_items_recoded)
moca_polychoric_cor <- polychoric(MoCA_items_recoded)

# Polychoric alpha
sage_alpha_poly <- psych::alpha(sage_polychoric_cor$rho, n.obs = nrow(SAGE_items_total_v2))
moca_alpha_poly <- psych::alpha(moca_polychoric_cor$rho, n.obs = nrow(MoCA_items_total))

sage_alpha_poly$total$raw_alpha
moca_alpha_poly$total$raw_alpha

# McDonald's omega - general-factor-saturation index
omega_sage <- omega(sage_polychoric_cor$rho,
                    nfactors = 1,
                    n.obs = nrow(SAGE_items_recoded))
omega_moca <- omega(moca_polychoric_cor$rho,
                    nfactors = 1,
                    n.obs = nrow(MoCA_items_total))


print(omega_sage)
print(omega_moca)

# CI - alpha
# SAGE
sage_alpha_est <- sage_alpha_full$total$raw_alpha
sage_alpha_se  <- sage_alpha_full$total$ase
sage_alpha_lower <- sage_alpha_est - 1.96 * sage_alpha_se
sage_alpha_upper <- sage_alpha_est + 1.96 * sage_alpha_se

# MoCA
moca_alpha_est <- moca_alpha_full$total$raw_alpha
moca_alpha_se  <- moca_alpha_full$total$ase
moca_alpha_lower <- moca_alpha_est - 1.96 * moca_alpha_se
moca_alpha_upper <- moca_alpha_est + 1.96 * moca_alpha_se


strip_zero <- function(x) sub("^(-?)0\\.", "\\1.", x)

#--- Tables ---
#Item-level table
sage_table <- data.frame(
  Item = rownames(sage_alpha_full$item.stats),
  r_full = sage_alpha_full$item.stats$r.drop, 
  alpha_drop_full = sage_alpha_full$alpha.drop$raw_alpha
) %>%
  mutate(
    flag_full = ifelse(r_full < .20, "⚠", "")
  )


moca_table <- data.frame(
  Item = rownames(moca_alpha_full$item.stats),
  r_full = moca_alpha_full$item.stats$r.drop,
  alpha_drop_full = moca_alpha_full$alpha.drop$raw_alpha
) %>%
  mutate(
    flag_full = ifelse(r_full < .20, "⚠", "")
  )


format_item_table <- function(df) {
  df %>%
    gt() %>%
    fmt(
      columns = c(r_full, alpha_drop_full),
      fns = function(x) strip_zero(sprintf("%.2f", x))
    ) %>%
    fmt(
      columns = c(r_full, alpha_drop_full),
      fns = function(x) strip_zero(sprintf("%.2f", x))
    ) %>%
    cols_label(
      Item = "Item",
      r_full = md("Corrected item-total *r*"),
      alpha_drop_full = md("α if item is deleted"),
      flag_full = "Flag"
    ) %>%
    tab_footnote(
      footnote = md("⚠ indicates corrected item-total *r* < .20")
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
      
      data_row.padding = gt::px(4),
      table.font.names = "Times New Roman",
      table.font.size = gt::px(12),
      table.width = gt::pct(60)
    ) 
}

format_item_table(sage_table)
format_item_table(moca_table)


#Internal consistency table
alpha_summary <- data.frame(
  Instrument = c("SAGE", "MoCA"),
  Alpha_Pearson = c(
    sprintf("%s [%s, %s]",
            strip_zero(sprintf("%.2f", sage_alpha_est)),
            strip_zero(sprintf("%.2f", sage_alpha_lower)),
            strip_zero(sprintf("%.2f", sage_alpha_upper))),
    sprintf("%s [%s, %s]",
            strip_zero(sprintf("%.2f", moca_alpha_est)),
            strip_zero(sprintf("%.2f", moca_alpha_lower)),
            strip_zero(sprintf("%.2f", moca_alpha_upper)))
)
)

alpha_summary <- alpha_summary %>%
  mutate(
    Alpha_Polychoric = c(
      strip_zero(sprintf("%.2f", sage_alpha_poly$total$raw_alpha)),
      strip_zero(sprintf("%.2f", moca_alpha_poly$total$raw_alpha))
    )
  )

alpha_omega_summary <- alpha_summary %>%
  mutate(
  Omega_Total = c(
    strip_zero(sprintf("%.2f", omega_sage$omega.tot)),
    strip_zero(sprintf("%.2f", omega_moca$omega.tot))
  )
)

alpha_omega_summary %>%
  gt() %>%
  cols_label(
    Instrument = "Instrument",
    Alpha_Pearson = md("Cronbach's α (Pearson correlations) [95% CI]"),
    Alpha_Polychoric = md("Cronbach's α (Polychoric correlations)"),
    Omega_Total = md("McDonald's ω")
  ) %>%
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
    
    table_body.border.bottom.style = "none",
    table_body.hlines.style = "none",
    table_body.vlines.style = "none",
    stub.border.style = "none",
    heading.border.bottom.style = "none",
    
    data_row.padding = gt::px(4),
    table.font.names = "Times New Roman",
    table.font.size = gt::px(12),
    table.width = gt::pct(60)
  )

#Plots
# SAGE
plot_sage <- ggplot(
  sage_table, 
  aes(x = r_full, y = reorder(Item, r_full))
) +
  geom_point(size = 3) +
  geom_text(
    aes(label = strip_zero(sprintf("%.2f", r_full))),
    hjust = -0.3
  ) +
  geom_vline(xintercept = 0.20, linetype = "dashed") +
  scale_x_continuous(
    labels = function(x) strip_zero(sprintf("%.2f", x)),
    expand = expansion(mult = c(0.05, 0.20))
  ) +
  labs(
    x = "Corrected item-total 𝑟",
    y = "Items"
  ) +
  theme_minimal()

# MoCA
plot_moca <- ggplot(
  moca_table, 
  aes(x = r_full, y = reorder(Item, r_full))
) +
  geom_point(size = 3) +
  geom_text(
    aes(label = strip_zero(sprintf("%.2f", r_full))),
    hjust = -0.3
  ) +
  geom_vline(xintercept = 0.20, linetype = "dashed") +
  scale_x_continuous(
    labels = function(x) strip_zero(sprintf("%.2f", x)),
    expand = expansion(mult = c(0.05, 0.20))
  ) +
  labs(
    x = "Corrected item-total 𝑟",
    y = "Items"
  ) +
  theme_minimal()

plot_grid(plot_sage, plot_moca)