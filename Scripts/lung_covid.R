# SETTING PREREQUISITES ========================================================
library(tidyverse)
library(ggpubr)
library(rstatix)
library(patchwork)
library(gt)
library(gtsummary)
library(broom)
library(magick)
library(grid)

output_dir  <- "output"
figure_dir  <- file.path(output_dir, "figures")
table_dir   <- file.path(output_dir, "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

antigen_order <- c("Wildtype", "Alpha", "Beta", "Gamma", "Delta", "Omicron")

cb_antigen <- c(
  "Wildtype" = "#000000",
  "Alpha"    = "#E69F00",
  "Beta"     = "#56B4E9",
  "Gamma"    = "#009E73",
  "Delta"    = "#0072B2",
  "Omicron"  = "#D55E00"
)

hiv_palette <- c("HIV-" = "#1f77b4", "HIV+" = "#d62728")
comp_palette <- c("BAL" = "#1b9e77", "PB" = "#7570b3")
iso_palette <- c("IgA" = "#1b9e77", "IgG" = "#d95f02")

format_p <- function(p) {
  ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
}

publication_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0),
      plot.title.position = "plot",
      strip.text = element_text(face = "bold", size = 18),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 18),
      panel.background = element_rect(fill = "white", colour = "black"),
      panel.grid.major = element_line(color = "grey70", linewidth = 0.1),
      panel.grid.minor = element_line(color = "grey85", linewidth = 0.08)
    )
}

# TABLE 1 ======================================================================
table1_data <- read_csv("data/Table1_data.csv") %>%
  dplyr::mutate(
    `HIV status` = factor(
      `HIV status`,
      levels = c("HIV-uninfected", "PLHIV ART>1yr")
    ),
    Sex = factor(Sex),
    `Vaccination status` = factor(`Vaccination status`),
    `Number of doses` = factor(`Number of doses`),
    `Primary vaccine type (Dose 1)` = factor(`Primary vaccine type (Dose 1)`),
    `Primary vaccine type (Dose 2)` = factor(`Primary vaccine type (Dose 2)`),
    `HIV viral load` = factor(
      `HIV viral load`,
      levels = c("Not Detected", "Detected")
    )
  )

Table1_final <- table1_data %>%
  select(
    `HIV status`,
    Age,
    Sex,
    `Vaccination status`,
    `Number of doses`,
    `Primary vaccine type (Dose 1)`,
    `Primary vaccine type (Dose 2)`,
    `Time since last vaccination (months)`,
    `Absolute CD4`,
    `HIV viral load`
  ) %>%
  tbl_summary(
    by = `HIV status`,
    type = list(`HIV viral load` ~ "categorical"),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    percent = "column",
    digits = list(all_categorical() ~ 0),
    missing = "no",
    label = list(
      Age ~ "Age (years)",
      Sex ~ "Sex",
      `Vaccination status` ~ "Vaccination status",
      `Number of doses` ~ "Number of doses",
      `Primary vaccine type (Dose 1)` ~ "Primary vaccine type (Dose 1)",
      `Primary vaccine type (Dose 2)` ~ "Primary vaccine type (Dose 2)",
      `Time since last vaccination (months)` ~
        "Time since last COVID-19 vaccination (months)",
      `Absolute CD4` ~ "Absolute CD4 cell count (cells/µL)",
      `HIV viral load` ~ "HIV viral load"
    )
  ) %>%
  add_p(
    include = c(
      Age,
      Sex,
      `Vaccination status`,
      `Number of doses`,
      `Primary vaccine type (Dose 1)`,
      `Primary vaccine type (Dose 2)`,
      `Time since last vaccination (months)`
    ),
    test = list(
      Age ~ "wilcox.test",
      `Time since last vaccination (months)` ~ "wilcox.test",
      Sex ~ "fisher.test",
      `Vaccination status` ~ "chisq.test",
      `Number of doses` ~ "chisq.test",
      `Primary vaccine type (Dose 1)` ~ "fisher.test",
      `Primary vaccine type (Dose 2)` ~ "fisher.test"
    ),
    test.args = all_tests("wilcox.test") ~ list(exact = FALSE),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        p.value = ifelse(
          label %in% c(
            "Absolute CD4 cell count (cells/µL)",
            "HIV viral load"
          ),
          NA_real_,
          p.value
        )
      )
  ) %>%
  bold_labels() %>%
  modify_header(label = "**Characteristic**")

Table1_gt <- Table1_final %>%
  as_gt() %>%
  tab_options(
    table.font.size = 12,
    data_row.padding = px(4)
  )

gtsave(Table1_gt, filename = file.path(table_dir, "Table1.html"))
Table1_gt

# TABLE 2 ======================================================================

regression_df_std <- read_csv("data/Table2_data.csv")
regression_df <- read_csv(file.path("data/Table2_data.csv"))

# Standardize continuous variables
regression_df_std <- regression_df %>%
  mutate(across(
    c(BAL_IgG, BAL_IgA, Serum_IgG, Serum_IgA, NLF_IgG, NLF_IgA,
      Serum_Neutralization, NLF_Neutralization, Age),
    ~ as.numeric(scale(.))
  ))

model_NLF <- lm(
  NLF_Neutralization ~ NLF_IgA + NLF_IgG + Serum_IgG + Serum_IgA + BAL_IgG + BAL_IgA + Age + Sex + HIV_Status + Vacc_Status,
  data = regression_df_std
)

model_serum <- lm(
  Serum_Neutralization ~ Serum_IgG + Serum_IgA + BAL_IgG + BAL_IgA + NLF_IgG + NLF_IgA + Age + Sex + HIV_Status + Vacc_Status,
  data = regression_df_std
)

model_BAL <- lm(
  BAL_IgG ~ Serum_IgG + Serum_IgA + BAL_IgA + NLF_IgG + NLF_IgA + Age + Sex + HIV_Status + Vacc_Status,
  data = regression_df_std
)

make_model_table <- function(model, estimate_name, p_name) {
  tidy(model, conf.int = TRUE) %>%
    mutate(
      "{estimate_name}" := sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
      "{p_name}" := ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    ) %>%
    select(term, all_of(c(estimate_name, p_name)))
}

tbl_NLF <- make_model_table(model_NLF, "NLF", "p_NLF")
tbl_serum <- make_model_table(model_serum, "Serum", "p_Serum")
tbl_BAL <- make_model_table(model_BAL, "BAL", "p_BAL")

table_all <- full_join(tbl_NLF, tbl_serum, by = "term") %>%
  full_join(tbl_BAL, by = "term") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Predictor = case_when(
      term == "Age" ~ "Age",
      term == "SexM" ~ "Sex (Male)",
      term == "HIV_StatusHIV+" ~ "HIV+ status",
      term == "Vacc_StatusVaccinated" ~ "Vaccinated",
      term == "Serum_IgG" ~ "Serum IgG binding breadth",
      term == "Serum_IgA" ~ "Serum IgA binding breadth",
      term == "BAL_IgG" ~ "BAL IgG binding breadth",
      term == "BAL_IgA" ~ "BAL IgA binding breadth",
      term == "NLF_IgG" ~ "NLF IgG binding breadth",
      term == "NLF_IgA" ~ "NLF IgA binding breadth",
      TRUE ~ term
    ),
    Group = case_when(
      Predictor %in% c("Age", "Sex (Male)", "HIV+ status", "Vaccinated") ~ "Participant characteristics",
      Predictor %in% c("Serum IgG binding breadth", "Serum IgA binding breadth") ~ "Systemic markers",
      Predictor %in% c("BAL IgG binding breadth", "BAL IgA binding breadth") ~ "Lower airway markers",
      Predictor %in% c("NLF IgG binding breadth", "NLF IgA binding breadth") ~ "Upper airway markers",
      TRUE ~ "Other"
    )
  ) %>%
  select(Group, Predictor, NLF, p_NLF, Serum, p_Serum, BAL, p_BAL) %>%
  mutate(across(c(NLF, p_NLF, Serum, p_Serum, BAL, p_BAL), ~replace_na(.x, "—")))

# Logical order
row_order <- c(
  "Age", "Sex (Male)", "HIV+ status", "Vaccinated",
  "Serum neutralization breadth", "Serum IgG binding breadth", "Serum IgA binding breadth",
  "BAL IgG binding breadth", "BAL IgA binding breadth",
  "NLF neutralization breadth", "NLF IgG binding breadth", "NLF IgA binding breadth"
)

table_all <- table_all %>%
  slice(match(row_order, Predictor)) %>%
  mutate(across(c(NLF, p_NLF, Serum, p_Serum, BAL, p_BAL), ~ ifelse(is.na(.), "—", .))) # Replace NAs with dash for cleaner reading

final_gt <- table_all %>%
  gt(groupname_col = "Group") %>%
  cols_label(
    Predictor = "Predictor",
    NLF = "NLF neutralization breadth\nβ (95% CI)",
    p_NLF = "p",
    Serum = "Serum neutralization breadth\nβ (95% CI)",
    p_Serum = "p",
    BAL = "BAL IgG binding breadth\nβ (95% CI)",
    p_BAL = "p"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = c(NLF, p_NLF), rows = suppressWarnings(as.numeric(sub("<", "", p_NLF))) < 0.05)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = c(Serum, p_Serum), rows = suppressWarnings(as.numeric(sub("<", "", p_Serum))) < 0.05)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = c(BAL, p_BAL), rows = suppressWarnings(as.numeric(sub("<", "", p_BAL))) < 0.05)
  ) %>%
  tab_source_note(
    source_note = paste0(
      "Multivariable linear regression models fitted on n = ",
      nrow(regression_df_std),
      " participants. β values are standardized coefficients."
    )
  ) %>%
  tab_options(table.font.size = 12, data_row.padding = px(4))

gtsave(final_gt, filename = file.path(table_dir, "Table2.html"))
final_gt

# FIGURE 1 =====================================================================
# Fig1a Panel A: binding antibody isotypes
fig1a_data <- read_csv("data/Fig1A_data.csv") %>%
  dplyr::mutate(
    Compartment = factor(Compartment, levels = c("Serum", "BAL", "NLF")),
    Isotype = factor(Isotype, levels = c("IgA", "IgG"))
  )

fig1a_p <- fig1a_data %>%
  group_by(Compartment, Antigens) %>%
  wilcox_test(
    Log10_Value ~ Isotype,
    paired = TRUE
  ) %>%
  ungroup() %>%
  adjust_pvalue(method = "BH") %>%
  add_xy_position(x = "Isotype") %>%
  mutate(
    p_label = ifelse(
      p.adj < 0.001,
      "p < 0.001",
      sprintf("p = %.3f", p.adj)
    )
  )

igg_line <- tibble(
  Antigens = "Nucleocapsid",
  Compartment = factor("Serum", levels = c("Serum", "BAL", "NLF")),
  yintercept = 1
)

Fig1A <- ggplot(fig1a_data, aes(Isotype, Log10_Value, fill = Isotype)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 1.5) +
  geom_hline(
    data = igg_line,
    aes(yintercept = yintercept),
    inherit.aes = FALSE,
    colour = "#d95f02",
    linetype = "dotted",
    linewidth = 1
  ) +
  facet_grid(Antigens ~ Compartment) +
  stat_pvalue_manual(fig1a_p, label = "p_label", tip.length = 0, size = 6) +
  scale_fill_manual(values = iso_palette) +
  scale_y_continuous(breaks = 0:5, limits = c(0, 5), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Isotype", y = "Log10 Antibody Concentrations (WHO BAU/ml)", title = "A.") +
  publication_theme() +
  theme(legend.position = "none")
Fig1A

#Fig1b: Panel B: breadth of binding antibody responses
fig1b_data <- read_csv("data/Fig1B_data.csv") %>%
  dplyr::mutate(
    Antigen = factor(Antigens, levels = antigen_order),
    Compartment = factor(Compartment, levels = c("Serum", "BAL", "NLF")),
    Isotype = factor(Isotype, levels = c("IgA", "IgG"))
  )

fig1b_p <- fig1b_data %>%
  group_by(Compartment, Isotype) %>%
  wilcox_test(Log10_Value ~ Antigen, paired = TRUE, ref.group = "Wildtype") %>%
  ungroup() %>%
  adjust_pvalue(method = "BH") %>%
  add_xy_position(x = "Antigen", step.increase = 0.25) %>%
  mutate(p_label = ifelse(p.adj < 0.001, "p < 0.001", sprintf("p = %.3f", p.adj)))

Fig1B <- ggplot(fig1b_data, aes(Antigen, Log10_Value, fill = Antigen)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.18, alpha = 0.5, color = "black", size = 1.5) +
  facet_grid(Isotype ~ Compartment) +
  stat_pvalue_manual(
    fig1b_p,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0,
    size = 6,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = cb_antigen) +
  scale_y_continuous(breaks = seq(-1, 7, 1), expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(-1, 9)) +
  labs(x = "Antigen", y = "Log10 Antibody Concentrations (MSD AU/ml)", title = "B.") +
  publication_theme() +
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1), legend.position = "none")
Fig1B

# Fig1c: cross-compartment antibody correlations
y_order <- rev(c(
  "NLF IgA",
  "NLF IgG",
  "BAL IgA",
  "BAL IgG",
  "Serum IgA"
))

x_order <- c(
  "Serum IgG",
  "Serum IgA",
  "BAL IgG",
  "BAL IgA",
  "NLF IgG"
)

fig1c_data <- read_csv("data/Fig1C_data.csv") %>%
mutate(
    sig = case_when(
      is.na(p.adj) ~ "",
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = if_else(is.na(r), "", paste0(round(r, 2), "\n", sig)),
    Antigens = factor(Antigens, levels = antigen_order),
    var1 = factor(var1, levels = x_order),
    var2 = factor(var2, levels = y_order)
  )

Fig1C <- ggplot(fig1c_data, aes(var1, var2, fill = r)) +
  geom_tile() +
  geom_text(aes(label = na_if(label, "")), size = 4) +
  scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, low = "#4575b4", high = "#d73027", na.value = "grey90") +
  coord_equal() +
  facet_wrap(~Antigens, ncol = 6) +
  labs(x = NULL, y = NULL, fill = "r", title = "C.") +
  publication_theme(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18))
Fig1C

# Fig1d: ACE2 inhibition across serum and NLF
fig1d_data <- read_csv("data/Fig1D_data.csv") %>%
  dplyr::mutate(
    Antigens = factor(Antigens, levels = antigen_order),
    Compartment = factor(Compartment, levels = c("Serum", "NLF"))
  )

fig1d_p <- fig1d_data %>%
  group_by(Antigens) %>%
  wilcox_test(Neutralization ~ Compartment, paired = TRUE) %>%
  ungroup() %>%
  adjust_pvalue(method = "BH") %>%
  mutate(p_label = ifelse(p.adj < 0.001, "p < 0.001", sprintf("FDR = %.3f", p.adj))) %>%
  left_join(
    fig1d_data %>% group_by(Antigens) %>% summarise(y.position = max(Neutralization, na.rm = TRUE) + 2, .groups = "drop"),
    by = "Antigens"
  ) %>%
  mutate(y.position = pmin(y.position, 100))

Fig1D <- ggplot(fig1d_data, aes(Compartment, Neutralization, fill = Compartment)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "grey20") +
  facet_wrap(~Antigens, nrow = 3) +
  stat_pvalue_manual(fig1d_p, label = "p_label", y.position = "y.position", tip.length = 0, size = 6) +
  scale_fill_manual(values = c("NLF" = "#0072B2", "Serum" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(20, 100, 20), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Compartment", y = "ACE2 inhibition (%)", title = "D.") +
  publication_theme() +
  theme(legend.position = "none")
Fig1D

# Fig1e: binding antibodies versus ACE2 inhibition
fig1e_data <- read_csv("data/Fig1E_data.csv") %>%
  dplyr::mutate(
    Antigens = factor(Antigens, levels = antigen_order),
    Compartment = factor(Compartment, levels = c("Serum", "NLF")),
    Isotype = factor(Isotype, levels = c("IgG", "IgA"))
  )

fig1e_cor <- fig1e_data %>%
  group_by(Antigens, Compartment, Isotype) %>%
  summarise(
    cor = cor(Log10_Value, Neutralization, method = "spearman", use = "complete.obs"),
    p = suppressWarnings(cor.test(Log10_Value, Neutralization, method = "spearman", exact = FALSE)$p.value),
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    p_label = ifelse(p.adj < 0.001, "p < 0.001", sprintf("p = %.3f", p.adj)),
    y_pos = if_else(Isotype == "IgG", 118, 108),
    label = paste0("r = ", round(cor, 2), ", ", p_label)
  )

Fig1E <- ggplot(fig1e_data, aes(Log10_Value, Neutralization, color = Isotype)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Compartment + Antigens, nrow = 2, ncol = 6) +
  geom_text(
    data = fig1e_cor,
    aes(x = -Inf, y = y_pos, label = label, color = Isotype),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.2,
    size = 5
  ) +
  scale_color_manual(values = c("IgG" = "#d95f02", "IgA" = "#1b9e77")) +
  scale_y_continuous(breaks = seq(0, 120, 20), expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(-1, 119)) +
  labs(x = "Log10 Antibody Level (MSD AU/ml)", y = "ACE2 Inhibition (%)", title = "E.") +
  publication_theme()
Fig1E

# FIGURE 2 =====================================================================
fig2a_data <- read_csv("data/Fig2A_data.csv") %>%
  dplyr::mutate(
    Compartment = factor(Compartment, levels = c("Serum", "BAL", "NLF")),
    Antigens = factor(Antigens, levels = c("Nucleocapsid", "RBD", "Wildtype"))
  )

fig2a_p <- fig2a_data %>%
  group_by(Antigens, Compartment) %>%
  wilcox_test(Log10_Value ~ HIV_Status) %>%
  ungroup() %>%
  adjust_pvalue(method = "BH") %>%
  add_xy_position(x = "HIV_Status") %>%
  mutate(p_label = ifelse(p.adj < 0.001, "p < 0.001", sprintf("p = %.3f", p.adj)))

Fig2A <- ggplot(fig2a_data, aes(HIV_Status, Log10_Value, fill = HIV_Status)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 1.5) +
  facet_grid(Antigens ~ Compartment) +
  stat_pvalue_manual(fig2a_p, label = "p_label", xmin = "group1", xmax = "group2", y.position = "y.position", tip.length = 0, size = 6, inherit.aes = FALSE) +
  scale_fill_manual(values = hiv_palette) +
  coord_cartesian(ylim = c(0, 6)) +
  labs(x = "HIV Status", y = "Log10 IgA Antibody Concentrations (WHO BAU/ml)", title = "A.") +
  publication_theme() + theme(legend.position = "none")
Fig2A

fig2b_data <- read_csv("data/Fig2B_data.csv") %>%
  dplyr::mutate(
    Compartment = factor(Compartment, levels = c("Serum", "BAL", "NLF")),
    Antigens = factor(Antigens, levels = c("Nucleocapsid", "RBD", "Wildtype"))
  )
fig2b_p <- fig2b_data %>%
  group_by(Antigens, Compartment) %>%
  wilcox_test(Log10_Value ~ HIV_Status) %>%
  ungroup() %>%
  adjust_pvalue(method = "BH") %>%
  add_xy_position(x = "HIV_Status") %>%
  mutate(p_label = ifelse(p.adj < 0.001, "p < 0.001", sprintf("p = %.3f", p.adj)))

igg_line <- tibble(
  Antigens = "Nucleocapsid",
  Compartment = factor("Serum", levels = c("Serum", "BAL", "NLF")),
  yintercept = 1
)

Fig2B <- ggplot(fig2b_data, aes(HIV_Status, Log10_Value, fill = HIV_Status)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 1.5) +
  geom_hline(data = igg_line, aes(yintercept = yintercept), inherit.aes = FALSE, colour = "#d95f02", linetype = "dotted", linewidth = 1) +
  facet_grid(Antigens ~ Compartment) +
  stat_pvalue_manual(fig2b_p, label = "p_label", xmin = "group1", xmax = "group2", y.position = "y.position", tip.length = 0, size = 6, inherit.aes = FALSE) +
  scale_fill_manual(values = hiv_palette) +
  coord_cartesian(ylim = c(0, 6)) +
  labs(x = "HIV Status", y = "Log10 IgG Antibody Concentrations (WHO BAU/ml)", title = "B.") +
  publication_theme() + theme(legend.position = "none")
Fig2B


fig2c_data <- read_csv("data/Fig2C_data.csv") %>%
  dplyr::mutate(
    Antigens = factor(Antigens, levels = antigen_order),
    Compartment = factor(Compartment, levels = c("Serum", "NLF"))
  )

fig2c_p <- fig2c_data %>%
  group_by(Antigens, Compartment) %>%
  wilcox_test(Neutralization ~ HIV_Status, paired = FALSE) %>%
  ungroup() %>%
  adjust_pvalue(method = "BH") %>%
  mutate(
    p_label = ifelse(p.adj < 0.001, "p < 0.001", sprintf("p = %.3f", p.adj)),
    group1 = "HIV-",
    group2 = "HIV+"
  ) %>%
  left_join(
    fig2c_data %>% group_by(Antigens, Compartment) %>% summarise(y.position = max(Neutralization, na.rm = TRUE) + 2, .groups = "drop"),
    by = c("Antigens", "Compartment")
  ) %>%
  mutate(y.position = pmin(y.position, 100))

Fig2C <- ggplot(fig2c_data, aes(HIV_Status, Neutralization, fill = HIV_Status)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "grey20") +
  facet_grid(Compartment ~ Antigens) +
  stat_pvalue_manual(fig2c_p, label = "p_label", xmin = "group1", xmax = "group2", y.position = "y.position", tip.length = 0, size = 6, inherit.aes = FALSE) +
  scale_fill_manual(values = hiv_palette) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(20, 100, 20)) +
  labs(x = "HIV Status", y = "ACE2 inhibition (%)", title = "C.") +
  publication_theme() + theme(legend.position = "none")
Fig2C

# FIGURE 3 =====================================================================

fig3b_data <- read_csv("data/Fig3B_data.csv") %>%
  dplyr::mutate(
    Stimulant = factor(Stimulant, levels = c("PMA & Ionomycin", "WT Spike", "WT Nucleocapsid")),
    Compartment = factor(Compartment, levels = c("BAL", "PB"))
  )

fig3b_p <- fig3b_data %>%
  group_by(Stimulant) %>%
  wilcox_test(Log_CD4Total ~ Compartment) %>%
  mutate(p_label = format_p(p)) %>%
  add_y_position(fun = "max", step.increase = 0.08)


Fig3B <- ggplot(fig3b_data, aes(Compartment, Log_CD4Total, fill = Compartment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~Stimulant) +
  stat_pvalue_manual(fig3b_p, label = "p_label", tip.length = 0, size = 5) +
  scale_fill_manual(values = comp_palette) +
  labs(title = "B.", y = "Log10[Total CD4 T cell cytokine-producing cells]", x = "Compartment") +
  labs(title = "B.", y = "Log10[Total CD4 T cell \n cytokine-producing cells]", x = "Compartment") +
  publication_theme() + theme(legend.position = "none")
Fig3B

fig3c_data <- read_csv("data/Fig3C_data.csv") %>%
  mutate(
    Stimulant = factor(Stimulant, levels = c("PMA & Ionomycin", "WT Spike", "WT Nucleocapsid")),
    Compartment = factor(Compartment, levels = c("BAL", "PB"))
  )

fig3c_p <- fig3c_data %>%
  group_by(Stimulant) %>%
  wilcox_test(Log_CD8Total ~ Compartment) %>%
  mutate(p_label = format_p(p)) %>%
  add_y_position(fun = "max", step.increase = 0.08)

Fig3C <- ggplot(fig3c_data, aes(Compartment, Log_CD8Total, fill = Compartment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~Stimulant) +
  stat_pvalue_manual(fig3c_p, label = "p_label", tip.length = 0, size = 5) +
  scale_fill_manual(values = comp_palette) +
  labs(title = "C.", y = "Log10[Total CD8 T cell cytokine-producing cells]", x = "Compartment") +
  labs(title = "C.", y = "Log10[Total CD8 T cell \n cytokine-producing cells]", x = "Compartment") +
  publication_theme() + theme(legend.position = "none")
Fig3C

fig3d_data <- read_csv("data/Fig3D_data.csv") %>%
  mutate(
    sig = case_when(
      is.na(p) ~ "",
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = if_else(is.na(r), "", paste0(round(r, 2), "\n", sig))
  )

Fig3D <- ggplot(fig3d_data, aes(var1, var2, fill = r)) +
  geom_tile() +
  geom_text(aes(label = na_if(label, "")), size = 5) +
  scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, low = "#4575b4", high = "#d73027", na.value = "grey90") +
  coord_equal() +
  facet_grid(. ~ Compartment + Stimulant) +
  labs(x = NULL, y = NULL, fill = "r", title = "D.") +
  publication_theme(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16), axis.text.y = element_text(size = 16))
Fig3D

# FIGURE 4 =====================================================================

fig4a_data <- read_csv("data/Fig4A_data.csv") %>%
  mutate(
    Stimulant = factor(Stimulant, levels = c("PMA & Ionomycin", "WT Spike", "WT Nucleocapsid")),
    Compartment = factor(Compartment, levels = c("BAL", "PB"))
  )

fig4a_p <- fig4a_data %>%
  group_by(Stimulant, Compartment) %>%
  wilcox_test(Log_CD4Total ~ HIV_status) %>%
  mutate(p_label = format_p(p)) %>%
  add_y_position(fun = "max", step.increase = 0.08)

Fig4A <- ggplot(fig4a_data, aes(HIV_status, Log_CD4Total, fill = Compartment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~Compartment + Stimulant, ncol = 3) +
  stat_pvalue_manual(fig4a_p, label = "p_label", tip.length = 0, size = 5) +
  scale_fill_manual(values = comp_palette) +
  labs(title = "A.", y = "Log10[Total CD4 T cell cytokine-producing cells]", x = "HIV status") +
  publication_theme() + theme(legend.position = "none")
Fig4A

fig4b_data <- read_csv("data/Fig4B_data.csv") %>%
  mutate(
    Stimulant = factor(Stimulant, levels = c("PMA & Ionomycin", "WT Spike", "WT Nucleocapsid")),
    Compartment = factor(Compartment, levels = c("BAL", "PB"))
  )

fig4b_p <- fig4b_data %>%
  group_by(Stimulant, Compartment) %>%
  wilcox_test(Log_CD8Total ~ HIV_status) %>%
  mutate(p_label = format_p(p)) %>%
  add_y_position(fun = "max", step.increase = 0.08)

Fig4B <- ggplot(fig4b_data, aes(HIV_status, Log_CD8Total, fill = Compartment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~Compartment + Stimulant, ncol = 3) +
  stat_pvalue_manual(fig4b_p, label = "p_label", tip.length = 0, size = 5) +
  scale_fill_manual(values = comp_palette) +
  labs(title = "B.", y = "Log10[Total CD8 T cell cytokine-producing cells]", x = "HIV status") +
  publication_theme() + theme(legend.position = "none")
Fig4B

fig4c_data <- read_csv("data/Fig4C_data.csv")
fig4c_cor <- fig4c_data %>%
  group_by(Stimulant) %>%
  summarise(r = cor(BAL, PB, method = "spearman"), p = cor.test(BAL, PB, method = "spearman", exact = FALSE)$p.value, .groups = "drop") %>%
  mutate(label = paste0("r = ", round(r, 2), ", ", format_p(p)))

Fig4C <- ggplot(fig4c_data, aes(PB, BAL, colour = HIV_status)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  facet_wrap(~Stimulant) +
  geom_text(data = fig4c_cor, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.1, vjust = 1.5, size = 6) +
  scale_colour_manual(values = hiv_palette) +
  labs(title = "C.", x = "PB CD4 response (log10)", y = "BAL CD4 response (log10)", colour = "HIV status") +
  publication_theme()
Fig4C

fig4d_data <- read_csv("data/Fig4D_data.csv")
fig4d_cor <- fig4d_data %>%
  group_by(Stimulant) %>%
  summarise(r = cor(BAL, PB, method = "spearman"), p = cor.test(BAL, PB, method = "spearman", exact = FALSE)$p.value, .groups = "drop") %>%
  mutate(label = paste0("r = ", round(r, 2), ", ", format_p(p)))

Fig4D <- ggplot(fig4d_data, aes(PB, BAL, colour = HIV_status)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  facet_wrap(~Stimulant) +
  geom_text(data = fig4d_cor, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.1, vjust = 1.5, size = 6) +
  scale_colour_manual(values = hiv_palette) +
  labs(title = "D.", x = "PB CD8 response (log10)", y = "BAL CD8 response (log10)", colour = "HIV status") +
  publication_theme()
Fig4D

# FIGURE 5 =====================================================================

#Figure 5 is not generated analytically from data within this document. 
#The individual phenotyping panels were generated externally in SPICE and PESTEL. This section only imports the completed panels, adds the manuscript panel labels, assembles them, and exports the final Figure 5. 
#The underlying Figure 5 data are supplied separately and are not an input to this build step.

Fig5A <- image_read(file.path(image_dir, "Figure 5B-02.png")) |> image_trim()
Fig5B <- image_read(file.path(image_dir, "Figure 5A-02.png")) |> image_trim()

target_width <- max(image_info(Fig5A)$width, image_info(Fig5B)$width)
Fig5A <- image_scale(Fig5A, paste0(target_width, "x"))
Fig5B <- image_scale(Fig5B, paste0(target_width, "x"))

add_panel_label <- function(img, label, width) {
  img <- image_extent(
    img,
    geometry = paste0(width, "x", image_info(img)$height + 100),
    gravity = "south",
    color = "white"
  )
  image_annotate(
    img,
    label,
    size = 70,
    size = 50,
    weight = 700,
    gravity = "northwest",
    location = "+20+10",
    color = "black"
  )
}

Fig5A <- add_panel_label(Fig5A, "A.", target_width)
Fig5B <- add_panel_label(Fig5B, "B.", target_width)
Figure5 <- image_append(c(Fig5A, Fig5B), stack = TRUE)

image_write(Figure5, path = file.path(figure_dir, "Figure5.png"), format = "png")
Figure5