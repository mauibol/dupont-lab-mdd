# Differential cell type proportion in MDD vs CTRL neurogenic clusters

library(Seurat)
library(dplyr)
library(lme4)
library(miloR)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)
library(BiocParallel)
library(tidyr)
library(lme4)         # glmer (binomial)
library(glmmTMB)      # betabinomial (optional)
library(broom.mixed)  # tidy mixed model outputs
library(purrr)
library(stringr)


##################### ON FULL #################################################
# Cell-level metadata
meta <- multi@meta.data

meta$RIN_w <- dplyr::case_when(
  meta$RIN < 7 ~ 7,
  .default = meta$RIN
)

meta$PMI_w <- dplyr::case_when(
  meta$PMI > 23 ~ 23,
  .default = meta$PMI
)

expected_cols <- c("official_annotation_short","Study.ID","sample","Diagnosis","Age","batch","RIN_w","PMI_w","Sex")
missing_cols  <- setdiff(expected_cols, names(meta))
if (length(missing_cols)) {
  message("Note: Missing columns skipped: ", paste(missing_cols, collapse = ", "))
}

# Coerce key factors; set NONE as reference if present
meta <- meta %>%
  dplyr::mutate(
    official_annotation_short = factor(.data$official_annotation_short),
    Study.ID         = factor(.data$Study.ID),
    sample      = factor(.data$sample),
    Age         = as.numeric(.data$Age),
    RIN_w         = as.numeric(.data$RIN_w),
    PMI_w         = as.numeric(.data$PMI_w),
    Sex         =  factor(.data$Sex),
    batch       = factor(.data$batch),
    Diagnosis   = if ("Diagnosis" %in% names(.)) factor(.data$Diagnosis) else factor("NONE"),
    Diagnosis   = if ("NONE" %in% levels(Diagnosis)) relevel(Diagnosis, ref = "NONE") else Diagnosis,
  )

# Per-sample cell-type counts
cell_counts <- meta %>%
  dplyr::group_by(Study.ID, sample, Diagnosis, official_annotation_short, Age, batch, RIN_w, PMI_w,RIN,PMI, Sex) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

# Total cells per (donor, sample)
totals <- meta %>%
  dplyr::group_by(Study.ID, sample, Diagnosis) %>%
  dplyr::summarise(total_cells = n(), .groups = "drop")

df <- left_join(cell_counts, totals, by = c("Study.ID", "sample", "Diagnosis"))
if ("NONE" %in% levels(df$Diagnosis)) df$Diagnosis <- relevel(df$Diagnosis, ref = "NONE")

# Get all cell types
cell_types <- unique(df$official_annotation_short)

results <- map_dfr(cell_types, function(ct) {
  dat <- df %>% filter(official_annotation_short == ct)
  
  # Skip if only one diagnosis group or no data
  if (length(unique(dat$Diagnosis)) < 2 || sum(dat$n_cells) == 0) return(NULL)
  
  fit <- glmmTMB(
    cbind(n_cells, total_cells - n_cells) ~ Diagnosis + scale(Age) + scale(RIN) + scale(PMI) + Sex + batch + (1 | Study.ID),
    family = betabinomial(link = "logit"),
    data = dat
  )
  
  tidy(fit, effects = "fixed") %>%
    mutate(
      annotations = ct,
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error)
    ) %>%
    dplyr::select(annotations, term, estimate, std.error, OR, CI_low, CI_high, p.value)
  
})
# Adjust p-values (FDR)
results <- results %>%
  mutate(adj_p_BH = p.adjust(p.value, method = "BH")) %>%
  arrange(adj_p_BH)

print(results %>% filter(term != '(Intercept)'), n=200)


plot_data <- df %>%
  group_by(sample, Study.ID, Diagnosis, official_annotation_short) %>%
  dplyr::summarise(prop = sum(n_cells) / sum(total_cells), .groups = "drop")

summary_data <- plot_data %>%
  group_by(official_annotation_short, Diagnosis) %>%
  dplyr::summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop   = sd(prop, na.rm = TRUE),
    .groups = "drop"
  )

summary_data
#################### On imgns_sub ##########################################
library(dplyr)
library(tidyr)
library(lme4)
library(broom.mixed)
library(glmmTMB)


# Cell-level metadata
meta <- neuro_sub@meta.data

meta$RIN_w <- dplyr::case_when(
  meta$RIN < 7 ~ 7,
  .default = meta$RIN
)

meta$PMI_w <- dplyr::case_when(
  meta$PMI > 23 ~ 23,
  .default = meta$PMI
)

expected_cols <- c("annotations","ID4","sample","Diagnosis","Age","batch","RIN_w","PMI_w","Sex")
missing_cols  <- setdiff(expected_cols, names(meta))
if (length(missing_cols)) {
  message("Note: Missing columns skipped: ", paste(missing_cols, collapse = ", "))
}

# Coerce key factors; set NONE as reference if present
meta <- meta %>%
  dplyr::mutate(
    annotations = factor(.data$annotations),
    ID4         = factor(.data$ID4),
    sample      = factor(.data$sample),
    Age         = as.numeric(.data$Age),
    RIN_w         = as.numeric(.data$RIN_w),
    PMI_w         = as.numeric(.data$PMI_w),
    Sex         =  factor(.data$Sex),
    batch       = factor(.data$batch),
    Diagnosis   = if ("Diagnosis" %in% names(.)) factor(.data$Diagnosis) else factor("NONE"),
    Diagnosis   = if ("NONE" %in% levels(Diagnosis)) relevel(Diagnosis, ref = "NONE") else Diagnosis,
  )

# Per-sample cell-type counts
cell_counts <- meta %>%
  dplyr::group_by(ID4, sample, Diagnosis, annotations, Age, batch, RIN_w, PMI_w,RIN,PMI, Sex) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

# Total cells per (donor, sample)
totals <- meta %>%
  dplyr::group_by(ID4, sample, Diagnosis) %>%
  dplyr::summarise(total_cells = n(), .groups = "drop")

df <- left_join(cell_counts, totals, by = c("ID4", "sample", "Diagnosis"))
if ("NONE" %in% levels(df$Diagnosis)) df$Diagnosis <- relevel(df$Diagnosis, ref = "NONE")

# Get all cell types
cell_types <- unique(df$annotations)

# Run model for each cell type and bind results
results_aic <- map_dfr(cell_types, function(ct) {
  dat <- df %>% filter(annotations == ct)
  
  # Skip if only one diagnosis group or no data
  if (length(unique(dat$Diagnosis)) < 2 || sum(dat$n_cells) == 0) return(NULL)
  
  # log normalized linear
  fit_log <- lmer(
    prop_log ~ Diagnosis + scale(Age) + scale(RIN) + scale(PMI) + Sex + batch +(1 | ID4),
    data = dat
  )
  
  # arcsine linear
  fit_arc <- lmer(
    prop_arcsine ~ Diagnosis + scale(Age) + scale(RIN) + scale(PMI) + Sex + batch +(1 | ID4),
    data = dat
  )
  
  # geta binomial
  fit_betabin <- glmmTMB(
    cbind(n_cells, total_cells - n_cells) ~ Diagnosis + scale(Age) + scale(RIN) + scale(PMI) + Sex + batch + (1 | ID4),
    family = betabinomial(link = "logit"),
    data = dat
  )
  
  
  aic_log <- AIC(fit_log)
  aic_betabin <- AIC(fit_betabin)
  aic_arc <- AIC(fit_arc)
  
  delta_aic <- aic_log - aic_arc
  
  tibble(
    annotations = ct,
    AIC_log = aic_log,
    AIC_arc = aic_arc,
    delta_AIC = delta_aic,
    better_model = ifelse(delta_aic > 2, "log arc",
                          ifelse(delta_aic < -2, "log linear", "similar"))
  )
})

results_aic %>%
  arrange(AIC_arc)


# Run beta binomial model
results <- map_dfr(cell_types, function(ct) {
  dat <- df %>% filter(annotations == ct)
  
  # Skip if only one diagnosis group or no data
  if (length(unique(dat$Diagnosis)) < 2 || sum(dat$n_cells) == 0) return(NULL)
  
  fit <- glmmTMB(
    cbind(n_cells, total_cells - n_cells) ~ Diagnosis + scale(Age) + scale(RIN) + scale(PMI) + Sex + batch  + (1 | ID4),
    family = betabinomial(link = "logit"),
    data = dat
  )
  
  tidy(fit, effects = "fixed") %>%
    mutate(
      annotations3 = ct,
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error)
    ) %>%
    dplyr::select(annotations3, term, estimate, std.error, OR, CI_low, CI_high, p.value)
  
})
# Adjust p-values (FDR)
results <- results %>%
  mutate(adj_p_BH = p.adjust(p.value, method = "BH")) %>%
  arrange(adj_p_BH)

print(results %>% filter(term != '(Intercept)'), n=30)



### PLOT DATA #####

library(dplyr)
library(ggplot2)
library(scales)
library(stringr)
library(dplyr)
library(ggplot2)


plot_data<- df %>%
  dplyr::filter(annotations %in% keep_clusters) %>% 
  group_by(sample, ID4, Diagnosis, annotations) %>%
  dplyr::summarise(prop = sum(n_cells) / sum(total_cells), .groups = "drop")


summary_data <- plot_data %>%
  group_by(annotations, Diagnosis) %>%
  dplyr::summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop   = sd(prop, na.rm = TRUE),
    .groups = "drop"
  )

top_samples <- plot_data %>%
  group_by(annotations) %>%
  slice_max(order_by = prop, n = 6, with_ties = FALSE) %>%
  ungroup()


plot_data$experimenter <- substr(plot_data$sample, 1,2)


#Plot
p <- results %>% filter(term %in% c("DiagnosisMDD"))  %>% filter(annotations %in% c('NSC.a','NSC.b','INP','NB','ImGC.1','ImGC.2','mGC')) %>% 
  mutate(sig = adj_p_BH < 0.05) %>%
  ggplot(aes(x = reorder(annotations, OR), y = OR, color = sig)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  scale_y_log10() +
  coord_flip() +
  scale_color_manual(values = c("grey50", "firebrick")) +
  labs(
    y = "Odds ratio (MDD / CTRL)",
    x=''
  ) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")


# Per donor cell type proportions
plot_data %>%
  ggplot(aes(x = Diagnosis, y = prop, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 1.5) +
  facet_wrap(~official_annotation_short, scales = "free_y") +
  labs(
    y = "Proportion per sample",
    x = "",
    title = "Per-sample cell-type proportions"
  ) +
  theme_bw(base_size = 12)

ggplot(plot_data, aes(x = reorder(official_annotation_short, prop, FUN = mean), y = prop, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, width = 0.6, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(alpha = 0.6, size = 1.5, position = position_dodge(width = 0.8)) +
  labs(
    y = "Proportion per sample",
    x = "Anatomical Region",
    title = "Per-sample cell-type proportions"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# Main plot
plot_data%>%
  filter(annotations %in% c("NSC.a","NSC.b","INP","NB","ImGC.1","ImGC.2")) %>%
  ggplot(aes(x = Diagnosis, y = prop, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 1.5) +
  
  # Add labels for top samples
  geom_text(
    data = top_samples_sub,
    aes(label = sample),
    vjust = -0.6,
    hjust = 0.5,
    size = 3,
    color = "black"
  ) +
  
  facet_wrap(~annotations, scales = "free_y") +
  geom_text(
    data = results_full %>% filter(term %in% c("DiagnosisMDD"))%>%
      filter(adj_p_BH < 0.1) %>%
      mutate(sig = case_when(
        adj_p_BH < 0.001 ~ "***",
        adj_p_BH < 0.01 ~ "**",
        adj_p_BH < 0.05 ~ "*"
      )),
    aes(x = 1.5, y = Inf, label = sig),
    vjust = 1.2, size = 5, color = "firebrick", inherit.aes = FALSE
  ) +
  labs(
    y = "Proportion per sample",
    x = "",
    title = "Per-sample cell-type proportions"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11)
  )

