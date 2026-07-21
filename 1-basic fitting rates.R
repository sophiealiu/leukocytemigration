# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Purpose: simple estimation of CD4 PBS forward and reverse rates
#        - these will be manually put into MATLAB file (a bit inelegant).
# -----------------------------------------------------------------------------

library(broom)
library(dplyr)
library(ggplot2)
library(janitor)
library(purrr)
library(readxl)
library(tidyverse)


# explicit path, set to working directory
leukdatadir <- "I:/Hu Lab/Sophie/2. Leukocyte migration/data"
dfmig <- read_excel(file.path(leukdatadir, "estim_0720.xlsx"), sheet = "summary")


# -----------------------------------------------------------------------------
# simple assumption u(t) = u0exp[-krev*t], reverse direction rate constant
# once labeled it cannot lose its label
# making it into format R nls understands
df_long <- dfmig %>%
  pivot_longer(
    cols = starts_with("t"),
    names_to = "time",
    values_to = "percent"
  ) %>%
  mutate(time = as.numeric(sub("^t", "", time)))

tissues <- unique(df_long$sample_type)
tissues <- tissues[tissues != "blood"]

# fitting nonlinear least sq, pct labeled response
k_revs <- vector("list", length(tissues))
names(k_revs) <- tissues

for (i in seq_along(tissues)) {
  tissue <- tissues[i]
  model <- nls(
    percent ~ exp(-k * time),         # treating as fractions u(0) = 1
    data = df_long %>% 
      filter(sample_type == tissue),
    start = list(k = 0.001))          # required to have initial guess. starting small based on inspection
  k_revs[[i]] <- coef(model)["k"]
}
print(k_revs)


# -----------------------------------------------------------------------------
# monte-carlo for singular tissue
set.seed(42)                          # the answer to life, the universe, and everything
plots <- vector("list", length(tissues))
names(plots) <- tissues

for (i in seq_along(tissues)) {
  
  tissue <- tissues[i]
  
  df_sub <- df_long %>%
    filter(sample_type == tissue)
  
  model <- nls(
    percent ~ exp(-k * time),
    data = df_sub,
    start = list(k = 0.001)
  )
  
  k_hat <- coef(model)["k"]
  k_se  <- summary(model)$parameters["k", "Std. Error"]
  
  k_sim <- rnorm(
    nsim,
    mean = k_hat,
    sd = k_se
  )
  
  times <- seq(0, max(df_sub$time), length.out = 200)
  
  pred_mat <- outer(
    k_sim,
    times,
    function(k, t) exp(-k * t)
  )
  
  model2 <- exp(-k_hat * times)
  
  lower <- apply(pred_mat, 2, quantile, 0.025)
  upper <- apply(pred_mat, 2, quantile, 0.975)
  
  montC_df <- data.frame(
    time = times,
    fit = model2,
    lower = lower,
    upper = upper
  )
  
  plots[[i]] <- ggplot() +
    geom_ribbon(
      data = montC_df,
      aes(time, ymin = lower, ymax = upper),
      alpha = 0.2
    ) +
    geom_line(
      data = montC_df,
      aes(time, fit),
      linewidth = 1.2
    ) +
    geom_point(
      data = df_sub,
      aes(time, percent),
      alpha = 0.7
    ) +
    labs(title = tissue)
}


# -----------------------------------------------------------------------------
# kfor = krev(counts tissue/counts blood) due to conservation of volume
k_fors <- vector("list", length(tissues))
names(k_fors) <- tissues

for (i in seq_along(tissues)) {
  tissue <- tissues[i]
  ratio <- dfmig %>%
    filter(sample_type == tissue) %>%                     # averaged across all mice
    summarise(mean_ratio = mean(cell_count / blood, na.rm = TRUE))
  
  calc <- k_revs[[i]] * ratio
  k_fors[[i]] <- calc
}
print(k_fors)


# -----------------------------------------------------------------------------
# visualization
order <- c(
  "bone_marrow", 
  "peyer", 
  "axillary_LN", 
  "iliac_LN", 
  "ing_LN", 
  "mes_LN", 
  "spleen_white", 
  "spleen_red"
)

df_plot <- data.frame(
  tissue = names(k_fors),
  k_value = unlist(k_fors),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
  row.names = NULL
)

df_plot$tissue <- factor(df_plot$tissue, levels = order)

ggplot(df_plot, aes(x = k_value, y = tissue)) +
  geom_col() +
  labs(title = "rate constants of entry for different tissues",
       x = "forward rate",
       y = "tissue")

