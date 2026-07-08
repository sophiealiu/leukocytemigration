# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 07/08/2026 3:47pm CDT
# Purpose: Estimating forward and reverse rates
# -----------------------------------------------------------------------------

library(broom)
library(dplyr)
library(ggplot2)
library(janitor)
library(purrr)
library(readxl)
library(tidyverse)


leukdatadir <- "I:/Hu Lab/Sophie/2. Leukocyte migration/data"
dfmig <- read_excel(file.path(leukdatadir, "estim_0708.xlsx"), sheet = "summary")


# -----------------------------------------------------------------------------
# u(t) = u0exp[-krev*t], reverse direction rate constant
# making it into format R nls understands
df_long <- dfmig %>%
  pivot_longer(
    cols = c(`0`, `10`, `20`, `30`, `40`),
    names_to = "time",
    values_to = "percent"
  ) %>%
  mutate(time = as.numeric(time))

tissues <- unique(df_long$sample_type)
tissues <- tissues[tissues != "blood"]

# fitting nonlinear least sq, pct labeled response
k_revs <- vector("list", length(tissues))
names(k_revs) <- tissues

for (i in seq_along(tissues)) {
  tissue <- tissues[i]
  model <- nls(
    percent ~ exp(-k * time),
    data = df_long %>% 
      filter(sample_type == tissue),
    start = list(k = 0.001))          # required to have initial guess. starting small based on inspection
  k_revs[[i]] <- coef(model)["k"]
}


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
blood_counts <- dfmig %>%
  filter(sample_type == "blood") %>%
  pull(cell_count)

dfmig2 <- dfmig %>%
  filter(sample_type != "blood") %>%
  group_by(sample_type) %>%
  mutate(
    mouse = row_number(),
    blood_count = blood_counts[mouse],
    ratio = cell_count / blood_count
  ) %>%
  ungroup()

k_fors <- vector("list", length(tissues))
names(k_fors) <- tissues

for (i in seq_along(tissues)) {
  tissue <- tissues[i]
  
  ratio <- dfmig2 %>%
    filter(sample_type == tissue) %>%
    summarise(mean_ratio = mean(ratio)) %>%
    pull(mean_ratio)
  
  calc <- k_revs[[i]] * ratio
  k_fors[[i]] <- calc
}
