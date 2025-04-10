
# FIGURES ----------------------------------------------------------------------
adm1_alm_res <- readRDS("results/estimates/adm1_alm_res.rds")
adm2_alm_res <- readRDS("results/estimates/adm2_alm_res.rds")
adm1_alm_no_cov_res <- readRDS("results/estimates/adm1_alm_no_cov_res.rds")
adm2_alm_no_cov_res <- readRDS("results/estimates/adm2_alm_no_cov_res.rds")

bin_adm1_res <- readRDS("results/estimates/bin_adm1_res.rds")
bin_adm2_res <- readRDS("results/estimates/bin_adm2_res.rds")
bin_nested_adm2_res <- readRDS("results/estimates/bin_nested_adm2_res.rds")
bin_no_cov_adm1_res <- readRDS("results/estimates/bin_no_cov_adm1_res.rds")
bin_no_cov_adm2_res <- readRDS("results/estimates/bin_no_cov_adm2_res.rds")

bbin_adm1_res <- readRDS("results/estimates/bbin_adm1_res.rds")
bbin_adm2_res <- readRDS("results/estimates/bbin_adm2_res.rds")

bbin_nested_adm2_res <- readRDS("results/estimates/bbin_nested_adm2_res.rds")
bbin_no_cov_adm1_res <- readRDS("results/estimates/bbin_no_cov_adm1_res.rds")
bbin_no_cov_adm2_res <- readRDS("results/estimates/bbin_no_cov_adm2_res.rds")
bbin_simframe_adm1_res <- readRDS("results/estimates/bbin_simframe_adm1_res.rds")
bbin_simframe_adm2_res <- readRDS("results/estimates/bbin_simframe_adm2_res.rds")

lono_adm1_res <- readRDS("results/estimates/lono_adm1_res.rds")
lono_adm2_res <- readRDS("results/estimates/lono_adm2_res.rds")
lono_no_cov_adm1_res <- readRDS("results/estimates/lono_no_cov_adm1_res.rds")
lono_no_cov_adm2_res <- readRDS("results/estimates/lono_no_cov_adm2_res.rds")
lono_simframe_adm1_res <- readRDS("results/estimates/lono_simframe_adm1_res.rds")
lono_simframe_adm2_res <- readRDS("results/estimates/lono_simframe_adm2_res.rds")

bin_LGM_res <- readRDS("results/estimates/bin_LGM_res.rds")
lono_LGM_res <- readRDS("results/estimates/lono_LGM_res.rds")
bbin_LGM_res <- readRDS("results/estimates/bbin_LGM_res.rds")
bin_LGM_no_cov_res <- readRDS("results/estimates/bin_LGM_no_cov_res.rds")
lono_LGM_no_cov_res <- readRDS("results/estimates/lono_LGM_no_cov_res.rds")
bbin_LGM_no_cov_res <- readRDS("results/estimates/bbin_LGM_no_cov_res.rds")
lono_LGM_simframe_res <- readRDS("results/estimates/lono_LGM_simframe_res.rds")
bbin_LGM_simframe_res <- readRDS("results/estimates/bbin_LGM_simframe_res.rds")

adm2_est_table <- bind_rows(
  adm2_alm_res$direct.est,
  adm2_alm_res$bym2.model.est |> mutate(method = "Fay-Herriot BYM2"),
  adm2_alm_no_cov_res$bym2.model.est |> mutate(method = "Fay-Herriot BYM2: No cov."),
  bin_adm2_res$bym2.model.est |> mutate(method = "Binomial BYM2"),
  bin_no_cov_adm2_res$bym2.model.est |> mutate(method = "Binomial BYM2: No cov."),
  bin_nested_adm2_res$bym2.model.est,
  lono_adm2_res$lonobinomial.bym2.model.est |> mutate(method = "Lonobinomial BYM2"),
  lono_no_cov_adm2_res$lonobinomial.bym2.model.est |> mutate(method = "Lonobinomial BYM2: No cov."),
  bbin_adm2_res$betabinomial.bym2.model.est |> mutate(method = "Betabinomial BYM2"),
  bbin_nested_adm2_res$betabinomial.bym2.model.est,
  
  bbin_no_cov_adm2_res$betabinomial.bym2.model.est |> mutate(method = "Betabinomial BYM2: No cov."),
  bin_LGM_res$binomial.spde.lgm.est.subdomain |> mutate(method = "Binomial GRF"),
  lono_LGM_res$binomial.spde.lgm.est.subdomain |> mutate(method = "Lonobinomial GRF"),
  bbin_LGM_res$betabinomial.spde.lgm.est.subdomain |> mutate(method = "Betabinomial GRF"),
  bin_LGM_no_cov_res$binomial.spde.lgm.est.subdomain |> mutate(method = "Binomial GRF: No cov."),
  lono_LGM_no_cov_res$binomial.spde.lgm.est.subdomain |> mutate(method = "Lonobinomial GRF: No cov."),
  bbin_LGM_no_cov_res$betabinomial.spde.lgm.est.subdomain |> mutate(method = "Betabinomial GRF: No cov."),
  lono_simframe_adm2_res$lonobinomial.bym2.model.est |> mutate(method = "Lonobinomial BYM2 Sim Frame"),
  bbin_simframe_adm2_res$betabinomial.bym2.model.est |> mutate(method = "Betabinomial BYM2 Sim Frame"),
  lono_LGM_simframe_res$binomial.spde.lgm.est.subdomain |> mutate(method = "Lonobinomial GRF Sim Frame"),
  bbin_LGM_simframe_res$betabinomial.spde.lgm.est.subdomain |> mutate(method = "Lonobinomial GRF Sim Frame"),
) |>
  left_join(adm2_alm_res$direct.est |>
              select(domain, median, var) |>
              rename(direct.est = median, 
                     direct.est.var = var),
            by = "domain") |>
  mutate(CI = paste0("(", round(lower, 3), ", ", round(upper, 3), ")"))
write.csv(adm2_est_table, "results/estimates/adm2_est.csv")
write.csv(adm2_est_table, "internal/adm2_est.csv")


adm2_est_table <- read.csv("internal/adm2_est.csv")
adm1_est_table <- read.csv("internal/adm1_est.csv")
selected_methods <- 
  c("Direct",
    "Fay-Herriot BYM2",
    "Binomial BYM2",
    "Betabinomial BYM2",
    "Binomial GRF",
    "Betabinomial GRF")
adm1_est_table <- adm1_est_table |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods))
adm2_est_table <- adm2_est_table |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods))


my_theme <- 
  theme(legend.position="bottom",
        strip.text = element_text(size=20),
        legend.text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(size=20),
        legend.key.width = unit(1.5,"cm"))

#### ADMIN 1 ####
plot_sf <- poly_adm1 |>
  left_join(adm1_est_table, by = "domain")
reordered_methods <- 
  c("Direct",
    "Binomial BYM2",
    "Binomial GRF",
    "Fay-Herriot BYM2",
    "Betabinomial BYM2",
    "Betabinomial GRF")
plot_sf <- plot_sf |>
  mutate(method = factor(method, levels = reordered_methods)) |>
  mutate(CV = sqrt(var) / median * 100)
reordered_methods_long <- 
  c("Direct",
    "Binomial BYM2",
    "Betabinomial BYM2",
    "Fay-Herriot BYM2",
    "Binomial GRF",
    "Betabinomial GRF")


#### **** Estimates **** ####
ggplot(plot_sf) + 
  geom_sf(aes(fill = median)) +
  facet_wrap(~method) +
  scale_fill_whitebox_c(palette = "viridi", direction = -1, limits = c(0, .33)) +
  theme_minimal()  + my_theme 
ggsave("results/figures/Zambia_adm1_est_map.png",
       width = 12.5, height = 8)
#### **** CVs **** ####
plot_sf <- plot_sf |>
  mutate(CV = sqrt(var) / median * 100)
ggplot(plot_sf) + 
  geom_sf(aes(fill = CV)) +
  facet_wrap(~method) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "CV (%)") +
  theme_minimal()  + my_theme 
ggsave("results/figures/Zambia_adm1_cv_map.png",
       width = 12.5, height = 8)


cv_table <- adm1_est_table |>
  filter(method == "Direct") |>
  mutate(CV = sqrt(var) / median * 100) |>
  dplyr::select(domain, median, CV) |>
  arrange(CV) 

cv_table |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3) |>
  writeLines(paste0("results/figures/zmb_hiv_adm1_cv_table.tex"))




#### **** Interval lengths **** ####
ggplot(plot_sf) + 
  geom_sf(aes(fill = upper - lower)) +
  facet_wrap(~method) +
  scale_fill_whitebox_c(palette = "bl_yl_rd") +
  labs(fill = "CI length") + 
  theme_minimal()  + my_theme
ggsave("results/figures/Zambia_adm1_int_length_map.png",
       width = 12.5, height = 8)


#### **** Credible intervals **** ####
direct_est <- adm1_est_table |> dplyr::filter(method == "Direct") |> dplyr::filter(!is.na(mean))
sorted_levels <- direct_est$domain[order(direct_est$mean, decreasing = T)] 
sorted_levels <- c(head(sorted_levels, n = 5), tail(sorted_levels, n = 5))
group_label <- c(rep("Highest direct estimates", 5), rep("Lowest direct estimates", 5))
# top 5 and bottom 5 
int_table <- adm1_est_table |>
  filter(domain %in% sorted_levels) |>
  mutate(group = group_label[match(domain, sorted_levels)])
int_table$domain <- factor(int_table$domain, 
                           levels = sorted_levels)
int_table <- int_table |> dplyr::arrange(domain, method)
ggplot2::ggplot(int_table, ggplot2::aes(x = domain, y = median, color = method)) +
  ggplot2::facet_wrap(~group, scales="free_x", ncol = 1) + 
  ggplot2::geom_point(position = ggplot2::position_dodge(width = .6)) + 
  ggplot2::geom_linerange(ggplot2::aes(x = domain, ymin = lower, ymax = upper), 
                          position = ggplot2::position_dodge(width = .6)) + 
  ggplot2::scale_color_discrete(name = "Method") + 
  ggplot2::ylab("Small Area Estimate") + ggplot2::xlab("Domain") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1),
                 legend.position="bottom",
                 legend.key.width = unit(.3,"cm"),
                 strip.text = element_text(size=14),
                 legend.text = element_text(size=9.5)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("results/figures/Zambia_adm1_interval_plots.png",
       width = 6, height = 6)
#### **** Ranking plots **** ####

# code for sampling from asymptotic distribution of direct estimates
# asymptotic variance
direct_est_sample <- matrix(rnorm(10 * 1000), nrow = 10, ncol = 1000)
direct_est_sample <-
  diag(sqrt(adm1_alm_res$direct.est$var)) %*% direct_est_sample 
direct_est_sample <- direct_est_sample + adm1_alm_res$direct.est$mean
compareEstimates(adm1_alm_res, direct_est_sample, "Asymptotic distribution of direct estimates") 
ggsave("results/figures/Zambia_adm1_direct_est_heatmap.png",
       width = 7, height = 8)

compareEstimates(adm1_alm_res, adm1_alm_res$bym2.model.sample, "Area level model (BYM2) posterior")
ggsave("results/figures/Zambia_adm1_alm_bym2_heatmap.png",
       width = 7, height = 8)
compareEstimates(adm1_alm_res, bbin_adm1_res$betabinomial.bym2.model.sample, "Betabinomial model (BYM2) posterior")
ggsave("results/figures/Zambia_adm1_bbin_bym2_heatmap.png",
       width = 7, height = 8)
compareEstimates(adm1_alm_res, bbin_LGM_res$betabinomial.spde.lgm.sample, "Betabinomial model (SPDE) posterior")
ggsave("results/figures/Zambia_adm1_bbin_spde_heatmap.png",
       width = 7, height = 8)


#### **** Estimate tables **** ####

out_table <- adm1_est_table |>
  filter(complete.cases(median)) |>
  mutate(est = paste0(round(median, 2), " (", round(lower, 2),", ", round(upper, 2), ")")) |>
  dplyr::select(est, median, method, domain) |>
  pivot_wider(values_from = c("median", "est"),
              names_from = method,
              names_glue = "{method}_{.value}") |>
  arrange(Direct_median) |>
  dplyr::select(domain, contains("est")) 
names(out_table) <- 
  c("Province",
    stringr::str_replace(names(out_table)[-1], "_est", ""))

out_table |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3) |>
  writeLines(paste0("results/figures/zmb_hiv_adm1_est_table.tex"))




#### **** Scatter Plots **** ####

ggplot(adm1_est_table, aes(x = direct.est, y = mean, color = method)) +
  geom_point(size = 3) +
  geom_abline(slope = 1) +
  xlab("Direct Estimate") +
  ylab("Model Estimate") +
  theme_minimal()  + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position="bottom",
    panel.grid.minor = element_blank(),
    legend.key.size = unit(.5, "cm"),
    legend.text = element_text(size=14),
    legend.key.width = unit(.5,"cm"))
ggsave("results/figures/Zambia_adm1_scatter.png",
       width = 10, height = 8)

ggplot(adm1_est_table, aes(x = direct.est, y = mean, color = method)) +
  geom_point(size = 3, show.legend = F) +
  geom_abline(slope = 1) +
  geom_linerange(aes(ymin = lower, ymax = upper), show.legend = F) + 
  facet_wrap(~method) + 
  xlab("Direct Estimate") +
  ylab("Model Estimate") +
  theme_minimal()  +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position="bottom",
    panel.grid.minor = element_blank(),
    legend.key.size = unit(.5, "cm"),
    legend.text = element_text(size=14),
    legend.key.width = unit(.5,"cm"))
ggsave("results/figures/Zambia_adm1_scatter_bars.png",
       width = 10, height = 8)

#### **** Binomial vs. Betabinomial **** ####

bbin_grf_comp <- adm1_est_table |>
  filter(method %in% c("Binomial GRF", "Betabinomial GRF")) |>
  mutate(int_length = upper - lower,
         method = str_replace(method, " GRF", "")) |>
  select(domain, method, int_length) |>
  pivot_wider(names_from = method, values_from = int_length) |>
  mutate(model = "GRF")

bbin_bym2_comp <- adm1_est_table |>
  filter(method %in% c("Binomial BYM2", "Betabinomial BYM2")) |>
  mutate(int_length = upper - lower,
         method = str_replace(method, " BYM2", "")) |>
  mutate(int_length = upper - lower) |>
  select(domain, method, int_length) |>
  pivot_wider(names_from = method, values_from = int_length) |>
  mutate(model = "BYM2")

bbin_comp <- bind_rows(bbin_bym2_comp, bbin_grf_comp)

ggplot(bbin_comp, aes(x = Binomial, y = Betabinomial)) + 
  facet_wrap(~model) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  ggtitle("Interval lengths (Admin-1)") +
  theme_minimal(base_size = 20)
ggsave("results/figures/Zambia_bin_vs_bbin_int_length_adm1.pdf",
       width = 10, height = 6)
bbin_comp |> group_by(model) |>
  summarize(mean_int_length_adm1_binomial = mean(Binomial),
            mean_int_length_adm1_betabinomial = mean(Betabinomial))


#### ADMIN 2 ####
plot_sf <- poly_adm2 |>
  left_join(adm2_est_table, by = "domain") 
plot_sf <- plot_sf |>
  mutate(method = factor(method, levels = reordered_methods))
#### **** Estimates **** ####
ggplot(plot_sf) + 
  geom_sf(aes(fill = median)) +
  facet_wrap(~method,nrow = 2) +
  scale_fill_whitebox_c(palette = "viridi", direction = -1, limits = c(0, .33)) +
  theme_minimal()  + my_theme
ggsave("results/figures/Zambia_adm2_est_map.png",
       width = 12.5, height = 8)
#### **** CVs **** ####
plot_sf <- plot_sf |>
  mutate(CV = sqrt(var) / median * 100)
ggplot(plot_sf) + 
  geom_sf(aes(fill = CV)) +
  facet_wrap(~method, nrow = 2) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "CV (%)", na.value="white") +
  theme_minimal()  + my_theme 
ggsave("results/figures/Zambia_adm2_cv_map.png",
       width = 12.5, height = 8)
cv_table <- adm2_est_table |>
  filter(method == "Direct") |>
  mutate(CV = sqrt(var) / median * 100) |>
  dplyr::select(domain, median, CV) |>
  arrange(CV) 

cv_table |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3) |>
  writeLines(paste0("results/figures/zmb_hiv_adm2_cv_table.tex"))
#### **** Interval lengths **** ####
ggplot(plot_sf) + 
  geom_sf(aes(fill = upper - lower)) +
  facet_wrap(~method) +
  scale_fill_whitebox_c(palette = "bl_yl_rd") +
  labs(fill = "CI length") + 
  theme_minimal()  + my_theme
ggsave("results/figures/Zambia_adm2_int_length_map.png",
       width = 12.5, height = 8)
#### **** Credible intervals **** ####
direct_est <- adm2_est_table |> dplyr::filter(method == "Direct") |> dplyr::filter(!is.na(mean))
sorted_levels <- direct_est$domain[order(direct_est$mean, decreasing = T)] 
sorted_levels <- c(head(sorted_levels, n = 5), tail(sorted_levels, n = 5))
group_label <- c(rep("Highest direct estimates", 5), rep("Lowest direct estimates", 5))
# top 5 and bottom 5 
int_table <- adm2_est_table |>
  filter(domain %in% sorted_levels) |>
  mutate(group = group_label[match(domain, sorted_levels)])
int_table$domain <- factor(int_table$domain, 
                           levels = sorted_levels)
int_table <- int_table |> dplyr::arrange(domain, method)
ggplot2::ggplot(int_table, ggplot2::aes(x = domain, y = median, color = method)) +
  ggplot2::facet_wrap(~group, scales="free_x", ncol = 1) + 
  ggplot2::geom_point(position = ggplot2::position_dodge(width = .6)) + 
  ggplot2::geom_linerange(ggplot2::aes(x = domain, ymin = lower, ymax = upper), 
                          position = ggplot2::position_dodge(width = .6)) + 
  ggplot2::scale_color_discrete(name = "Method") + 
  ggplot2::ylab("Small Area Estimate") + ggplot2::xlab("Domain") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1),
                 legend.position="bottom",
                 legend.key.width = unit(.3,"cm"),
                 strip.text = element_text(size=14),
                 legend.text = element_text(size=9.5)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE))

ggsave("results/figures/Zambia_adm2_interval_plots.png",
       width = 6, height = 6)
#### **** Ranking plots **** ####
#### **** Estimate tables **** ####

out_table <- adm2_est_table |>
  filter(complete.cases(median)) |>
  mutate(est = paste0(round(median, 2), " (", round(lower, 2),", ", round(upper, 2), ")")) |>
  dplyr::select(est, median, method, domain) |>
  pivot_wider(values_from = c("median", "est"),
              names_from = method,
              names_glue = "{method}_{.value}") |>
  arrange(Direct_median) |>
  dplyr::select(domain, contains("est")) 
names(out_table) <- 
  c("District",
    stringr::str_replace(names(out_table)[-1], "_est", ""))
out_table |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3) |>
  writeLines(paste0("results/figures/zmb_hiv_adm2_est_table.tex"))


#### **** Scatter Plots **** ####


ggplot(adm2_est_table, aes(x = direct.est, y = mean, color = method)) +
  geom_point(size = 3) +
  geom_abline(slope = 1) +
  xlab("Direct Estimate") +
  ylab("Model Estimate") +
  theme_minimal()  + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position="bottom",
    panel.grid.minor = element_blank(),
    legend.key.size = unit(.5, "cm"),
    legend.text = element_text(size=14),
    legend.key.width = unit(.5,"cm"))
ggsave("results/figures/Zambia_adm2_scatter.png",
       width = 10, height = 8)

ggplot(adm2_est_table, aes(x = direct.est, y = mean, color = method)) +
  geom_point(size = 3, show.legend = F) +
  geom_abline(slope = 1) +
  geom_linerange(aes(ymin = lower, ymax = upper), show.legend = F) + 
  facet_wrap(~method) + 
  xlab("Direct Estimate") +
  ylab("Model Estimate") +
  theme_minimal()  +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position="bottom",
    panel.grid.minor = element_blank(),
    legend.key.size = unit(.5, "cm"),
    legend.text = element_text(size=14),
    legend.key.width = unit(.5,"cm"))
ggsave("results/figures/Zambia_adm2_scatter_bars.png",
       width = 10, height = 8)

#### **** Binomial vs. Betabinomial **** ####

bbin_grf_comp <- adm2_est_table |>
  filter(method %in% c("Binomial GRF", "Betabinomial GRF")) |>
  mutate(int_length = upper - lower,
         method = str_replace(method, " GRF", "")) |>
  select(domain, method, int_length) |>
  pivot_wider(names_from = method, values_from = int_length) |>
  mutate(model = "GRF")

bbin_bym2_comp <- adm2_est_table |>
  filter(method %in% c("Binomial BYM2", "Betabinomial BYM2")) |>
  mutate(int_length = upper - lower,
         method = str_replace(method, " BYM2", "")) |>
  mutate(int_length = upper - lower) |>
  select(domain, method, int_length) |>
  pivot_wider(names_from = method, values_from = int_length) |>
  mutate(model = "BYM2")

bbin_comp <- bind_rows(bbin_bym2_comp, bbin_grf_comp)
bbin_comp |> group_by(model) |>
  summarize(mean_int_length_adm2_binomial = mean(Binomial),
            mean_int_length_adm2_betabinomial = mean(Betabinomial))
ggplot(bbin_comp, aes(x = Binomial, y = Betabinomial)) + 
  facet_wrap(~model) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  ggtitle("Interval lengths (Admin-2)") +
  theme_minimal(base_size = 20) 
ggsave("results/figures/Zambia_bin_vs_bbin_int_length_adm2.pdf",
       width = 10, height = 6)


#### **** Admin-1 vs. Admin-2 Variability **** ####
combined_table <- bind_rows(adm1_est_table |> mutate(level = "Admin-1"),
                            adm2_est_table |> mutate(level = "Admin-2"))
ggplot(combined_table, aes(x = median, y = level, color = method))+
  facet_wrap(~method) + 
  geom_boxplot(show.legend = F) +
  theme_minimal()  +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position="bottom",
    panel.grid.minor = element_blank(),
    legend.key.size = unit(.5, "cm"),
    legend.text = element_text(size=14),
    legend.key.width = unit(.5,"cm"))
ggsave("results/figures/Zambia_adm1_vs_adm2_boxplots.png",
       width = 10, height = 8)

ggplot(combined_table, aes(x = median, y = level, color = method))+
  facet_wrap(~method) + 
  geom_violin(show.legend = F) +
  theme_minimal()  +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position="bottom",
    panel.grid.minor = element_blank(),
    legend.key.size = unit(.5, "cm"),
    legend.text = element_text(size=14),
    legend.key.width = unit(.5,"cm"))
ggsave("results/figures/Zambia_adm1_vs_adm2_violins.png",
       width = 10, height = 8)


#### CROSS-VALIDATION ####
holdout_res <- readRDS("results/estimates/holdout_res.rds")
holdout_res |>
  group_by(method) |>
  summarize(rmse = sqrt(mean((median - direct.est) ^ 2)) * 100,
            mae = mean(abs(median - direct.est)) * 100,
            cov = mean(lower < direct.est & upper > direct.est))
holdout_res$method[holdout_res$method == "Area level model: BYM2"] <-
  "Fay-Herriot BYM2"
holdout_res$method[holdout_res$method == "Binomial SPDE LGM"] <-
  "Binomial GRF"
holdout_res$method[holdout_res$method == "Betabinomial SPDE LGM"] <-
  "Betabinomial GRF"
selected_methods <- 
  c("Direct",
    "Fay-Herriot BYM2",
    "Binomial BYM2",
    "Betabinomial BYM2",
    "Binomial GRF",
    "Betabinomial GRF")
holdout_res <- holdout_res |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods)) |>
  mutate(error = direct.est - median) 


# logit(direct.est) is normally distributed with some mean and V_i
# from model, we get modeled logit of risk 

holdout_res |>
  mutate(CV = sqrt(direct.est.var) / direct.est * 100) |>
  filter(CV < 10) |>
  group_by(method) |>
  summarize(rmse = sqrt(mean((median - direct.est) ^ 2)) * 100,
            mae = mean(abs(median - direct.est)) * 100,
            cov = mean(lower < direct.est & upper > direct.est),
            log_cpo = mean(dnorm(SUMMER::logit(direct.est),
                                 mean = SUMMER::logit(median), 
                                 sd = sqrt(direct.est.var / direct.est^2 / (1-direct.est) ^2), log = T)),
            is = mean((upper - lower) + 
                        2 / 0.2 * ifelse(direct.est < lower, lower - direct.est, 0) +
                        2 / 0.2 * ifelse(direct.est > upper, direct.est - upper, 0))) |>
  rename("Method" = method,
         "RMSE (x 100)" = rmse,
         "MAE (x 100)" = mae,
         "80% Coverage" = cov,
         "Mean log(CPO)" = log_cpo,
         "Int. Score" = is)

holdout_res |>
  group_by(method) |>
  summarize(rmse = sqrt(mean((median - direct.est) ^ 2)) * 100,
            mae = mean(abs(median - direct.est)) * 100,
            cov = mean(lower < direct.est & upper > direct.est),
            log_cpo = mean(dnorm(SUMMER::logit(direct.est),
                                 mean = SUMMER::logit(median), 
                                 sd = sqrt(direct.est.var / direct.est^2 / (1-direct.est) ^2), log = T)),
            is = mean((upper - lower) + 
                        2 / 0.2 * ifelse(direct.est < lower, lower - direct.est, 0) +
                        2 / 0.2 * ifelse(direct.est > upper, direct.est - upper, 0))) |>
  rename("Method" = method,
         "RMSE (x 100)" = rmse,
         "MAE (x 100)" = mae,
         "80% Coverage" = cov,
         "Mean log(CPO)" = log_cpo,
         "Int. Score" = is) |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3) |>
  writeLines(paste0("results/figures/Zambia_adm2_loocv_tbl.tex"))


plot_sf <- holdout_res |> 
  group_split(method) |>
  map(function(x) poly_adm2 |>
        left_join(x, by = "domain") |>
        mutate(method = x$method[1])) |>
  bind_rows()


ggplot(plot_sf) + geom_sf(aes(fill = error)) +  facet_wrap(~method) + 
  scale_fill_gradientn(
    limits  = c(-.2, .2),
    colours = c("tomato", "white", "dodgerblue"),
    values  = c(0, .5, 1),
  ) + 
  theme_minimal() +
  my_theme
ggsave("results/figures/Zambia_adm2_loocv_map.png",
       width = 12.5, height = 8)








holdout_res <- readRDS("results/estimates/holdout_res_adm1.rds")
holdout_res$method[holdout_res$method == "Area level model: BYM2"] <- "Fay-Herriot BYM2"
selected_methods <- 
  c("Direct",
    "Fay-Herriot BYM2",
    "Binomial BYM2",
    "Betabinomial BYM2",
    "Binomial GRF",
    "Betabinomial GRF")


holdout_res <- holdout_res |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods)) |>
  mutate(error = direct.est - median)
holdout_res <- holdout_res |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods)) |>
  mutate(error = direct.est - median) 
holdout_res |>
  group_by(method) |>
  summarize(rmse = sqrt(mean((median - direct.est) ^ 2)) * 100,
            mae = mean(abs(median - direct.est)) * 100,
            cov = mean(lower < direct.est & upper > direct.est)) |>
  rename("Method" = method,
         "RMSE (x 100)" = rmse,
         "MAE (x 100)" = mae,
         "50% Coverage" = cov) |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3) |>
  writeLines(paste0("results/figures/Zambia_adm1_loocv_tbl.tex"))

plot_sf <- holdout_res |> 
  group_split(method) |>
  map(function(x) poly_adm1 |>
        left_join(x, by = "domain") |>
        mutate(method = x$method[1])) |>
  bind_rows()


ggplot(plot_sf) + geom_sf(aes(fill = error)) +  facet_wrap(~method) + 
  scale_fill_gradientn(
    limits  = c(-.1, .1),
    colours = c("tomato", "white", "dodgerblue"),
    values  = c(0, .5, 1),
  ) + 
  theme_minimal() +
  my_theme
ggsave("results/figures/Zambia_adm1_loocv_map.png",
       width = 12.5, height = 8)


holdout_res |> 
  filter(method %in% c("Binomial BYM2", "Betabinomial BYM2", "Binomial GRF", "Betabinomial GRF")) |>
  ggplot(aes(x = direct.est, y = median, color = method)) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  geom_abline(slope = 1) + 
  facet_wrap(~method)

#### COVARIATE MODELING COMPARISON ####
adm1_est_table <- read.csv("results/estimates/adm1_est.csv")
adm2_est_table <- read.csv("results/estimates/adm2_est.csv")

selected_methods <- 
  c("Fay-Herriot BYM2",
    "Binomial BYM2",
    "Binomial GRF",
    "Fay-Herriot BYM2: No cov.",
    "Binomial BYM2: No cov.",
    "Binomial GRF: No cov.")
adm1_est_table <- adm1_est_table |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods))
plot_sf <- poly_adm1 |>
  left_join(adm1_est_table, by = "domain") |>
  mutate(method = factor(method, levels = selected_methods)) |>
  mutate(CV = sqrt(var) / median * 100)
ggplot(plot_sf) + 
  geom_sf(aes(fill = median)) +
  facet_wrap(~method) +
  scale_fill_whitebox_c(palette = "viridi", direction = -1, limits = c(0, .33)) +
  theme_minimal()  + my_theme 
ggsave("results/figures/Zambia_adm1_covariate_comp_map.png",
       width = 12.5, height = 8)
plot_sf <- plot_sf |>
  mutate(CV = sqrt(var) / median * 100)
ggplot(plot_sf) + 
  geom_sf(aes(fill = CV)) +
  facet_wrap(~method) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "CV (%)") +
  theme_minimal()  + my_theme 
ggsave("results/figures/Zambia_adm1_covariate_comp_CV_map.png",
       width = 12.5, height = 8)
adm2_est_table <- adm2_est_table |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods))
plot_sf <- poly_adm2 |>
  left_join(adm2_est_table, by = "domain") |>
  mutate(method = factor(method, levels = selected_methods)) |>
  mutate(CV = sqrt(var) / median * 100)
ggplot(plot_sf) + 
  geom_sf(aes(fill = median)) +
  facet_wrap(~method) +
  scale_fill_whitebox_c(palette = "viridi", direction = -1, limits = c(0, .33)) +
  theme_minimal()  + my_theme 
ggsave("results/figures/Zambia_adm2_covariate_comp_map.png",
       width = 12.5, height = 8)
ggplot(plot_sf) + 
  geom_sf(aes(fill = CV)) +
  facet_wrap(~method) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "CV (%)") +
  theme_minimal()  + my_theme 
ggsave("results/figures/Zambia_adm2_covariate_comp_CV_map.png",
       width = 12.5, height = 8)

holdout_res <- readRDS("results/estimates/holdout_res.rds")
holdout_res |>
  group_by(method) |>
  summarize(rmse = sqrt(mean((median - direct.est) ^ 2)) * 100,
            mae = mean(abs(median - direct.est)) * 100,
            cov = mean(lower < direct.est & upper > direct.est))
holdout_res$method[holdout_res$method == "Area level model: BYM2"] <- "Fay-Herriot BYM2"
holdout_res$method[holdout_res$method == "Area level model: BYM2: No cov."] <- "Fay-Herriot BYM2: No cov."
selected_methods <- 
  c("Fay-Herriot BYM2",
    "Binomial BYM2",
    "Binomial GRF",
    "Fay-Herriot BYM2: No cov.",
    "Binomial BYM2: No cov.",
    "Binomial GRF: No cov.")
holdout_res <- holdout_res |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods)) |>
  mutate(error = direct.est - median) 
holdout_res |>
  group_by(method) |>
  summarize(rmse = sqrt(mean((median - direct.est) ^ 2)) * 100,
            mae = mean(abs(median - direct.est)) * 100,
            cov = mean(lower < direct.est & upper > direct.est)) |>
  rename("Method" = method,
         "RMSE (x 100)" = rmse,
         "MAE (x 100)" = mae,
         "80% Coverage" = cov) |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3) |>
  writeLines(paste0("results/figures/Zambia_adm2_covariate_comp_loocv_tbl.tex"))


plot_sf <- holdout_res |> 
  group_split(method) |>
  map(function(x) poly_adm2 |>
        left_join(x, by = "domain") |>
        mutate(method = x$method[1])) |>
  bind_rows()


ggplot(plot_sf) + geom_sf(aes(fill = error)) +  facet_wrap(~method) + 
  scale_fill_gradientn(
    limits  = c(-.2, .2),
    colours = c("tomato", "white", "dodgerblue"),
    values  = c(0, .5, 1),
  ) + 
  theme_minimal() +
  my_theme
ggsave("results/figures/Zambia_adm2_covariate_comp_loocv_map.png",
       width = 12.5, height = 8)

holdout_res <- readRDS("results/estimates/holdout_res_adm1.rds")
holdout_res |>
  group_by(method) |>
  summarize(rmse = sqrt(mean((median - direct.est) ^ 2)) * 100,
            mae = mean(abs(median - direct.est)) * 100,
            cov = mean(lower < direct.est & upper > direct.est))
holdout_res$method[holdout_res$method == "Area level model: BYM2"] <- "Fay-Herriot BYM2"
holdout_res$method[holdout_res$method == "Area level model: BYM2: No cov."] <- "Fay-Herriot BYM2: No cov."
selected_methods <- 
  c("Fay-Herriot BYM2",
    "Binomial BYM2",
    "Binomial GRF",
    "Fay-Herriot BYM2: No cov.",
    "Binomial BYM2: No cov.",
    "Binomial GRF: No cov.")
holdout_res <- holdout_res |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods)) |>
  mutate(error = direct.est - median) 
holdout_res |>
  group_by(method) |>
  summarize(rmse = sqrt(mean((median - direct.est) ^ 2)) * 100,
            mae = mean(abs(median - direct.est)) * 100,
            cov = mean(lower < direct.est & upper > direct.est)) |>
  rename("Method" = method,
         "RMSE (x 100)" = rmse,
         "MAE (x 100)" = mae,
         "50% Coverage" = cov) |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3) |>
  writeLines(paste0("results/figures/Zambia_adm1_covariate_comp_loocv_tbl.tex"))


plot_sf <- holdout_res |> 
  group_split(method) |>
  map(function(x) poly_adm1 |>
        left_join(x, by = "domain") |>
        mutate(method = x$method[1])) |>
  bind_rows()


ggplot(plot_sf) + geom_sf(aes(fill = error)) +  facet_wrap(~method) + 
  scale_fill_gradientn(
    limits  = c(-.2, .2),
    colours = c("tomato", "white", "dodgerblue"),
    values  = c(0, .5, 1),
  ) + 
  theme_minimal() +
  my_theme
ggsave("results/figures/Zambia_adm1_covariate_comp_loocv_map.png",
       width = 12.5, height = 8)


adm1_est_table <- read.csv("results/estimates/adm1_est.csv")
adm2_est_table <- read.csv("results/estimates/adm2_est.csv")

selected_methods <- 
  c("Direct",
    "Fay-Herriot BYM2",
    "Fay-Herriot BYM2: No cov.",
    "Binomial BYM2",
    "Binomial BYM2: No cov.",
    "Binomial GRF",
    "Binomial GRF: No cov.")
adm1_est_table <- adm1_est_table |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods))
direct_est <- adm1_est_table |> dplyr::filter(method == "Direct") |> dplyr::filter(!is.na(mean))
sorted_levels <- direct_est$domain[order(direct_est$mean, decreasing = T)] 
sorted_levels <- c(head(sorted_levels, n = 5), tail(sorted_levels, n = 5))
group_label <- c(rep("Highest direct estimates", 5), rep("Lowest direct estimates", 5))
# top 5 and bottom 5 
int_table <- adm1_est_table |>
  filter(domain %in% sorted_levels) |>
  mutate(group = group_label[match(domain, sorted_levels)])
int_table$domain <- factor(int_table$domain, 
                           levels = sorted_levels)
int_table <- int_table |> dplyr::arrange(domain, method)
ggplot2::ggplot(int_table, ggplot2::aes(x = domain, y = median, color = method)) +
  ggplot2::facet_wrap(~group, scales="free_x", ncol = 1) + 
  ggplot2::geom_point(position = ggplot2::position_dodge(width = .6)) + 
  ggplot2::geom_linerange(ggplot2::aes(x = domain, ymin = lower, ymax = upper), 
                          position = ggplot2::position_dodge(width = .6)) + 
  ggplot2::scale_color_discrete(name = "Method") + 
  ggplot2::ylab("Small Area Estimate") + ggplot2::xlab("Domain") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1),
                 legend.position="bottom",
                 legend.key.width = unit(.3,"cm"),
                 strip.text = element_text(size=14),
                 legend.text = element_text(size=9.5)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("results/figures/Zambia_adm1_covariate_comp_interval_plots.png",
       width = 7.5, height = 6)
adm2_est_table <- adm2_est_table |>
  filter(method %in% selected_methods) |>
  mutate(method = factor(method, levels = selected_methods))
direct_est <- adm2_est_table |> dplyr::filter(method == "Direct") |> dplyr::filter(!is.na(mean))
sorted_levels <- direct_est$domain[order(direct_est$mean, decreasing = T)] 
sorted_levels <- c(head(sorted_levels, n = 5), tail(sorted_levels, n = 5))
group_label <- c(rep("Highest direct estimates", 5), rep("Lowest direct estimates", 5))
# top 5 and bottom 5 
int_table <- adm2_est_table |>
  filter(domain %in% sorted_levels) |>
  mutate(group = group_label[match(domain, sorted_levels)])
int_table$domain <- factor(int_table$domain, 
                           levels = sorted_levels)
int_table <- int_table |> dplyr::arrange(domain, method)
ggplot2::ggplot(int_table, ggplot2::aes(x = domain, y = median, color = method)) +
  ggplot2::facet_wrap(~group, scales="free_x", ncol = 1) + 
  ggplot2::geom_point(position = ggplot2::position_dodge(width = .6)) + 
  ggplot2::geom_linerange(ggplot2::aes(x = domain, ymin = lower, ymax = upper), 
                          position = ggplot2::position_dodge(width = .6)) + 
  ggplot2::scale_color_discrete(name = "Method") + 
  ggplot2::ylab("Small Area Estimate") + ggplot2::xlab("Domain") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1),
                 legend.position="bottom",
                 legend.key.width = unit(.3,"cm"),
                 strip.text = element_text(size=14),
                 legend.text = element_text(size=9.5)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("results/figures/Zambia_adm2_covariate_comp_interval_plots.png",
       width = 7.5, height = 6)
#### MODEL SUMMARIES ####

adm1_alm_res <- readRDS("results/estimates/adm1_alm_res.rds")
adm2_alm_res <- readRDS("results/estimates/adm2_alm_res.rds")

bin_adm1_res <- readRDS("results/estimates/bin_adm1_res.rds")
bin_adm2_res <- readRDS("results/estimates/bin_adm2_res.rds")

bbin_adm1_res <- readRDS("results/estimates/bbin_adm1_res.rds")
bbin_adm2_res <- readRDS("results/estimates/bbin_adm2_res.rds")

bin_LGM_fit <- readRDS("results/estimates/bin_LGM_fit.rds")
bbin_LGM_fit <- readRDS("results/estimates/bbin_LGM_fit.rds")

unit_level_cov_adm1 <- bind_rows(
  list("Binomial BYM2" = bin_adm1_res$bym2.model.fit$summary.fixed,
       "Betabinomial BYM2" =  bbin_adm1_res$betabinomial.bym2.model.fit$summary.fixed,
       "Binomial GRF" = bin_LGM_fit$summary.fixed,
       "Betabinomial GRF" = bbin_LGM_fit$summary.fixed) |>
    map(function(x) tibble::rownames_to_column(x, "Variable")),
  .id = "method"
) |>
  mutate(Estimate = paste0(round(mean, 2), " (", round(sd, 2), ")")) |>
  mutate(Variable = stringr::str_replace(Variable, "intercept", "(Intercept)")) |>
  select(method, Variable, Estimate) |>
  pivot_wider(names_from = "method", values_from = "Estimate") |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3,
                caption = "Fixed effects estimates for unit-level covariate models (Admin-1)") |>
  writeLines(paste0("results/figures/zmb_hiv_adm1_unit_cov.tex"))

area_level_cov_adm1 <- bind_rows(
  list("Fay-Herriot BYM2" = adm1_alm_res$bym2.model.fit$summary.fixed) |>
    map(function(x) tibble::rownames_to_column(x, "Variable")),
  .id = "method"
) |>
  mutate(Estimate = paste0(round(mean, 2), " (", round(sd, 2), ")")) |>
  mutate(Variable = stringr::str_replace(Variable, "intercept", "(Intercept)")) |>
  select(method, Variable, Estimate) |>
  pivot_wider(names_from = "method", values_from = "Estimate") |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3,
                caption = "Fixed effects estimates for area-level covariate models (Admin-1)") |>
  writeLines(paste0("results/figures/zmb_hiv_adm1_area_cov.tex"))

unit_level_cov_adm2 <- bind_rows(
  list("Binomial BYM2" = bin_adm2_res$bym2.model.fit$summary.fixed,
       "Betabinomial BYM2" =  bbin_adm2_res$betabinomial.bym2.model.fit$summary.fixed,
       "Binomial GRF" = bin_LGM_fit$summary.fixed,
       "Betabinomial GRF" = bbin_LGM_fit$summary.fixed) |>
    map(function(x) tibble::rownames_to_column(x, "Variable")),
  .id = "method"
) |>
  mutate(Estimate = paste0(round(mean, 2), " (", round(sd, 2), ")")) |>
  mutate(Variable = stringr::str_replace(Variable, "intercept", "(Intercept)")) |>
  select(method, Variable, Estimate) |>
  pivot_wider(names_from = "method", values_from = "Estimate") |>
  knitr::kable(format = "latex", booktabs = T, linesep = "", digits = 3,
               caption = "Fixed effects estimates for unit-level covariate models (Admin-2)") |>
  writeLines(paste0("results/figures/zmb_hiv_adm2_unit_cov.tex"))

area_level_cov_adm2 <- bind_rows(
  list("Fay-Herriot BYM2" = adm2_alm_res$bym2.model.fit$summary.fixed) |>
    map(function(x) tibble::rownames_to_column(x, "Variable")),
  .id = "method"
) |>
  mutate(Estimate = paste0(round(mean, 2), " (", round(sd, 2), ")")) |>
  mutate(Variable = stringr::str_replace(Variable, "intercept", "(Intercept)")) |>
  select(method, Variable, Estimate) |>
  pivot_wider(names_from = "method", values_from = "Estimate") |>
  knitr::kable( format = "latex", booktabs = T, linesep = "", digits = 3,
                caption = "Fixed effects estimates for area-level covariate models (Admin-2)") |>
  writeLines(paste0("results/figures/zmb_hiv_adm2_area_cov.tex"))
hyperpar_list <- list(
  "Fay-Herriot BYM2 Admin-1" = adm1_alm_res$bym2.model.fit$summary.hyperpar,
  "Fay-Herriot BYM2 Admin-2" = adm2_alm_res$bym2.model.fit$summary.hyperpar,
  "Binomial BYM2 Admin-1" = bin_adm1_res$bym2.model.fit$summary.hyperpar,
  "Betabinomial BYM2 Admin-1" = bbin_adm1_res$betabinomial.bym2.model.fit$summary.hyperpar,
  "Binomial BYM2 Admin-2" = bin_adm2_res$bym2.model.fit$summary.hyperpar,
  "Betabinomial BYM2 Admin-2" = bbin_adm2_res$betabinomial.bym2.model.fit$summary.hyperpar,
  "Binomial GRF" = bin_LGM_fit$summary.hyperpar,
  "Betabinomial GRF" = bbin_LGM_fit$summary.hyperpar
)


f <- file("results/figures/zmb_hiv_hyperpar.tex", open = "wt")
for (i in 1:length(hyperpar_list)) {
  hyperpar_list[[i]] |> tibble::rownames_to_column("Variable") |>
    mutate(Estimate = paste0(round(mean, 2), " (", round(sd, 2), ")")) |>
    select(Variable, Estimate) |>
    knitr::kable(format = "latex", booktabs = T, linesep = "", digits = 3,
                 caption = paste0("Hyperparameter estimates for ", names(hyperpar_list)[i], " model")) |>
    write(f, append=TRUE)
}
close(f)

#### EXCESS BINOMIAL VARIATION ####
median_cl_size <- svy_dat |> 
  st_set_geometry(NULL) |>
  group_by(cluster) |>
  summarize(n = n()) |>
  summarize(n = median(n)) |>
  pull(n)
rho <- bbin_adm1_res$betabinomial.bym2.model.fit$summary.hyperpar[1, 3:5]
1 + (median_cl_size - 1) * rho

rho <- bbin_adm2_res$betabinomial.bym2.model.fit$summary.hyperpar[1, 3:5]
1 + (median_cl_size - 1) * rho

rho <- bbin_LGM_fit$summary.hyperpar[1, 3:5]
1 + (median_cl_size - 1) * rho


rho <- bbin_nested_adm2_res$betabinomial.bym2.model.fit$summary.hyperpar[1, 3:5]
1 + (median_cl_size - 1) * rho


bin_adm1_res$bym2.model.fit$summary.hyperpar
bin_adm2_res$bym2.model.fit$summary.hyperpar

#### URBAN PROPORTION ####

sample_des <- svydesign(id = ~cluster + hshold,
                        strata = ~stratum, nest = T,
                        weights = ~wt, data = svy_dat)
f <- file("results/figures/zmb_hiv_urban_only_fe.tex", open = "wt")
summary(svyglm(hiv ~ urban, design = sample_des))$coefficients |>
  knitr::kable(format = "latex", booktabs = T, linesep = "", digits = 3,
               caption = paste0("Parameter estimates for urban-only model")) |>
  write(f, append=TRUE)
close(f)
urban_props <- svy_dat |> 
  st_set_geometry(NULL) |> 
  group_by(admin1_name) |>
  summarize(hiv_urban_prop = mean(urban)) |>
  mutate(frame_urban_prop = c(.277, .789, .138, .181, .821, .172, .184, .222, .278, .143))
ggplot(urban_props, aes(x = frame_urban_prop, y = hiv_urban_prop)) +
  geom_point() +
  geom_abline(slope = 1) +
  ylab("Urban proportion in sample") + 
  xlab("Urban proportion in frame") + 
  theme_minimal()
ggsave("results/figures/Zambia_adm1_urban_oversampling.png",
       width = 7.5, height = 6)


#### METHOD X METHOD scatter ####
methods_grid <- expand.grid(
  method1 = selected_methods,
  method2 = selected_methods
)
pdf("results/figures/Zambia_adm1_scatter_matrix.pdf",
    width = 10, height = 10)
par(mfrow = c(6, 6), mar = c(4.1, 4.1, 2.1, 1.1))
for (i in 1:nrow(methods_grid)) {
  m1 = methods_grid$method1[i]
  m2 = methods_grid$method2[i]
  plot(adm1_est_table |> filter(method == m1) |> pull(median),
       adm1_est_table |> filter(method == m2) |> pull(median),
       xlab = m1, ylab = m2, pch = 16,
       xlim = c(0, max(adm1_est_table$median)),
       ylim = c(0, max(adm1_est_table$median)))
  abline(0, 1)
}
dev.off()
methods_grid <- expand.grid(
  method1 = selected_methods,
  method2 = selected_methods
)
adm2_est_table <- adm2_est_table |> arrange(domain)
pdf("results/figures/Zambia_adm2_scatter_matrix.pdf",
    width = 14, height = 14)
par(mfrow = c(8, 8), mar = c(4.1, 4.1, 2.1, 1.1))
for (i in 1:nrow(methods_grid)) {
  m1 = methods_grid$method1[i]
  m2 = methods_grid$method2[i]
  plot(adm2_est_table |> filter(method == m1) |> pull(median),
       adm2_est_table |> filter(method == m2) |> pull(median),
       xlab = m1, ylab = m2, pch = 16,
       xlim = c(0, max(na.omit(adm2_est_table$median))),
       ylim = c(0, max(na.omit(adm2_est_table$median))))
  abline(0, 1)
}
dev.off()
adm2_est_table <- adm2_est_table |>
  mutate(CV = sqrt(var) / median * 100)
pdf("results/figures/Zambia_adm2_scatter_cv_matrix.pdf",
    width = 14, height = 14)
par(mfrow = c(8, 8), mar = c(4.1, 4.1, 2.1, 1.1))
for (i in 1:nrow(methods_grid)) {
  m1 = methods_grid$method1[i]
  m2 = methods_grid$method2[i]
  plot(adm2_est_table |> filter(method == m1) |> pull(CV),
       adm2_est_table |> filter(method == m2) |> pull(CV),
       xlab = m1, ylab = m2, pch = 16,
       xlim = c(0, max(na.omit(adm2_est_table$CV))),
       ylim = c(0, max(na.omit(adm2_est_table$CV))))
  abline(0, 1)
}
dev.off()


