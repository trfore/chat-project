# evoked_ap_analysis.R

library(matools)
library(tidyverse)

# Project Data Processing ----
data_path <- paste(getwd(), "data/cell_attached_grc_evoked_ap", sep = "/")
data_folder <- paste(data_path, "data", sep = "/")
data_parameters <- paste(data_path, "chat_evoked_ap_parameters.csv", sep = "/")

data_collection <- analysis_evoked_ap(
  data_folder_path = data_folder,
  data_parameters = data_parameters
)

# Figure 8C - D: Spike Raster ----
for (name in data_collection$cell_id) {
  df_events <-
    data_collection %>%
    dplyr::filter(cell_id == name) %>%
    dplyr::select(cell_id, experiment_id, data_events) %>%
    tidyr::unnest(data_events)

  plot_raster <-
    plot_spike_raster(
      df = df_events,
      filename = df_events$cell_id[1],
      id = df_events$experiment_id[1],
      time_of_stim = unique(na.omit(df_events$stim_times)),
      condition_names = unique(na.omit(df_events$condition)),
      auto_window = TRUE,
      auto_window_size = 30,
      auto_window_symmetric = FALSE,
      manual_condition_color = c("black", "#00B0F6", "#00BF7D", "#F8766D", "#E76BF3"),
      plot_conditions_only = TRUE
    )

  plot(plot_raster)

  ggplot2::ggsave(paste(df_events$cell_id[1], "spike_rasterplot.pdf", sep = "_"),
    plot = plot_raster,
    width = 4.95, height = 3.17, dpi = 600, units = "in",
    device = "pdf",
    useDingbats = FALSE
  )
}

# Figure 8G - H: Probability ----
df_extracted_histo <-
  data_collection %>%
  select(cell_id, experiment_id, data_histogram) %>%
  unnest(data_histogram)

df_summary <-
  df_extracted_histo %>%
  group_by(condition) %>%
  select(experiment_id, condition, h_bin_centers, h_prob_norm_max, probability) %>%
  unnest(c(h_bin_centers, h_prob_norm_max, probability)) %>%
  group_by(experiment_id, condition, h_bin_centers) %>%
  mutate(
    avg_prob = mean(probability),
    avg_prob_sd = sd(probability),
    avg_prob_sem = sd(probability) / n(),
    avg_prob_norm = mean(h_prob_norm_max),
    avg_prob_norm_sd = sd(h_prob_norm_max),
    avg_prob_norm_sem = sd(h_prob_norm_max) / n(),
    n = n()
  ) %>%
  select(
    experiment_id, condition, h_bin_centers, avg_prob, avg_prob_sd, avg_prob_sem,
    avg_prob_norm, avg_prob_norm_sd, avg_prob_norm_sem, n
  ) %>%
  distinct(avg_prob_norm, .keep_all = TRUE)

# Figure 8G
df_base_inc <-
  df_summary %>%
  filter(experiment_id == "GrC_eAP_inc" & condition %in% c("atropine", "muscarine", "control"))
df_inlay_inc <-
  df_summary %>%
  filter(experiment_id == "GrC_eAP_inc" & condition == "GABAzine")
grc_increase <-
  plot_probability(
    df_base_inc,
    plot_stim_times = c(0, 10, 20),
    plot_value = "avg_prob_norm",
    ymax = 2,
    y_label = "spike probability (norm)",
    df_inset = df_inlay_inc,
    sd_value = "avg_prob_norm_sem",
    manual_condition_color = c("black", "#00B0F6", "#00BF7D", "#F8766D")
  )
plot(grc_increase)

ggsave(
  filename = "grc_prob_inc.pdf",
  grc_increase,
  width = 1.7, height = 1.8,
  dpi = 600, units = "in", device = "pdf"
)

# Figure 8H
df_base_dec <-
  df_summary %>%
  filter(experiment_id == "GrC_eAP_dec" & condition %in% c("control", "muscarine", "atropine"))
df_inlay_dec <-
  df_summary %>%
  filter(experiment_id == "GrC_eAP_dec" & condition == "GABAzine")
grc_decrease <-
  plot_probability(
    df_base_dec,
    plot_stim_times = c(0, 10, 20),
    plot_value = "avg_prob_norm",
    ymax = 1,
    y_label = "spike probability (norm)",
    df_inset = df_inlay_dec,
    sd_value = "avg_prob_norm_sem",
    manual_condition_color = c("black", "#00B0F6", "#00BF7D", "#F8766D")
  )
plot(grc_decrease)

ggsave(
  filename = "grc_prob_dec.pdf",
  grc_decrease,
  width = 1.7, height = 1.8,
  dpi = 600, units = "in", device = "pdf"
)

# Figure 8J: First Spike Latency/Jitter ----
df_events <-
  data_collection %>%
  dplyr::mutate(
    data_collection,
    experiment_id = forcats::fct_recode(experiment_id, "GrC_eAP_inc" = "GrC_eAP_NS_inc")
  ) %>%
  tidyr::unnest_wider("data_events") %>%
  dplyr::filter(experiment_id == "GrC_eAP_inc") %>%
  dplyr::group_by(cell_id, experiment_id) %>%
  dplyr::select(condition, stimulus, jitter, event_index) %>%
  tidyr::unnest(condition, stimulus, jitter, event_index) %>%
  dplyr::filter(condition != "NBQX") %>%
  dplyr::mutate(condition_factor = forcats::fct_relevel(condition, c("control", "muscarine", "atropine", "GABAzine")))

grc_jitter_plot <-
  plot_jitter(
    df = df_events,
    group_value = "condition_factor",
    ymax = 7,
    sd_value = "sd",
    scale_color_manual(values = c("black", "#00B0F6", "#00BF7D", "#F8766D"))
  )

plot(grc_jitter_plot)

ggsave(
  filename = "grc_spike_jitter.pdf",
  grc_jitter_plot,
  width = 1.7, height = 1.8,
  dpi = 600, units = "in", device = "pdf"
)

# Figure 8I: Probability by condition ----
data_collection <- analysis_evoked_ap(
  data_folder_path = data_folder,
  data_parameters = data_parameters,
  histo_bin_size = 10,
  histo_window_size = 50,
  histo_one_event_per_bin = TRUE
)

df_extracted_histo <-
  data_collection %>%
  select(cell_id, experiment_id, data_histogram) %>%
  unnest(data_histogram)

grc_max_probability <-
  plot_max_probability(
    df = df_extracted_histo,
    x_axis_condition = "control",
    y_axis_condition = "muscarine",
    metric = "probability",
    x_label = "control peak spike probability",
    y_label = "muscarine peak spike probability"
  ) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::scale_color_manual(values = c("orange", "black", "grey"))

plot(grc_max_probability)

ggsave(
  filename = "grc_max_probability.pdf",
  grc_max_probability,
  width = 1.7, height = 1.8,
  dpi = 600, units = "in", device = "pdf"
)
