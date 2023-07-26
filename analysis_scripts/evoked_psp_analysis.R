# evoked_psp_analysis.R

library(matools)
library(tidyverse)

# Script Set Project Directory ----
here::i_am("analysis_scripts/evoked_psp_analysis.R")

# Project Data Processing ----
data_path <- paste(getwd(), "data/evoked_psp", sep = "/")
data_folder <- paste(data_path, "data", sep = "/")
data_parameters <- paste(data_path, "chat_evoked_psp_parameters.csv", sep = "/")

data_collection <- analysis_evoked_psp(
  data_folder_path = data_folder,
  data_parameters = data_parameters
)

# Project Data Processing for Plots ----
#' @title evoked_data_group
#' @description aggregate amplitude (pA & normalized), failure rate data, ppr (stimulus_num = 2)
#'
#' @param df tibble of aggregated data
#' @param experiment_type str, ex: "GrC_eIPSC"
#' @param stimulus_num int, default 1
evoked_data_group <- function(df, experiment_type, stimulus_num = 1) {
  data <-
    df %>%
    dplyr::filter(experiment_id == experiment_type) %>%
    dplyr::select(cell_id, experiment_id, data_events) %>%
    tidyr::unnest_wider("data_events") %>%
    dplyr::select(cell_id, experiment_id, condition, amplitude, amplitude_normalized, stimulus, event_index, ppr) %>%
    tidyr::unnest(cell_id, experiment_id, condition, amplitude, amplitude_normalized, stimulus, event_index, ppr) %>%
    dplyr::filter(condition != "NA", stimulus == stimulus_num, event_index %in% c(1, NA)) %>%
    dplyr::mutate(
      amplitude = ifelse(is.na(amplitude), 0, amplitude),
      event_index = ifelse(is.na(event_index), 1, event_index),
      amplitude_normalized = ifelse(is.na(amplitude_normalized), 0, amplitude_normalized),
      ppr = ifelse(is.na(ppr), 0, ppr),
      condition_factor = forcats::fct_relevel(condition, c("control", "muscarine", "atropine", "NBQX"))
    ) %>%
    dplyr::group_by(experiment_id, cell_id, condition_factor) %>%
    dplyr::summarise(
      n = n(),
      amp_avg = mean(amplitude),
      amp_sd = sd(amplitude),
      amp_sem = sd(amplitude) / sqrt(n()),
      amp_norm_avg = mean(amplitude_normalized),
      failure_rate = (sum(amplitude == 0) / n()) * 100,
      potency_avg = mean(amplitude[amplitude > 0]),
      potency_normalized_avg = mean(amplitude_normalized[amplitude_normalized > 0]),
      ppr_avg = mean(ppr)
    )
  return(data)
}

#' @title evoked_data_single
#'
#' @param df tibble of aggregated data
#' @param experiment_type str, ex: "GrC_eIPSC"
#' @param recording_name str, file name without extension
#' @param stimulus_num int, default to 1
evoked_data_single <- function(df, experiment_type, recording_name, stimulus_num = 1) {
  df_events <-
    df %>%
    dplyr::filter(experiment_id == experiment_type & cell_id == recording_name) %>%
    dplyr::select(cell_id, experiment_id, data_events) %>%
    tidyr::unnest_wider("data_events") %>%
    dplyr::select(cell_id, experiment_id, sweep, condition, amplitude, stimulus, event_index) %>%
    tidyr::unnest(cell_id, experiment_id, sweep, condition, amplitude, stimulus, event_index) %>%
    dplyr::filter(condition != "NA", stimulus == stimulus_num, event_index %in% c(1, NA)) %>%
    dplyr::mutate(
      sweep = as.numeric(sweep),
      amplitude = ifelse(is.na(amplitude), 0, amplitude),
      event_index = ifelse(is.na(event_index), 1, event_index),
      condition = forcats::fct_relevel(condition, c("control", "muscarine", "atropine", "NBQX"))
    ) %>%
    dplyr::filter(condition %in% c("control", "muscarine", "atropine"))
}

# Figure 5C: Granule Cell Evoked IPSC Scatter Plot ----
experiments <- data_collection %>%
  dplyr::filter(experiment_id == "GrC_eIPSC") %>%
  dplyr::pull(cell_id)

for (cell in experiments) {
  df_events <-
    evoked_data_single(
      data_collection,
      experiment_type = "GrC_eIPSC",
      recording_name = cell,
      stimulus_num = 1
    )

  scatter_plot <-
    plot_scatterplot_amplitude(
      df = df_events,
      filename = cell,
      experiment_id = "GrC Evoke IPSC",
      sweep_duration = 5,
      monochrome = TRUE,
      truncate_x = TRUE,
      ymax = 150,
      y_label = "IPSC (pA)"
    )
  plot(scatter_plot)

  ggplot2::ggsave(paste(cell, "scatter_ipsc_amplitude.pdf", sep = "_"),
    plot = scatter_plot,
    width = 3.71, height = 2.38, dpi = 600, units = "in",
    device = "pdf",
    useDingbats = FALSE
  )
}

# Figure 5D - E: Granule Cell Evoked IPSC EPSC Event Plots ----
df_events_ipsc <- evoked_data_group(data_collection, experiment_type = "GrC_eIPSC")

plot_evoked_ipsc_amp <-
  plot_event_avg(
    df = df_events_ipsc %>% dplyr::filter(condition_factor != "NBQX"),
    x_axis = condition_factor,
    y_axis = amp_avg,
    plot_grouping = "cell_id",
    plot_mean = TRUE,
    sd_value = "sem",
    experiment_name = "GrC eIPSC",
    x_label = "condition",
    y_label = "amplitude (pA)",
    ymax = 100
  )
plot(plot_evoked_ipsc_amp)

ggplot2::ggsave(
  filename = "grc_evoked_ipsc_amp.pdf",
  plot_evoked_ipsc_amp,
  width = 1.2, height = 2.5,
  dpi = 600, units = "in", device = "pdf",
  useDingbats = FALSE
)

# Figure 5D Normalized Amplitude
plot_evoked_ipsc_amp_norm <-
  plot_event_avg(
    df = df_events_ipsc %>% dplyr::filter(condition_factor != "NBQX"),
    x_axis = condition_factor,
    y_axis = amp_norm_avg,
    plot_grouping = "cell_id",
    plot_mean = TRUE,
    sd_value = "sem",
    experiment_name = "GrC eIPSC",
    x_label = "condition",
    y_label = "IPSC (norm)",
    ymax = 2
  )
plot(plot_evoked_ipsc_amp_norm)

ggplot2::ggsave(
  filename = "grc_evoked_ipsc_amp_norm.pdf",
  plot_evoked_ipsc_amp_norm,
  width = 1.2, height = 2.5,
  dpi = 600, units = "in", device = "pdf",
  useDingbats = FALSE
)

# Figure 5E Failure Rate
plot_evoked_ipsc_fr <-
  plot_event_avg(
    df = df_events_ipsc,
    x_axis = condition_factor,
    y_axis = failure_rate,
    plot_grouping = "cell_id",
    plot_mean = TRUE,
    sd_value = "sem",
    experiment_name = "GrC eIPSC",
    x_label = "condition",
    y_label = "failure rate",
    ymax = 100
  )
plot(plot_evoked_ipsc_fr)

ggplot2::ggsave(
  filename = "grc_evoked_ipsc_fr.pdf",
  plot_evoked_ipsc_fr,
  width = 1.2, height = 2.5,
  dpi = 600, units = "in", device = "pdf",
  useDingbats = FALSE
)

# Figure 6C: Granule Cell Evoked EPSC Scatter Plot ----
experiments <- data_collection %>%
  dplyr::filter(experiment_id == "GrC_eEPSC") %>%
  dplyr::pull(cell_id)

for (cell in experiments) {
  df_events <-
    evoked_data_single(
      data_collection,
      experiment_type = "GrC_eEPSC",
      recording_name = cell,
      stimulus_num = 1
    )

  scatter_plot <-
    plot_scatterplot_amplitude(
      df = df_events,
      filename = cell,
      experiment_id = "GrC Evoke EPSC",
      sweep_duration = 5,
      monochrome = TRUE,
      truncate_x = TRUE,
      y_label = "EPSC (pA)"
    )
  plot(scatter_plot)

  ggplot2::ggsave(paste(cell, "scatter_epsc_amplitude.pdf", sep = "_"),
    plot = scatter_plot,
    width = 3.71, height = 2.38, dpi = 600, units = "in",
    device = "pdf",
    useDingbats = FALSE
  )
}

# Figure 6D - F: Granule Cell Evoked EPSC Event Plots ----
df_events_grc_epsc <- evoked_data_group(data_collection, experiment_type = "GrC_eEPSC")

# Figure 6D Normalized Amplitude
plot_evoked_epsc_amp_norm <-
  plot_event_avg(
    df = df_events_grc_epsc %>% dplyr::filter(condition_factor != "NBQX"),
    x_axis = condition_factor,
    y_axis = amp_norm_avg,
    plot_grouping = "cell_id",
    plot_mean = TRUE,
    sd_value = "sem",
    experiment_name = "GrC_eEPSC",
    x_label = "condition",
    y_label = "EPSC (norm)",
    ymax = 2
  )
plot(plot_evoked_epsc_amp_norm)

ggplot2::ggsave(
  filename = "grc_evoked_epsc_amp_norm.pdf",
  plot_evoked_epsc_amp_norm,
  width = 1.2, height = 2.5,
  dpi = 600, units = "in", device = "pdf",
  useDingbats = FALSE
)

# Figure 6E Failure Rate
plot_evoked_epsc_fr <-
  plot_event_avg(
    df = df_events_grc_epsc,
    x_axis = condition_factor,
    y_axis = failure_rate,
    plot_grouping = "cell_id",
    plot_mean = TRUE,
    sd_value = "sem",
    experiment_name = "GrC eEPSC",
    x_label = "condition",
    y_label = "failure rate"
  )
plot(plot_evoked_epsc_fr)

ggplot2::ggsave(
  filename = "grc_evoked_epsc_fr.pdf",
  plot_evoked_epsc_fr,
  width = 1.2, height = 2.5,
  dpi = 600, units = "in", device = "pdf",
  useDingbats = FALSE
)

# Figure 6F Potency
plot_evoked_epsc_potency <-
  plot_event_avg(
    df = df_events_grc_epsc,
    x_axis = condition_factor,
    y_axis = potency_normalized_avg,
    plot_grouping = "cell_id",
    plot_mean = TRUE,
    sd_value = "sem",
    experiment_name = "GrC eEPSC",
    x_label = "condition",
    y_label = "potency (norm)",
    ymax = 2
  )
plot(plot_evoked_epsc_potency)

ggplot2::ggsave(
  filename = "grc_evoked_epsc_potency.pdf",
  plot_evoked_epsc_potency,
  width = 1.2, height = 2.5,
  dpi = 600, units = "in", device = "pdf",
  useDingbats = FALSE
)

# Granule Cell Evoked EPSC PPR (25 Hz) ----
df_events <- evoked_data_group(data_collection, experiment_type = "GrC_eEPSC", stimulus_num = 2)

plot_evoked_grc_ppr <-
  plot_event_avg(
    df = df_events,
    x_axis = condition_factor,
    y_axis = ppr_avg,
    plot_grouping = "cell_id",
    plot_mean = TRUE,
    sd_value = "sem",
    experiment_name = "GrC PPR",
    x_label = "condition",
    y_label = "PPR",
    ymax = 2
  )
plot(plot_evoked_grc_ppr)

ggplot2::ggsave(
  filename = "grc_evoked_epsc_ppr.pdf",
  plot_evoked_grc_ppr,
  width = 1.2, height = 2.5,
  dpi = 600, units = "in", device = "pdf",
  useDingbats = FALSE
)

# Figure 7C: Golgi Cell Evoked EPSC Scatter Plot ----
experiments <- data_collection %>%
  dplyr::filter(experiment_id == "GoC_eEPSC") %>%
  dplyr::pull(cell_id)

for (cell in experiments) {
  df_events <-
    evoked_data_single(
      data_collection,
      experiment_type = "GoC_eEPSC",
      recording_name = cell,
      stimulus_num = 1
    )

  scatter_plot <-
    plot_scatterplot_amplitude(
      df = df_events,
      filename = cell,
      experiment_id = "Golgi Cell Evoke EPSC",
      sweep_duration = 5,
      monochrome = TRUE,
      truncate_x = TRUE,
      y_label = "EPSC (pA)"
    )
  plot(scatter_plot)

  ggplot2::ggsave(paste(cell, "scatter_epsc_amplitude.pdf", sep = "_"),
    plot = scatter_plot,
    width = 3.71, height = 2.38, dpi = 600, units = "in",
    device = "pdf",
    useDingbats = FALSE
  )
}

# Figure 7D - F: Golgi Cell Evoked EPSC Event Plots ----
df_events_goc_epsc <- evoked_data_group(data_collection, experiment_type = "GoC_eEPSC")

# Figure 7D Normalized Amplitude
plot_evoked_goc_epsc_amp_norm <-
  plot_event_avg(
    df = df_events_goc_epsc %>% dplyr::filter(condition_factor != "NBQX"),
    x_axis = condition_factor,
    y_axis = amp_norm_avg,
    plot_grouping = "cell_id",
    plot_mean = TRUE,
    sd_value = "sem",
    experiment_name = "GoC_eEPSC",
    x_label = "condition",
    y_label = "EPSC (norm)",
    ymax = 2
  )
plot(plot_evoked_goc_epsc_amp_norm)

ggplot2::ggsave(
  filename = "goc_evoked_epsc_amp_norm.pdf",
  plot_evoked_goc_epsc_amp_norm,
  width = 1.2, height = 2.5,
  dpi = 600, units = "in", device = "pdf",
  useDingbats = FALSE
)

# Figure 7E Failure Rate
plot_evoked_goc_epsc_fr <-
  plot_event_avg(
    df = df_events_goc_epsc,
    x_axis = condition_factor,
    y_axis = failure_rate,
    plot_grouping = "cell_id",
    plot_mean = TRUE,
    sd_value = "sem",
    experiment_name = "GoC eEPSC",
    x_label = "condition",
    y_label = "failure rate"
  )
plot(plot_evoked_goc_epsc_fr)

ggplot2::ggsave(
  filename = "goc_evoked_epsc_fr.pdf",
  plot_evoked_goc_epsc_fr,
  width = 1.2, height = 2.5,
  dpi = 600, units = "in", device = "pdf",
  useDingbats = FALSE
)

# Figure 7F Potency
plot_evoked_goc_epsc_potency <-
  plot_event_avg(
    df = df_events_goc_epsc,
    x_axis = condition_factor,
    y_axis = potency_normalized_avg,
    plot_grouping = "cell_id",
    plot_mean = TRUE,
    sd_value = "sem",
    experiment_name = "GoC eEPSC",
    x_label = "condition",
    y_label = "potency (norm)",
    ymax = 2
  )
plot(plot_evoked_goc_epsc_potency)

ggplot2::ggsave(
  filename = "goc_evoked_epsc_potency.pdf",
  plot_evoked_goc_epsc_potency,
  width = 1.2, height = 2.5,
  dpi = 600, units = "in", device = "pdf",
  useDingbats = FALSE
)

# Golgi Cell Evoked EPSC PPR (25 Hz) ----
df_events_goc_ppr <- evoked_data_group(data_collection, experiment_type = "GoC_eEPSC", stimulus_num = 2)

plot_evoked_goc_ppr <-
  plot_event_avg(
    df = df_events_goc_ppr,
    x_axis = condition_factor,
    y_axis = ppr_avg,
    plot_grouping = "cell_id",
    plot_mean = TRUE,
    sd_value = "sem",
    experiment_name = "GoC PPR",
    x_label = "condition",
    y_label = "PPR",
    ymax = 2
  )
plot(plot_evoked_goc_ppr)

ggplot2::ggsave(
  filename = "goc_evoked_epsc_ppr.pdf",
  plot_evoked_goc_ppr,
  width = 1.2, height = 2.5,
  dpi = 600, units = "in", device = "pdf",
  useDingbats = FALSE
)

# Golgi Cell Evoked EPSC PPR using GABAb Antagonist
df_events_goc_ppr_gabab <- evoked_data_group(data_collection, experiment_type = "GoC_GABAb_eEPSC", stimulus_num = 2)

plot_evoked_goc_ppr_gabab <-
  plot_event_avg(
    df = df_events_goc_ppr_gabab,
    x_axis = condition_factor,
    y_axis = ppr_avg,
    plot_grouping = "cell_id",
    plot_mean = TRUE,
    sd_value = "sem",
    experiment_name = "GoC PPR (GABAb)",
    x_label = "condition",
    y_label = "PPR",
    ymax = 2
  )
plot(plot_evoked_goc_ppr_gabab)

ggplot2::ggsave(
  filename = "goc_evoked_epsc_ppr_gabab.pdf",
  plot_evoked_goc_ppr_gabab,
  width = 1.2, height = 2.5,
  dpi = 600, units = "in", device = "pdf",
  useDingbats = FALSE
)
