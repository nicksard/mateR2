#' @title Plot MCMC Diagnostic Plots
#' @description Generates a 6-panel plot including histograms and trace plots for
#'   Np, SR, Mean Mates, and Log Probability based on MCMC output. These plots
#'   are essential for assessing the sampler's convergence and performance.
#' @param mcmc_output The list returned by `generate_map_table()` (must contain
#'   `$samples`, `$history`, and `$acceptance_rate`).
#' @param config_info The config_info data frame used for the MCMC run.
#' @param target_values The list of target values used (`Np_target`, `sr_target`,
#'   `mean_mates_target`).
#' @param mcmc_params The list of MCMC parameters used (needed for `burn_in`).
#' @param plot_filename Optional character string. If provided, saves the plot
#'   to this PNG file using ggsave. Default NULL.
#' @param plot_width Width of the PNG file in pixels (default 1000).
#' @param plot_height Height of the PNG file in pixels (default 1200).
#' @param plot_res Resolution of the PNG file in dpi (default 120).
#' @return A combined `patchwork` / `ggplot` object.
#' @import ggplot2
#' @importFrom patchwork plot_annotation
#' @export
plot_mcmc_diagnostics <- function(
    mcmc_output,
    config_info,
    target_values,
    mcmc_params,
    plot_filename = NULL,
    plot_width = 1000,
    plot_height = 1200,
    plot_res = 120
) {
  # --- Check Inputs ---
  if (!is.list(mcmc_output) || is.null(mcmc_output$samples) || is.null(mcmc_output$history) || is.null(mcmc_output$acceptance_rate)) {
    stop("Input 'mcmc_output' does not seem to be a valid MCMC results structure.")
  }
  if (!is.data.frame(config_info) || !all(c("Males", "Females") %in% names(config_info))) {
    stop("'config_info' must be a data frame with 'Males' and 'Females' columns.")
  }

  # --- Calculate Stats from Samples ---
  num_samples <- length(mcmc_output$samples)
  samples_df <- data.frame(Np = numeric(0), sr = numeric(0), mm = numeric(0))

  if (num_samples > 0) {
    message("Analyzing ", num_samples, " stored samples for plotting.")
    males_vec <- config_info$Males
    females_vec <- config_info$Females

    stats_list <- lapply(mcmc_output$samples, function(Counts) {
      if(length(Counts) != length(males_vec)) return(data.frame(Np=NA, sr=NA, mm=NA))
      Nm = sum(males_vec * Counts)
      Nf = sum(females_vec * Counts)
      Np = Nm + Nf
      sr = ifelse(Nf > 0, Nm / Nf, Inf)
      Total_Matings = sum(males_vec * females_vec * Counts)
      mm = ifelse(Np > 0, Total_Matings / Np, 0)
      return(data.frame(Np = Np, sr = sr, mm = mm))
    })

    samples_df <- do.call(rbind, stats_list)
    samples_df <- samples_df[is.finite(samples_df$Np) & is.finite(samples_df$sr) & is.finite(samples_df$mm), ]
  } else {
    message("No samples available for plotting distributions.")
  }

  # Extract History
  history_df <- as.data.frame(mcmc_output$history)
  have_samples <- nrow(samples_df) > 0
  have_history <- nrow(history_df) > 0 && all(c("iteration", "log_prob", "Np", "sr", "mean_mates") %in% names(history_df))

  # Base theme for consistency
  base_theme <- theme_minimal() + theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = c(0.98, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.85), color = "grey80", linewidth = 0.3),
    legend.title = element_blank(),
    legend.margin = margin(t = 2, r = 4, b = 2, l = 4)
  )

  empty_plot <- function(title) { ggplot() + theme_void() + ggtitle(paste(title, "(No Data)")) }

  # --- Row 1: Np ---
  if (have_samples) {
    p_hist_np <- ggplot(samples_df, aes(x = Np)) +
      geom_histogram(bins = 30, fill = "lightblue", color = "white") +
      geom_vline(aes(xintercept = target_values$Np_target, color = "Target", linetype = "Target"), linewidth = 1) +
      geom_vline(aes(xintercept = mean(Np), color = "Mean", linetype = "Mean"), linewidth = 1) +
      scale_color_manual(name = "", values = c("Target" = "red", "Mean" = "darkblue")) +
      scale_linetype_manual(name = "", values = c("Target" = "solid", "Mean" = "dashed")) +
      labs(title = "Distribution of Sampled Np", x = "Np", y = "Frequency") +
      base_theme
  } else { p_hist_np <- empty_plot("Np Histogram") }

  if (have_history) {
    p_trace_np <- ggplot(history_df, aes(x = iteration, y = Np)) +
      geom_line(color = "skyblue4") +
      geom_vline(aes(xintercept = mcmc_params$burn_in, color = "End Burn-in", linetype = "End Burn-in"), linewidth = 1) +
      geom_hline(aes(yintercept = target_values$Np_target, color = "Target", linetype = "Target"), linewidth = 1) +
      scale_color_manual(name = "", values = c("Target" = "blue", "End Burn-in" = "red")) +
      scale_linetype_manual(name = "", values = c("Target" = "solid", "End Burn-in" = "dashed")) +
      labs(title = "Trace: Np", x = "Iteration", y = "Np") +
      base_theme
  } else { p_trace_np <- empty_plot("Np Trace") }

  # --- Row 2: SR ---
  if (have_samples) {
    p_hist_sr <- ggplot(samples_df, aes(x = sr)) +
      geom_histogram(bins = 30, fill = "lightgreen", color = "white") +
      geom_vline(aes(xintercept = target_values$sr_target, color = "Target", linetype = "Target"), linewidth = 1) +
      geom_vline(aes(xintercept = mean(sr), color = "Mean", linetype = "Mean"), linewidth = 1) +
      scale_color_manual(name = "", values = c("Target" = "red", "Mean" = "darkgreen")) +
      scale_linetype_manual(name = "", values = c("Target" = "solid", "Mean" = "dashed")) +
      labs(title = "Distribution of Sampled SR", x = "Sex Ratio", y = "Frequency") +
      base_theme
  } else { p_hist_sr <- empty_plot("SR Histogram") }

  if (have_history) {
    p_trace_sr <- ggplot(history_df, aes(x = iteration, y = sr)) +
      geom_line(color = "seagreen4") +
      geom_vline(aes(xintercept = mcmc_params$burn_in, color = "End Burn-in", linetype = "End Burn-in"), linewidth = 1) +
      geom_hline(aes(yintercept = target_values$sr_target, color = "Target", linetype = "Target"), linewidth = 1) +
      scale_color_manual(name = "", values = c("Target" = "blue", "End Burn-in" = "red")) +
      scale_linetype_manual(name = "", values = c("Target" = "solid", "End Burn-in" = "dashed")) +
      labs(title = "Trace: SR", x = "Iteration", y = "Sex Ratio") +
      base_theme
  } else { p_trace_sr <- empty_plot("SR Trace") }

  # --- Row 3: Mean Mates & Log Prob ---
  if (have_samples) {
    p_hist_mm <- ggplot(samples_df, aes(x = mm)) +
      geom_histogram(bins = 30, fill = "lightcoral", color = "white") +
      geom_vline(aes(xintercept = target_values$mean_mates_target, color = "Target", linetype = "Target"), linewidth = 1) +
      geom_vline(aes(xintercept = mean(mm), color = "Mean", linetype = "Mean"), linewidth = 1) +
      scale_color_manual(name = "", values = c("Target" = "red", "Mean" = "darkred")) +
      scale_linetype_manual(name = "", values = c("Target" = "solid", "Mean" = "dashed")) +
      labs(title = "Distribution of Mean Mates", x = "Mean Mates", y = "Frequency") +
      base_theme
  } else { p_hist_mm <- empty_plot("Mean Mates Histogram") }

  if (have_history) {
    p_trace_lp <- ggplot(history_df, aes(x = iteration, y = log_prob)) +
      geom_line(color = "grey30") +
      geom_vline(aes(xintercept = mcmc_params$burn_in, color = "End Burn-in", linetype = "End Burn-in"), linewidth = 1) +
      scale_color_manual(name = "", values = c("End Burn-in" = "red")) +
      scale_linetype_manual(name = "", values = c("End Burn-in" = "dashed")) +
      labs(title = "Trace: Log Probability", x = "Iteration", y = "Log Prob") +
      base_theme
  } else { p_trace_lp <- empty_plot("Log Prob Trace") }

  # --- Combine Plots using patchwork ---
  max_m <- length(unique(config_info$Males))
  title_text <- paste0(
    "MCMC Analysis (Np=", target_values$Np_target,
    ", SR=", target_values$sr_target,
    ", MM=", target_values$mean_mates_target,
    ", MaxMates=", max_m, ")\n",
    "Overall Acceptance Rate: ", round(mcmc_output$acceptance_rate * 100, 2), "%"
  )

  combined_plot <- (p_hist_np | p_trace_np) /
    (p_hist_sr | p_trace_sr) /
    (p_hist_mm | p_trace_lp) +
    patchwork::plot_annotation(
      title = title_text,
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  # --- Save or Print ---
  if (!is.null(plot_filename)) {
    # Convert pixel dimensions to inches for ggsave based on resolution
    w_in <- plot_width / plot_res
    h_in <- plot_height / plot_res
    message("Saving ggplot to: ", plot_filename)
    ggsave(filename = plot_filename, plot = combined_plot, width = w_in, height = h_in, dpi = plot_res, bg = "white")
  } else {
    # This explicitly forces the plot to render in RStudio!
    suppressWarnings(print(combined_plot))
  }

  # Return invisibly so it doesn't print a bunch of text to the console
  return(invisible(combined_plot))

}
