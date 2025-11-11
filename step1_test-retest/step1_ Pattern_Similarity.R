# ======================================================================
# Test-Retest Scatter & Violin Plots
#   Networks: GWM-HFN and GM-GM
#   Datasets: SLIM (long-term: time1 vs time2), BNU-3 (short-term: time1/2/3)
#   Plots:
#     - Edge-wise average connectivity scatter (lower triangle of group mean)
#     - Subject-wise test-retest correlation (violin + half points + box)
# ======================================================================

# ------------------------- Data Import -------------------------
# (Ensure the working directory is set appropriately before running)
library(R.matlab)

# ----- BNU-3 (Short-term, three sessions) -----
# Data path example: E:/Neuroimage/MyProject/GMWM_Network/Data/BNU3
time1_GMWMGM <- readMat("../../Data/BNU3/GWM_HFN_time1.mat")
time2_GMWMGM <- readMat("../../Data/BNU3/GWM_HFN_time2.mat")
time3_GMWMGM <- readMat("../../Data/BNU3/GWM_HFN_time3.mat")

time1_GWMHFN_short <- time1_GMWMGM$time1.GWM.HFN
time2_GWMHFN_short <- time2_GMWMGM$time2.GWM.HFN
time3_GWMHFN_short <- time3_GMWMGM$time3.GWM.HFN

time1_GMGM <- readMat("../../Data/BNU3/GMGM_time1.mat")
time2_GMGM <- readMat("../../Data/BNU3/GMGM_time2.mat")
time3_GMGM <- readMat("../../Data/BNU3/GMGM_time3.mat")

time1_GMGM_short <- time1_GMGM$time1.GMGM
time2_GMGM_short <- time2_GMGM$time2.GMGM
time3_GMGM_short <- time3_GMGM$time3.GMGM

# ----- SLIM (Long-term, two sessions) -----
# Data path example: E:/Neuroimage/MyProject/GMWM_Network/Data/SLIM
time1_GMWMGM <- readMat("../../Data/SLIM/GWM_HFN_time1_common.mat")
time2_GMWMGM <- readMat("../../Data/SLIM/GWM_HFN_time2_common.mat")

time1_GWMHFN_long <- time1_GMWMGM$time1.GWM.HFN.common
time2_GWMHFN_long <- time2_GMWMGM$time2.GWM.HFN.common

time1_GMGM <- readMat("../../Data/SLIM/GMGM_time1_common.mat")
time2_GMGM <- readMat("../../Data/SLIM/GMGM_time2_common.mat")

time1_GMGM_long <- time1_GMGM$time1.GMGM.common
time2_GMGM_long <- time2_GMGM$time2.GMGM.common

# Clean temporary objects
rm(time1_GMWMGM, time2_GMWMGM, time3_GMWMGM, time1_GMGM, time2_GMGM, time3_GMGM)

# ------------------------- Libraries for Plotting -------------------------
library(ggplot2)
library(ggpointdensity)
library(tidyverse)
library(broom)
library(gghalves)

# ------------------------- Helper Functions (Original Style) -------------------------
# Average connectivity matrix across subjects (3D array: node x node x subject)
compute_average_matrix <- function(array_3d) {
  apply(array_3d, c(1, 2), mean)
}

# Extract strict lower triangle (excluding diagonal) as a vector
extract_lower_tri <- function(matrix_2d) {
  lower_tri_indices <- lower.tri(matrix_2d)
  matrix_2d[lower_tri_indices]
}

# Subject-wise test-retest correlation (lower-triangle edge vectors)
compute_subject_correlations <- function(array_time1, array_time2) {
  n_subjects <- dim(array_time1)[3]
  correlations <- numeric(n_subjects)
  for (i in 1:n_subjects) {
    v1 <- extract_lower_tri(array_time1[,,i])
    v2 <- extract_lower_tri(array_time2[,,i])
    correlations[i] <- cor(v1, v2)
  }
  correlations
}

# =============== GWM-HFN: Long-term Scatter (SLIM) ====================
average_matrix_time1 <- compute_average_matrix(time1_GWMHFN_long)
average_matrix_time2 <- compute_average_matrix(time2_GWMHFN_long)

lower_tri_time1 <- extract_lower_tri(average_matrix_time1)
lower_tri_time2 <- extract_lower_tri(average_matrix_time2)

plot_data <- data.frame(Time1 = lower_tri_time1,
                        Time2 = lower_tri_time2)

# Correlation (for record)
cor.test(plot_data$Time1, plot_data$Time2) %>% tidy()

plot_long <- ggplot(plot_data, aes(x = Time1, y = Time2)) +
  geom_pointdensity(size = 3) +
  scale_colour_gradient2(name = 'Counts',
                         low = "#4292C6",
                         high = "#F16913",
                         mid = '#FFF7BC',
                         midpoint = 80,
                         limits = c(1, 180)) +
  geom_abline(intercept = 0, slope = 1, size = 1.5, color = 'red') +
  labs(
    title = "Average GWM-HFN Across Participants",
    subtitle = "Long-term Test-Retest Reliability (SLIM)",
    x = "Average Connection Strengths (Time1)",
    y = "Average Connection Strengths (Time2)",
    color = "Density"
  ) +
  annotate("text",
           x = range(plot_data$Time1, plot_data$Time2)[1] + 0.01 * diff(range(plot_data$Time1, plot_data$Time2)),
           y = range(plot_data$Time1, plot_data$Time2)[2] - 0.01 * diff(range(plot_data$Time1, plot_data$Time2)),
           label = paste0("r = ", round(cor(plot_data$Time1, plot_data$Time2), 3), ", p < 0.001"),
           size = 10, hjust = 0, vjust = 1, fontface = "bold") +
  coord_cartesian(xlim = range(plot_data$Time1, plot_data$Time2),
                  ylim = range(plot_data$Time1, plot_data$Time2)) +
  theme_bw() +
  theme(legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = 'black'),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text  = element_text(size = 14),
        axis.text.x  = element_text(size = 20, face = 'bold'),
        axis.text.y  = element_text(size = 20, face = 'bold'),
        axis.title.x = element_text(size = 22, face = 'bold'),
        axis.title.y = element_text(size = 22, face = 'bold'),
        plot.title   = element_text(size = 25, face = 'bold'),
        plot.subtitle= element_text(size = 23))

print(plot_long)
ggsave("../visualization/Figure3_Test-Retest/Scatter_Plot_of_Average_GWM-HFN_Connectivity_Long-term.png",
       plot_long, width = 8, height = 8, dpi = 300)

# =============== GWM-HFN: Short-term Scatter (BNU-3) ==================

# ---- Time1 vs Time2 ----
average_matrix_time1 <- compute_average_matrix(time1_GWMHFN_short)
average_matrix_time2 <- compute_average_matrix(time2_GWMHFN_short)
lower_tri_time1 <- extract_lower_tri(average_matrix_time1)
lower_tri_time2 <- extract_lower_tri(average_matrix_time2)
plot_data <- data.frame(Time1 = lower_tri_time1, Time2 = lower_tri_time2)

plot_short12 <- ggplot(plot_data, aes(x = Time1, y = Time2)) +
  geom_pointdensity(size = 3) +
  scale_colour_gradient2(name = 'Counts',
                         low = "#4292C6",
                         high = "#F16913",
                         mid = '#FFF7BC',
                         midpoint = 80,
                         limits = c(1, 180)) +
  geom_abline(intercept = 0, slope = 1, size = 1.5, color = 'red') +
  labs(
    title = "Average GWM-HFN Across Participants",
    subtitle = "Short-term Test-Retest Reliability (BNU-3)",
    x = "Average Connection Strengths (Time1)",
    y = "Average Connection Strengths (Time2)",
    color = "Density"
  ) +
  annotate("text",
           x = range(plot_data$Time1, plot_data$Time2)[1] + 0.01 * diff(range(plot_data$Time1, plot_data$Time2)),
           y = range(plot_data$Time1, plot_data$Time2)[2] - 0.01 * diff(range(plot_data$Time1, plot_data$Time2)),
           label = paste0("r = ", round(cor(plot_data$Time1, plot_data$Time2), 3), ", p < 0.001"),
           size = 10, hjust = 0, vjust = 1, fontface = "bold") +
  coord_cartesian(xlim = range(plot_data$Time1, plot_data$Time2),
                  ylim = range(plot_data$Time1, plot_data$Time2)) +
  theme_bw() +
  theme(legend.position = c(0.99,0.01),
        legend.justification = c(1,0),
        legend.background = element_rect(colour = 'black'),
        legend.title = element_text(size = 14,face = 'bold'),
        legend.text  = element_text(size = 14),
        axis.text.x  = element_text(size = 20,face = 'bold'),
        axis.text.y  = element_text(size = 20,face = 'bold'),
        axis.title.x = element_text(size = 22,face = 'bold'),
        axis.title.y = element_text(size = 22,face = 'bold'),
        plot.title   = element_text(size = 25,face = 'bold'),
        plot.subtitle= element_text(size = 23))

print(plot_short12)
ggsave("../visualization/Figure3_Test-Retest/Scatter_Plot_of_Average_GWM-HFN_Connectivity_Short-term12.png",
       plot_short12, width = 8, height = 8, dpi = 300)

# ---- Time1 vs Time3 ----
average_matrix_time1 <- compute_average_matrix(time1_GWMHFN_short)
average_matrix_time3 <- compute_average_matrix(time3_GWMHFN_short)
lower_tri_time1 <- extract_lower_tri(average_matrix_time1)
lower_tri_time3 <- extract_lower_tri(average_matrix_time3)
plot_data <- data.frame(Time1 = lower_tri_time1, Time2 = lower_tri_time3)

plot_short13 <- ggplot(plot_data, aes(x = Time1, y = Time2)) +
  geom_pointdensity(size = 3) +
  scale_colour_gradient2(name = 'Counts',
                         low = "#4292C6",
                         high = "#F16913",
                         mid = '#FFF7BC',
                         midpoint = 80,
                         limits = c(1, 180)) +
  geom_abline(intercept = 0, slope = 1, size = 1.5, color = 'red') +
  labs(
    title = "Average GWM-HFN Across Participants",
    subtitle = "Short-term Test-Retest Reliability (BNU-3)",
    x = "Average Connection Strengths (Time1)",
    y = "Average Connection Strengths (Time3)",
    color = "Density"
  ) +
  annotate("text",
           x = range(plot_data$Time1, plot_data$Time2)[1] + 0.01 * diff(range(plot_data$Time1, plot_data$Time2)),
           y = range(plot_data$Time1, plot_data$Time2)[2] - 0.01 * diff(range(plot_data$Time1, plot_data$Time2)),
           label = paste0("r = ", round(cor(plot_data$Time1, plot_data$Time2), 3), ", p < 0.001"),
           size = 10, hjust = 0, vjust = 1, fontface = "bold") +
  coord_cartesian(xlim = range(plot_data$Time1, plot_data$Time2),
                  ylim = range(plot_data$Time1, plot_data$Time2)) +
  theme_bw() +
  theme(legend.position = c(0.99,0.01),
        legend.justification = c(1,0),
        legend.background = element_rect(colour = 'black'),
        legend.title = element_text(size = 14,face = 'bold'),
        legend.text  = element_text(size = 14),
        axis.text.x  = element_text(size = 20,face = 'bold'),
        axis.text.y  = element_text(size = 20,face = 'bold'),
        axis.title.x = element_text(size = 22,face = 'bold'),
        axis.title.y = element_text(size = 22,face = 'bold'),
        plot.title   = element_text(size = 25,face = 'bold'),
        plot.subtitle= element_text(size = 23))

print(plot_short13)
ggsave("../visualization/Figure3_Test-Retest/Scatter_Plot_of_Average_GWM-HFN_Connectivity_Short-term13.png",
       plot_short13, width = 8, height = 8, dpi = 300)

# ---- Time2 vs Time3 ----
average_matrix_time2 <- compute_average_matrix(time2_GWMHFN_short)
average_matrix_time3 <- compute_average_matrix(time3_GWMHFN_short)
lower_tri_time2 <- extract_lower_tri(average_matrix_time2)
lower_tri_time3 <- extract_lower_tri(average_matrix_time3)
plot_data <- data.frame(Time1 = lower_tri_time2, Time2 = lower_tri_time3)

plot_short23 <- ggplot(plot_data, aes(x = Time1, y = Time2)) +
  geom_pointdensity(size = 3) +
  scale_colour_gradient2(name = 'Counts',
                         low = "#4292C6",
                         high = "#F16913",
                         mid = '#FFF7BC',
                         midpoint = 105,
                         limits = c(1, 210)) +
  geom_abline(intercept = 0, slope = 1, size = 1.5, color = 'red') +
  labs(
    title = "Average GWM-HFN Across Participants",
    subtitle = "Short-term Test-Retest Reliability (BNU-3)",
    x = "Average Connection Strengths (Time2)",
    y = "Average Connection Strengths (Time3)",
    color = "Density"
  ) +
  annotate("text",
           x = range(plot_data$Time1, plot_data$Time2)[1] + 0.01 * diff(range(plot_data$Time1, plot_data$Time2)),
           y = range(plot_data$Time1, plot_data$Time2)[2] - 0.01 * diff(range(plot_data$Time1, plot_data$Time2)),
           label = paste0("r = ", round(cor(plot_data$Time1, plot_data$Time2), 3), ", p < 0.001"),
           size = 10, hjust = 0, vjust = 1, fontface = "bold") +
  coord_cartesian(xlim = range(plot_data$Time1, plot_data$Time2),
                  ylim = range(plot_data$Time1, plot_data$Time2)) +
  theme_bw() +
  theme(legend.position = c(0.99,0.01),
        legend.justification = c(1,0),
        legend.background = element_rect(colour = 'black'),
        legend.title = element_text(size = 14,face = 'bold'),
        legend.text  = element_text(size = 14),
        axis.text.x  = element_text(size = 20,face = 'bold'),
        axis.text.y  = element_text(size = 20,face = 'bold'),
        axis.title.x = element_text(size = 22,face = 'bold'),
        axis.title.y = element_text(size = 22,face = 'bold'),
        plot.title   = element_text(size = 25,face = 'bold'),
        plot.subtitle= element_text(size = 23))

print(plot_short23)
ggsave("../visualization/Figure3_Test-Retest/Scatter_Plot_of_Average_GWM-HFN_Connectivity_Short-term23.png",
       plot_short23, width = 8, height = 8, dpi = 300)

# =============== GWM-HFN: Violin (Subject-wise Correlations) ==========

# ---- Long-term (SLIM) ----
subject_correlations <- compute_subject_correlations(time1_GWMHFN_long, time2_GWMHFN_long)
violin_data <- data.frame(Subject = 1:length(subject_correlations),
                          Correlation = subject_correlations)

violin_plot_long <- ggplot(violin_data, aes(x = 1, y = Correlation)) +
  geom_half_violin(side = "l", fill = "#0073C2FF", color = "#00539CFF", alpha = 0.8) +
  geom_half_point(side = "r", alpha = 0.7, color = "#EFC000FF") +
  geom_boxplot(width = 0.1, fill = "#F8F8F8", color = "#4D4D4D",
               alpha = 0.9, outlier.color = NA) +
  labs(title = "SLIM",
       x = NULL,
       y = "Correlation Coefficients") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 22, face = "bold", hjust = 0))

print(violin_plot_long)
ggsave("../visualization/Figure3_Test-Retest/Violin_Plot_of_Average_GWM-HFN_Connectivity_Long-term.png",
       violin_plot_long, width = 2, height = 4, dpi = 300)

# ---- Short-term Time1 vs Time2 ----
subject_correlations <- compute_subject_correlations(time1_GWMHFN_short, time2_GWMHFN_short)
violin_data <- data.frame(Subject = 1:length(subject_correlations),
                          Correlation = subject_correlations)

violin_plot_short <- ggplot(violin_data, aes(x = 1, y = Correlation)) +
  geom_half_violin(side = "l", fill = "#0073C2FF", color = "#00539CFF", alpha = 0.8) +
  geom_half_point(side = "r", alpha = 0.7, color = "#EFC000FF") +
  geom_boxplot(width = 0.1, fill = "#F8F8F8", color = "#4D4D4D",
               alpha = 0.9, outlier.color = NA) +
  labs(title = "BNU-3",
       x = NULL,
       y = "Correlation Coefficients") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 22, face = "bold", hjust = 0))

print(violin_plot_short)
ggsave("../visualization/Figure3_Test-Retest/Violin_Plot_of_Average_GWM-HFN_Connectivity_Short-term12.png",
       violin_plot_short, width = 2, height = 4, dpi = 300)

# ---- Short-term Time1 vs Time3 ----
subject_correlations <- compute_subject_correlations(time1_GWMHFN_short, time3_GWMHFN_short)
violin_data <- data.frame(Subject = 1:length(subject_correlations),
                          Correlation = subject_correlations)

violin_plot_short <- ggplot(violin_data, aes(x = 1, y = Correlation)) +
  geom_half_violin(side = "l", fill = "#0073C2FF", color = "#00539CFF", alpha = 0.8) +
  geom_half_point(side = "r", alpha = 0.7, color = "#EFC000FF") +
  geom_boxplot(width = 0.1, fill = "#F8F8F8", color = "#4D4D4D",
               alpha = 0.9, outlier.color = NA) +
  labs(title = "BNU-3",
       x = NULL,
       y = "Correlation Coefficients") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 22, face = "bold", hjust = 0))

print(violin_plot_short)
ggsave("../visualization/Figure3_Test-Retest/Violin_Plot_of_Average_GWM-HFN_Connectivity_Short-term13.png",
       violin_plot_short, width = 2, height = 4, dpi = 300)

# ---- Short-term Time2 vs Time3 ----
subject_correlations <- compute_subject_correlations(time2_GWMHFN_short, time3_GWMHFN_short)
violin_data <- data.frame(Subject = 1:length(subject_correlations),
                          Correlation = subject_correlations)

violin_plot_short <- ggplot(violin_data, aes(x = 1, y = Correlation)) +
  geom_half_violin(side = "l", fill = "#0073C2FF", color = "#00539CFF", alpha = 0.8) +
  geom_half_point(side = "r", alpha = 0.7, color = "#EFC000FF") +
  geom_boxplot(width = 0.1, fill = "#F8F8F8", color = "#4D4D4D",
               alpha = 0.9, outlier.color = NA) +
  labs(title = "BNU-3",
       x = NULL,
       y = "Correlation Coefficients") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 22, face = "bold", hjust = 0))

print(violin_plot_short)
ggsave("../visualization/Figure3_Test-Retest/Violin_Plot_of_Average_GWM-HFN_Connectivity_Short-term23.png",
       violin_plot_short, width = 2, height = 4, dpi = 300)

# ======================= GM-GM: Long-term Scatter =====================
avg_gmgm_time1 <- compute_average_matrix(time1_GMGM_long)
avg_gmgm_time2 <- compute_average_matrix(time2_GMGM_long)
lt1 <- extract_lower_tri(avg_gmgm_time1)
lt2 <- extract_lower_tri(avg_gmgm_time2)
plot_data_gmgm_long <- data.frame(Time1 = lt1, Time2 = lt2)

plot_long_gmgm <- ggplot(plot_data_gmgm_long, aes(Time1, Time2)) +
  geom_pointdensity(size = 3) +
  scale_colour_gradient2(name = 'Counts',
                         low = "#4292C6",
                         high = "#F16913",
                         mid = '#FFF7BC',
                         midpoint = 1300,
                         limits = c(1, 2600)) +
  geom_abline(intercept = 0, slope = 1, size = 1.5, color = 'red') +
  labs(
    title = "Average GM-GM Across Participants",
    subtitle = "Long-term Test-Retest Reliability (SLIM)",
    x = "Average Connection Strengths (Time1)",
    y = "Average Connection Strengths (Time2)",
    color = "Density"
  ) +
  annotate("text",
           x = range(plot_data_gmgm_long$Time1, plot_data_gmgm_long$Time2)[1] + 0.01 * diff(range(plot_data_gmgm_long$Time1, plot_data_gmgm_long$Time2)),
           y = range(plot_data_gmgm_long$Time1, plot_data_gmgm_long$Time2)[2] - 0.01 * diff(range(plot_data_gmgm_long$Time1, plot_data_gmgm_long$Time2)),
           label = paste0("r = ", round(cor(plot_data_gmgm_long$Time1, plot_data_gmgm_long$Time2), 3), ", p < 0.001"),
           size = 10, hjust = 0, vjust = 1, fontface = "bold") +
  coord_cartesian(xlim = range(plot_data_gmgm_long$Time1, plot_data_gmgm_long$Time2),
                  ylim = range(plot_data_gmgm_long$Time1, plot_data_gmgm_long$Time2)) +
  theme_bw() +
  theme(legend.position = c(0.99,0.01),
        legend.justification = c(1,0),
        legend.background = element_rect(colour = 'black'),
        legend.title = element_text(size = 14,face = 'bold'),
        legend.text  = element_text(size = 14),
        axis.text.x  = element_text(size = 20,face = 'bold'),
        axis.text.y  = element_text(size = 20,face = 'bold'),
        axis.title.x = element_text(size = 22,face = 'bold'),
        axis.title.y = element_text(size = 22,face = 'bold'),
        plot.title   = element_text(size = 25,face = 'bold'),
        plot.subtitle= element_text(size = 23))

print(plot_long_gmgm)
ggsave("../visualization/Figure3_Test-Retest/Scatter_Plot_of_Average_GM-GM_Connectivity_Long-term.png",
       plot_long_gmgm, width = 8, height = 8, dpi = 300)


# ====================== GM-GM: Short-term Scatter =====================

# Prepare average matrices
avg_gmgm_short_1 <- compute_average_matrix(time1_GMGM_short)
avg_gmgm_short_2 <- compute_average_matrix(time2_GMGM_short)
avg_gmgm_short_3 <- compute_average_matrix(time3_GMGM_short)

# ---- Time1 vs Time2 ----
lt1 <- extract_lower_tri(avg_gmgm_short_1)
lt2 <- extract_lower_tri(avg_gmgm_short_2)
plot_data_gmgm <- data.frame(Time1 = lt1, Time2 = lt2)

plot_gmgm_short12 <- ggplot(plot_data_gmgm, aes(Time1, Time2)) +
  geom_pointdensity(size = 3) +
  scale_colour_gradient2(name='Counts',
                         low="#4292C6",
                         high="#F16913",
                         mid='#FFF7BC',
                         midpoint = 1300,
                         limits = c(1, 2600)) +
  geom_abline(intercept=0, slope=1, size=1.5, color='red') +
  labs(
    title = "Average GM-GM Across Participants",
    subtitle = "Short-term Test-Retest Reliability (BNU-3)",
    x = "Average Connection Strengths (Time1)",
    y = "Average Connection Strengths (Time2)",
    color = "Density"
  ) +
  annotate("text",
           x = range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)[1] + 0.01 * diff(range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)),
           y = range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)[2] - 0.01 * diff(range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)),
           label = paste0("r = ", round(cor(plot_data_gmgm$Time1, plot_data_gmgm$Time2), 3), ", p < 0.001"),
           size = 10, hjust=0, vjust=1, fontface="bold") +
  coord_cartesian(xlim = range(plot_data_gmgm$Time1, plot_data_gmgm$Time2),
                  ylim = range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)) +
  theme_bw() +
  theme(legend.position = c(0.99,0.01),
        legend.justification = c(1,0),
        legend.background = element_rect(colour = 'black'),
        legend.title = element_text(size = 14,face = 'bold'),
        legend.text  = element_text(size = 14),
        axis.text.x  = element_text(size = 20,face = 'bold'),
        axis.text.y  = element_text(size = 20,face = 'bold'),
        axis.title.x = element_text(size = 22,face = 'bold'),
        axis.title.y = element_text(size = 22,face = 'bold'),
        plot.title   = element_text(size = 25,face = 'bold'),
        plot.subtitle= element_text(size = 23))

print(plot_gmgm_short12)
ggsave("../visualization/Figure3_Test-Retest/Scatter_Plot_of_Average_GM-GM_Connectivity_Short-term12.png",
       plot_gmgm_short12, width=8, height=8, dpi=300)

# ---- Time1 vs Time3 ----
lt1 <- extract_lower_tri(avg_gmgm_short_1)
lt3 <- extract_lower_tri(avg_gmgm_short_3)
plot_data_gmgm <- data.frame(Time1 = lt1, Time2 = lt3)

plot_gmgm_short13 <- ggplot(plot_data_gmgm, aes(Time1, Time2)) +
  geom_pointdensity(size = 3) +
  scale_colour_gradient2(name='Counts',
                         low="#4292C6",
                         high="#F16913",
                         mid='#FFF7BC',
                         midpoint = 1300,
                         limits = c(1, 2600)) +
  geom_abline(intercept=0, slope=1, size=1.5, color='red') +
  labs(
    title = "Average GM-GM Across Participants",
    subtitle = "Short-term Test-Retest Reliability (BNU-3)",
    x = "Average Connection Strengths (Time1)",
    y = "Average Connection Strengths (Time3)",
    color = "Density"
  ) +
  annotate("text",
           x = range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)[1] + 0.01 * diff(range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)),
           y = range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)[2] - 0.01 * diff(range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)),
           label = paste0("r = ", round(cor(plot_data_gmgm$Time1, plot_data_gmgm$Time2), 3), ", p < 0.001"),
           size = 10, hjust=0, vjust=1, fontface="bold") +
  coord_cartesian(xlim = range(plot_data_gmgm$Time1, plot_data_gmgm$Time2),
                  ylim = range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)) +
  theme_bw() +
  theme(legend.position = c(0.99,0.01),
        legend.justification = c(1,0),
        legend.background = element_rect(colour = 'black'),
        legend.title = element_text(size = 14,face = 'bold'),
        legend.text  = element_text(size = 14),
        axis.text.x  = element_text(size = 20,face = 'bold'),
        axis.text.y  = element_text(size = 20,face = 'bold'),
        axis.title.x = element_text(size = 22,face = 'bold'),
        axis.title.y = element_text(size = 22,face = 'bold'),
        plot.title   = element_text(size = 25,face = 'bold'),
        plot.subtitle= element_text(size = 23))

print(plot_gmgm_short13)
ggsave("../visualization/Figure3_Test-Retest/Scatter_Plot_of_Average_GM-GM_Connectivity_Short-term13.png",
       plot_gmgm_short13, width=8, height=8, dpi=300)

# ---- Time2 vs Time3 ----
lt2 <- extract_lower_tri(avg_gmgm_short_2)
lt3 <- extract_lower_tri(avg_gmgm_short_3)
plot_data_gmgm <- data.frame(Time1 = lt2, Time2 = lt3)

plot_gmgm_short23 <- ggplot(plot_data_gmgm, aes(Time1, Time2)) +
  geom_pointdensity(size = 3) +
  scale_colour_gradient2(name='Counts',
                         low="#4292C6",
                         high="#F16913",
                         mid='#FFF7BC',
                         midpoint = 1300,
                         limits = c(1, 2600)) +
  geom_abline(intercept=0, slope=1, size=1.5, color='red') +
  labs(
    title = "Average GM-GM Across Participants",
    subtitle = "Short-term Test-Retest Reliability (BNU-3)",
    x = "Average Connection Strengths (Time2)",
    y = "Average Connection Strengths (Time3)",
    color = "Density"
  ) +
  annotate("text",
           x = range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)[1] + 0.01 * diff(range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)),
           y = range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)[2] - 0.01 * diff(range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)),
           label = paste0("r = ", round(cor(plot_data_gmgm$Time1, plot_data_gmgm$Time2), 3), ", p < 0.001"),
           size = 10, hjust=0, vjust=1, fontface="bold") +
  coord_cartesian(xlim = range(plot_data_gmgm$Time1, plot_data_gmgm$Time2),
                  ylim = range(plot_data_gmgm$Time1, plot_data_gmgm$Time2)) +
  theme_bw() +
  theme(legend.position = c(0.99,0.01),
        legend.justification = c(1,0),
        legend.background = element_rect(colour = 'black'),
        legend.title = element_text(size = 14,face = 'bold'),
        legend.text  = element_text(size = 14),
        axis.text.x  = element_text(size = 20,face = 'bold'),
        axis.text.y  = element_text(size = 20,face = 'bold'),
        axis.title.x = element_text(size = 22,face = 'bold'),
        axis.title.y = element_text(size = 22,face = 'bold'),
        plot.title   = element_text(size = 25,face = 'bold'),
        plot.subtitle= element_text(size = 23))

print(plot_gmgm_short23)
ggsave("../visualization/Figure3_Test-Retest/Scatter_Plot_of_Average_GM-GM_Connectivity_Short-term23.png",
       plot_gmgm_short23, width=8, height=8, dpi=300)

# ======================= GM-GM: Violin Plots ==========================

# ---- Long-term (SLIM) ----
subject_correlations <- compute_subject_correlations(time1_GMGM_long, time2_GMGM_long)
violin_data_gmgm <- data.frame(Subject = 1:length(subject_correlations),
                               Correlation = subject_correlations)

violin_plot_gmgm_long <- ggplot(violin_data_gmgm, aes(x = 1, y = Correlation)) +
  geom_half_violin(side = "l", fill = "#0073C2FF", color = "#00539CFF", alpha = 0.8) +
  geom_half_point(side = "r", alpha = 0.7, color = "#EFC000FF") +
  geom_boxplot(width = 0.1, fill = "#F8F8F8", color = "#4D4D4D",
               alpha = 0.9, outlier.color = NA) +
  labs(title = "SLIM",
       x = NULL,
       y = "Correlation Coefficients") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        plot.title   = element_text(size = 22, face = "bold", hjust = 0))

print(violin_plot_gmgm_long)
ggsave("../visualization/Figure3_Test-Retest/Violin_Plot_of_Average_GM-GM_Connectivity_Long-term.png",
       violin_plot_gmgm_long, width = 2, height = 4, dpi = 300)

# ---- Short-term Time1 vs Time2 ----
subject_correlations <- compute_subject_correlations(time1_GMGM_short, time2_GMGM_short)
violin_data_gmgm <- data.frame(Subject = 1:length(subject_correlations),
                               Correlation = subject_correlations)

violin_plot_gmgm_short12 <- ggplot(violin_data_gmgm, aes(x = 1, y = Correlation)) +
  geom_half_violin(side="l", fill="#0073C2FF", color="#00539CFF", alpha=0.8) +
  geom_half_point(side="r", alpha=0.7, color="#EFC000FF") +
  geom_boxplot(width=0.1, fill="#F8F8F8", color="#4D4D4D",
               alpha=0.9, outlier.color=NA) +
  labs(title = "BNU-3",
       x = NULL,
       y = "Correlation Coefficients") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=18, face="bold"),
        axis.title.y=element_text(size=20, face="bold"),
        plot.title  =element_text(size=22, face="bold", hjust=0))
print(violin_plot_gmgm_short12)
ggsave("../visualization/Figure3_Test-Retest/Violin_Plot_of_Average_GM-GM_Connectivity_Short-term12.png",
       violin_plot_gmgm_short12, width=2, height=4, dpi=300)

# ---- Short-term Time1 vs Time3 ----
subject_correlations <- compute_subject_correlations(time1_GMGM_short, time3_GMGM_short)
violin_data_gmgm <- data.frame(Subject = 1:length(subject_correlations),
                               Correlation = subject_correlations)

violin_plot_gmgm_short13 <- ggplot(violin_data_gmgm, aes(x = 1, y = Correlation)) +
  geom_half_violin(side="l", fill="#0073C2FF", color="#00539CFF", alpha=0.8) +
  geom_half_point(side="r", alpha=0.7, color="#EFC000FF") +
  geom_boxplot(width=0.1, fill="#F8F8F8", color="#4D4D4D",
               alpha=0.9, outlier.color=NA) +
  labs(title = "BNU-3",
       x = NULL,
       y = "Correlation Coefficients") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=18, face="bold"),
        axis.title.y=element_text(size=20, face="bold"),
        plot.title  =element_text(size=22, face="bold", hjust=0))
print(violin_plot_gmgm_short13)
ggsave("../visualization/Figure3_Test-Retest/Violin_Plot_of_Average_GM-GM_Connectivity_Short-term13.png",
       violin_plot_gmgm_short13, width=2, height=4, dpi=300)

# ---- Short-term Time2 vs Time3 ----
subject_correlations <- compute_subject_correlations(time2_GMGM_short, time3_GMGM_short)
violin_data_gmgm <- data.frame(Subject = 1:length(subject_correlations),
                               Correlation = subject_correlations)

violin_plot_gmgm_short23 <- ggplot(violin_data_gmgm, aes(x = 1, y = Correlation)) +
  geom_half_violin(side="l", fill="#0073C2FF", color="#00539CFF", alpha=0.8) +
  geom_half_point(side="r", alpha=0.7, color="#EFC000FF") +
  geom_boxplot(width=0.1, fill="#F8F8F8", color="#4D4D4D",
               alpha=0.9, outlier.color=NA) +
  labs(title = "BNU-3",
       x = NULL,
       y = "Correlation Coefficients") +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=18, face="bold"),
        axis.title.y=element_text(size=20, face="bold"),
        plot.title  =element_text(size=22, face="bold", hjust=0))
print(violin_plot_gmgm_short23)
ggsave("../visualization/Figure3_Test-Retest/Violin_Plot_of_Average_GM-GM_Connectivity_Short-term23.png",
       violin_plot_gmgm_short23, width=2, height=4, dpi=300)


# ======================= Summary Table (Save) =========================
# Recompute edge-wise (group-average) and subject-wise correlations
# for all required dataset/timepoint combinations (no reliance on
# overwritten variables above).

# ---- Long-term (SLIM) ----
# GWM-HFN
avg_long_gwmhfn_t1 <- apply(time1_GWMHFN_long, c(1,2), mean)
avg_long_gwmhfn_t2 <- apply(time2_GWMHFN_long, c(1,2), mean)
vec_long_gwmhfn_t1 <- avg_long_gwmhfn_t1[lower.tri(avg_long_gwmhfn_t1)]
vec_long_gwmhfn_t2 <- avg_long_gwmhfn_t2[lower.tri(avg_long_gwmhfn_t2)]
edge_r_long_gwmhfn <- cor(vec_long_gwmhfn_t1, vec_long_gwmhfn_t2)

# Subject-wise
n_subj_long_gwmhfn <- dim(time1_GWMHFN_long)[3]
subcorr_long_gwmhfn <- numeric(n_subj_long_gwmhfn)
for(i in 1:n_subj_long_gwmhfn){
  v1 <- time1_GWMHFN_long[,,i][lower.tri(time1_GWMHFN_long[,,i])]
  v2 <- time2_GWMHFN_long[,,i][lower.tri(time2_GWMHFN_long[,,i])]
  subcorr_long_gwmhfn[i] <- cor(v1,v2)
}

# GM-GM
avg_long_gmgm_t1 <- apply(time1_GMGM_long, c(1,2), mean)
avg_long_gmgm_t2 <- apply(time2_GMGM_long, c(1,2), mean)
vec_long_gmgm_t1 <- avg_long_gmgm_t1[lower.tri(avg_long_gmgm_t1)]
vec_long_gmgm_t2 <- avg_long_gmgm_t2[lower.tri(avg_long_gmgm_t2)]
edge_r_long_gmgm <- cor(vec_long_gmgm_t1, vec_long_gmgm_t2)

n_subj_long_gmgm <- dim(time1_GMGM_long)[3]
subcorr_long_gmgm <- numeric(n_subj_long_gmgm)
for(i in 1:n_subj_long_gmgm){
  v1 <- time1_GMGM_long[,,i][lower.tri(time1_GMGM_long[,,i])]
  v2 <- time2_GMGM_long[,,i][lower.tri(time2_GMGM_long[,,i])]
  subcorr_long_gmgm[i] <- cor(v1,v2)
}

# ---- Short-term (BNU-3) ----
# Precompute average matrices
avg_short_gwmhfn_1 <- apply(time1_GWMHFN_short, c(1,2), mean)
avg_short_gwmhfn_2 <- apply(time2_GWMHFN_short, c(1,2), mean)
avg_short_gwmhfn_3 <- apply(time3_GWMHFN_short, c(1,2), mean)

avg_short_gmgm_1 <- apply(time1_GMGM_short, c(1,2), mean)
avg_short_gmgm_2 <- apply(time2_GMGM_short, c(1,2), mean)
avg_short_gmgm_3 <- apply(time3_GMGM_short, c(1,2), mean)

# Helper inline block for edge-wise correlation
edge_cor_pair <- function(mA, mB){
  vA <- mA[lower.tri(mA)]
  vB <- mB[lower.tri(mB)]
  cor(vA,vB)
}

# GWM-HFN edge-wise
edge_r_short_gwmhfn_12 <- edge_cor_pair(avg_short_gwmhfn_1, avg_short_gwmhfn_2)
edge_r_short_gwmhfn_13 <- edge_cor_pair(avg_short_gwmhfn_1, avg_short_gwmhfn_3)
edge_r_short_gwmhfn_23 <- edge_cor_pair(avg_short_gwmhfn_2, avg_short_gwmhfn_3)

# GM-GM edge-wise
edge_r_short_gmgm_12 <- edge_cor_pair(avg_short_gmgm_1, avg_short_gmgm_2)
edge_r_short_gmgm_13 <- edge_cor_pair(avg_short_gmgm_1, avg_short_gmgm_3)
edge_r_short_gmgm_23 <- edge_cor_pair(avg_short_gmgm_2, avg_short_gmgm_3)

# Subject-wise short-term correlations
get_subcorr <- function(arrA, arrB){
  nS <- dim(arrA)[3]
  out <- numeric(nS)
  for(i in 1:nS){
    v1 <- arrA[,,i][lower.tri(arrA[,,i])]
    v2 <- arrB[,,i][lower.tri(arrB[,,i])]
    out[i] <- cor(v1,v2)
  }
  out
}

subcorr_short_gwmhfn_12 <- get_subcorr(time1_GWMHFN_short, time2_GWMHFN_short)
subcorr_short_gwmhfn_13 <- get_subcorr(time1_GWMHFN_short, time3_GWMHFN_short)
subcorr_short_gwmhfn_23 <- get_subcorr(time2_GWMHFN_short, time3_GWMHFN_short)

subcorr_short_gmgm_12 <- get_subcorr(time1_GMGM_short, time2_GMGM_short)
subcorr_short_gmgm_13 <- get_subcorr(time1_GMGM_short, time3_GMGM_short)
subcorr_short_gmgm_23 <- get_subcorr(time2_GMGM_short, time3_GMGM_short)

# ---- Assemble Summary Data Frame ----
summary_table <- rbind(
  data.frame(Dataset="SLIM", Network="GWM-HFN", TimePair="1-2",
             Edgewise_r=edge_r_long_gwmhfn,
             Subjectwise_r_mean=mean(subcorr_long_gwmhfn),
             Subjectwise_r_sd=sd(subcorr_long_gwmhfn),
             N_subjects=length(subcorr_long_gwmhfn)),
  data.frame(Dataset="SLIM", Network="GM-GM", TimePair="1-2",
             Edgewise_r=edge_r_long_gmgm,
             Subjectwise_r_mean=mean(subcorr_long_gmgm),
             Subjectwise_r_sd=sd(subcorr_long_gmgm),
             N_subjects=length(subcorr_long_gmgm)),
  data.frame(Dataset="BNU-3", Network="GWM-HFN", TimePair="1-2",
             Edgewise_r=edge_r_short_gwmhfn_12,
             Subjectwise_r_mean=mean(subcorr_short_gwmhfn_12),
             Subjectwise_r_sd=sd(subcorr_short_gwmhfn_12),
             N_subjects=length(subcorr_short_gwmhfn_12)),
  data.frame(Dataset="BNU-3", Network="GWM-HFN", TimePair="1-3",
             Edgewise_r=edge_r_short_gwmhfn_13,
             Subjectwise_r_mean=mean(subcorr_short_gwmhfn_13),
             Subjectwise_r_sd=sd(subcorr_short_gwmhfn_13),
             N_subjects=length(subcorr_short_gwmhfn_13)),
  data.frame(Dataset="BNU-3", Network="GWM-HFN", TimePair="2-3",
             Edgewise_r=edge_r_short_gwmhfn_23,
             Subjectwise_r_mean=mean(subcorr_short_gwmhfn_23),
             Subjectwise_r_sd=sd(subcorr_short_gwmhfn_23),
             N_subjects=length(subcorr_short_gwmhfn_23)),
  data.frame(Dataset="BNU-3", Network="GM-GM", TimePair="1-2",
             Edgewise_r=edge_r_short_gmgm_12,
             Subjectwise_r_mean=mean(subcorr_short_gmgm_12),
             Subjectwise_r_sd=sd(subcorr_short_gmgm_12),
             N_subjects=length(subcorr_short_gmgm_12)),
  data.frame(Dataset="BNU-3", Network="GM-GM", TimePair="1-3",
             Edgewise_r=edge_r_short_gmgm_13,
             Subjectwise_r_mean=mean(subcorr_short_gmgm_13),
             Subjectwise_r_sd=sd(subcorr_short_gmgm_13),
             N_subjects=length(subcorr_short_gmgm_13)),
  data.frame(Dataset="BNU-3", Network="GM-GM", TimePair="2-3",
             Edgewise_r=edge_r_short_gmgm_23,
             Subjectwise_r_mean=mean(subcorr_short_gmgm_23),
             Subjectwise_r_sd=sd(subcorr_short_gmgm_23),
             N_subjects=length(subcorr_short_gmgm_23))
)

# Save and print
summary_csv_path <- file.path(out_dir,
                              "Summary_Edge_Subject_Correlations_GWM-HFN_GM-GM.csv")
write.csv(summary_table, summary_csv_path, row.names = FALSE)

cat("\n=== Summary Table (Edge-wise & Subject-wise Correlations) ===\n")
print(summary_table)
cat("\nSaved summary table to:", summary_csv_path, "\n")



