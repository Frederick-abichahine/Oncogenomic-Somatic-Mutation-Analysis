#####################
# BQSR Plot Generator
#####################

# Setting working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading libraries
library(ggplot2)
library(ggthemes)

# Reading data from CSV file
df_control <- read.csv("./data/control.recal.csv", header = TRUE)
df_tumor <- read.csv("./data/tumor.recal.csv", header = TRUE)

# Defining a function to create the BQSR plot
plot_bqsr <- function(df, sample_type, color) {
  ggplot(df, aes(x = AverageReportedQuality, y = EmpiricalQuality)) +
    geom_point(size = 1.5, alpha = 0.6, color = color) +
    geom_smooth(method = "loess", se = FALSE, color = "black", linetype = "dashed", size = 0.7) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("BQSR: Empirical vs Reported Quality Scores (", sample_type, ")", sep = ""),
      x = "Average Reported Quality Score",
      y = "Empirical Quality Score"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    coord_cartesian(
      xlim = c(0, max(df$AverageReportedQuality, na.rm = TRUE)),
      ylim = c(0, max(df$EmpiricalQuality, na.rm = TRUE))
    )
}

# Generating the plots for control and tumor samples
p_control <- plot_bqsr(df_control, sample_type = "Control", color = "#1f77b4")
p_tumor <- plot_bqsr(df_tumor, sample_type = "Tumor", color = "#d62728")

# Saving the plots as images
ggsave("./figures/control_bqsr_plot.png", plot = p_control, width = 7, height = 5, dpi = 300, bg = "white")
ggsave("./figures/tumor_bqsr_plot.png", plot = p_tumor, width = 7, height = 5, dpi = 300, bg = "white")
