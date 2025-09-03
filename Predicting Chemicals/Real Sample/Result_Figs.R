# --- RStudio: Limonene SOA figures (pie + bar + combined) ---

# Load libraries
library(ggplot2)
library(gridExtra)

# Set output directory
setwd("C:/Users/91962/Predicting Chemicals/Real Sample/")

# Example values (edit if you have exact numbers)
total_compounds <- 500
pct_harmful <- 21
pct_safe <- 100 - pct_harmful

# Data for pie chart
pie_data <- data.frame(
  category = c(paste0("Predicted Harmful (~", pct_harmful, "%)"),
               paste0("Predicted Not Harmful (~", pct_safe, "%)")),
  value = c(pct_harmful, pct_safe)
)

# Create pie chart
pie_chart <- ggplot(pie_data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = paste0("Predicted Carcinogenicity of Limonene SOA (~", total_compounds, " compounds)")) +
  theme_void(base_size = 12) +
  theme(legend.title = element_blank())

# Data for bar chart (by molecular size/complexity)
bar_data <- data.frame(
  category = c("Monomers (<C10)", "Dimers", "Oligomers"),
  harmful_pct = c(10, 25, 30)  # adjust if needed
)

# Create bar chart
bar_chart <- ggplot(bar_data, aes(x = category, y = harmful_pct, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(harmful_pct, "%")), vjust = -0.5, size = 4) +
  ylim(0, 40) +
  labs(title = "Predicted Carcinogenicity by Molecular Size (Limonene SOA)",
       y = "Predicted Harmful (%)", x = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# Save pie and bar separately
ggsave("limonene_pie_chart.png", pie_chart, width = 6, height = 6, dpi = 300)
ggsave("limonene_bar_chart.png", bar_chart, width = 7, height = 5, dpi = 300)

# Combine both plots side by side
combined <- grid.arrange(pie_chart, bar_chart, ncol = 2)
ggsave("limonene_pie_bar_combined.png", combined, width = 12, height = 6, dpi = 300)
