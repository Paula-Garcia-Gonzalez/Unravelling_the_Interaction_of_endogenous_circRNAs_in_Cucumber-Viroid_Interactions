# Load required libraries
library(ggplot2)
library(dplyr)

# Load data
data <- read.table("filtrados_con_tamaño.tab", header = TRUE, sep = "\t")

# Create bar plot with IDs on X axis and Length on Y axis
plot <- ggplot(data, aes(x = ID, y = Length)) +
  geom_bar(stat = "identity", fill = "blue", color = "blue") +
  theme(axis.text.x = element_blank(),   # Remove X axis labels
        axis.ticks.x = element_blank(),  # Remove X axis ticks
        axis.title.x = element_blank(),  # Remove X axis title
        axis.title.y = element_text(size = 14)) +
  labs(y = "Length") +
  theme_minimal()

# Export the plot to a file (e.g., PNG)
ggsave("filtered_circulars_plot.png", plot, width = 10, height = 6, dpi = 300)


# Create plot only for aligned circRNAs
# Highlight differentially expressed ones
# Use dots instead of bars (density plot idea)

library(ggplot2)
library(dplyr)

# Load dataset
data <- read.table("filtrados_con_tamaño.tab", header=TRUE, sep="\t")

# Create line plot
line_plot <- data %>%
  ggplot(aes(x = factor(ID, levels = ID), y = Length, group = 1)) +
  geom_line(color = "#69b3a2", size = 1) +
  geom_point(color = "#69b3a2", size = 0.5) +
  labs(x = "Differentially expressed circRNA", y = "circRNA Length", title = "") +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank()
  ) +
  scale_x_discrete(
    breaks = data$ID[seq(1, nrow(data), length.out = 5)],
    labels = data$ID[seq(1, nrow(data), length.out = 5)]
  )

# Display plot
print(line_plot)

# Export if needed
ggsave("differential_density_plot.png", plot = line_plot, width = 10, height = 6, dpi = 300)


# SUPERIMPOSE 12 SIGNIFICANT IDs

library(ggplot2)
library(dplyr)

# Load main dataset
data <- read.table("filtered_with_size.tab", header=TRUE, sep="\t")

# Load significant IDs dataset
significant_data <- read.table("significant_sizes.tab", header=TRUE, sep="\t")

# Create base plot with all data
base_plot <- data %>%
  ggplot(aes(x = factor(ID, levels = ID), y = Length, group = 1)) +
  geom_line(color = "#69b3a2", size = 1) +
  geom_point(color = "#69b3a2", size = 2) +
  labs(x = "Differentially expressed circRNA", y = "circRNA Length", title = "") +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Superimpose significant IDs
superimposed_plot <- base_plot +
  geom_line(data = significant_data, aes(x = factor(ID, levels = data$ID), y = Length, group = 1), 
            color = "#FF5733", size = 1) +
  geom_point(data = significant_data, aes(x = factor(ID, levels = data$ID), y = Length), 
             color = "#FF5733", size = 3) +
  scale_x_discrete(
    breaks = unique(c(data$ID[seq(1, nrow(data), length.out = 5)], significant_data$ID)),
    labels = unique(c(data$ID[seq(1, nrow(data), length.out = 5)], significant_data$ID))
  )

# Display plot
print(superimposed_plot)

# Export if needed
ggsave("superimposed_plot.png", plot = superimposed_plot, width = 12, height = 6, dpi = 300)


# ALTERNATIVE VERSION

library(ggplot2)
library(dplyr)

# Load main dataset
data <- read.table("filtrados_con_tamaño.tab", header=TRUE, sep="\t")

# Load significant IDs dataset
significant_data <- read.table("tamaño_significativos.tab", header=TRUE, sep="\t")

# Create base area plot
base_plot <- data %>%
  ggplot(aes(x = factor(ID, levels = ID), y = Length, group = 1)) +
  geom_area(fill = "#69b3a2", color = "#69b3a2", alpha = 0.4) +
  geom_line(color = "#69b3a2", size = 1) +
  labs(x = "Differentially expressed circRNA", y = "circRNA Length", title = "") +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Superimpose significant IDs and adjust X labels
superimposed_plot <- base_plot +
  geom_area(data = significant_data, aes(x = factor(ID, levels = data$ID), y = Length), 
            fill = "#FF5733", color = "#FF5733", alpha = 0.6) +
  geom_line(data = significant_data, aes(x = factor(ID, levels = data$ID), y = Length), 
            color = "#FF5733", size = 1) +
  scale_x_discrete(
    breaks = significant_data$ID
  )

# Display plot
print(superimposed_plot)

# Export if needed
ggsave("superimposed_plot_with_labels.png", plot = superimposed_plot, width = 12, height = 6, dpi = 300)


# ALTERNATIVE VERSION 2 (marking on X axis)

library(ggplot2)
library(dplyr)

# Load main dataset
data <- read.table("filtrados_con_tamaño.tab", header=TRUE, sep="\t")

# Load significant IDs dataset
significant_data <- read.table("tamaño_significativos.tab", header=TRUE, sep="\t")

# Create base area plot
base_plot <- data %>%
  ggplot(aes(x = factor(ID, levels = ID), y = Length, group = 1)) +
  geom_area(fill = "#69b3a2", color = "#69b3a2", alpha = 0.4) +
  geom_line(color = "#69b3a2", size = 1) +
  labs(x = "Differentially expressed circRNA", y = "circRNA Length", title = "") +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Superimpose significant IDs and add X axis labels
superimposed_plot <- base_plot +
  geom_area(data = significant_data, aes(x = factor(ID, levels = data$ID), y = Length), 
            fill = "#FF5733", color = "#FF5733", alpha = 0.6) +
  geom_line(data = significant_data, aes(x = factor(ID, levels = data$ID), y = Length), 
            color = "#FF5733", size = 1) +
  scale_x_discrete(
    breaks = significant_data$ID,
    labels = significant_data$ID
  ) +
  geom_text(data = significant_data, aes(x = factor(ID, levels = data$ID), y = Length, label = ID),
            angle = 45, vjust = -0.5, hjust = 1, size = 3) +
  geom_segment(data = significant_data, aes(x = factor(ID, levels = data$ID), xend = factor(ID, levels = data$ID), 
                                            y = 0, yend = Length),
               color = "#FF5733", size = 0.5, linetype = "dashed") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(1, 2, 4, 1, "cm"),
    axis.title.x = element_text(margin = margin(t = 20))
  )

# Display plot
print(superimposed_plot)

# Export if needed
ggsave("superimposed_plot_with_labels.png", plot = superimposed_plot, width = 12, height = 6, dpi = 300)
