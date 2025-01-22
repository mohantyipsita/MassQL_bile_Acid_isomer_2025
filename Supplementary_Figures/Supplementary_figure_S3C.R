#Collision cross section verse retention time plot for conjugated bile acids
#Load libraries 
library(ggplot2)
library(extrafont)
library(svglite)
# Load fonts
loadfonts(device = "win")

# Load the data table containing collision cross section (CCS) and retention time (RT) 
df <- read.csv("IpsitaBA_bio_ecoded.csv", check.names = FALSE, fileEncoding = "UTF-8", encoding = "Latin-1")

# Display the first few rows of the data frame to verify
head(df)

# Create a color column based on the suffix of the Bile Acid names using pastel colors
df$color <- ifelse(grepl("-4-Aminophenol$", df$Samples), "#26355D",
                   ifelse(grepl("-3-Aminophenol$", df$Samples), "#AF47D2",
                          ifelse(grepl("-2-Aminophenol$", df$Samples), "#FF8F00",
                                 "#337357")))

# Create an alpha column for transparency
df$alpha <- ifelse(df$RT == 9, 0.3, 0.5)

# Define the box dimensions
# Adjust RT to the range 8.95 to 9.10
rt_lower_limit <- 8.985
rt_upper_limit <- 9.015

# Keep the CCS limits as they are
ccs_lower_limit <- 223.7
ccs_upper_limit <- 225.5

# Create scatter plot with CCS in the y-axis and RT in the X-axis
p <- ggplot(df, aes(x = RT, y = avCCS)) +
  geom_point(aes(color = color, alpha = alpha), size = 7) +  # Increase point size for better visibility
  scale_color_identity() +
  scale_alpha_identity() +
  labs(title = "Scatter Plot of CCS vs RT for Bile Acids",
       x = "Retention Time (RT)",
       y = "Collision Cross Section (CCS)") +
  geom_rect(aes(xmin = rt_lower_limit, xmax = rt_upper_limit,
                ymin = ccs_lower_limit, ymax = ccs_upper_limit),
            color = "black", fill = NA, size = 1, linetype = "solid") +
  scale_x_continuous(name = "Retention Time (RT)", limits = c(8.27, 9.1)) +
  scale_y_continuous(name = "Collision Cross Section (CCS)", limits = c(205, 230)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(face = "bold", family = "DejaVu Sans", size = 18),  # Increase title size
    axis.title.x = element_text(face = "bold", family = "DejaVu Sans", size = 16),  # Increase x-axis label size
    axis.title.y = element_text(face = "bold", family = "DejaVu Sans", size = 16),  # Increase y-axis label size
    axis.text.x = element_text(face = "bold", family = "DejaVu Sans", size = 16),   # Increase x-axis text size
    axis.text.y = element_text(face = "bold", family = "DejaVu Sans", size = 16)    # Increase y-axis text size
  ) +
  # Label only the specified samples at RT 9
  geom_text_repel(data = df[df$RT == 9 & df$Samples %in% c("3α12α-2-Aminophenol", "3α12β-2-Aminophenol"), ],
                  aes(label = Samples), size = 4.5,  # Adjust label size for visibility
                  nudge_x = 0.02, nudge_y = 0.2,    # Nudge labels slightly for better positioning
                  box.padding = 2,  # Add padding around the text box
                  point.padding = 2,  # Add padding around the points
                  force = 5,  # Adjust force to move labels to empty spaces
                  min.segment.length = 0.5,  # Set minimum length for connecting lines
                  show.legend = FALSE)               # Hide the legend for text

# Label all other samples
p <- p + geom_text_repel(data = df[df$RT != 9 | !(df$Samples %in% c("3α12α-2-Aminophenol", "3α12β-2-Aminophenol")), ],
                         aes(label = Samples), size = 4.5,  # Adjust label size for visibility
                         nudge_x = 0.02, nudge_y = 0.2,    # Nudge labels slightly for better positioning
                         box.padding = .5,  # Add padding around the text box
                         point.padding = .5,  # Add padding around the points
                         force = 1.5,  # Adjust force to move labels to empty spaces
                         min.segment.length = 0.5,  # Set minimum length for connecting lines
                         show.legend = FALSE)               # Hide the legend for text

# Print the plot before saving
print(p)
# Save the plot svg
ggsave("ScatterplotBA11_Ipsita.svg", plot = p, height = 9, width = 14, device = "svg")