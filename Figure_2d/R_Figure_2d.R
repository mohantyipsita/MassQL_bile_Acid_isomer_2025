######Script for Figure 2d - ITIS diest study in Rheumatoid Arthritis study#########
#GNPS2 final link: https://gnps2.org/status?task=1234bd8eb0f940feba6170da30a07e36

setwd("C:/Users/imoha/OneDrive - University of California, San Diego Health/Bile_Acid_Fragmentation_tree/CASE_guma_RA")

library(stringr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(data.table)
library(ggplot2)
library(tibble)
library(gridExtra)


transpose_and_format_df <- function(df) {
  # Transpose the dataframe
  df_t <- t(df)
  df_t <- as.data.frame(df_t)  # Ensure it's a dataframe
  
  # Use the first row to name the columns
  col_names <- as.character(unlist(df_t[1, ]))
  df_t <- df_t[-1, ]  # Remove the first row used for column names
  colnames(df_t) <- col_names
  
  # Make the row index the first column
  df_t <- cbind(filename = row.names(df_t), df_t)
  
  # Optionally, remove the row names
  row.names(df_t) <- NULL
  
  return(df_t)
}

#read the original BA libhit
df_libhit <- fread("Guma_libhit.tsv") 
df_libhit_req <- df_libhit[,c(2,15,32)]

#read the FBMN table from Audrey dataset
df_feature <- fread("fbmn_quant_Marta_RA.csv")
df_feature_req <- df_feature[,-c(2:3)]

#remove columns that are not needed
df_feature_req_filtered <- df_feature_req %>%
  dplyr::select(1, dplyr::matches("FD|PD"))

#substitute the name
df_feature_req_filtered <- df_feature_req_filtered %>%
  rename_with(~ ifelse(.x == names(df_feature_req_filtered)[1], .x, sub("\\.mzXML Peak area$", "", .x)), everything())

#metadata
df_metadata <- fread("GUMA_RA_metadata.csv")

#subset metadata for timepoint T2 and T3 -- the two timepoints closest to the diet intervention
df_metadata_T2T3 <- df_metadata %>%
  filter(timepoint %in% c("T2", "T3"))

#subset feature table for the DCA-2-Aminophenol -- feature #18901; #22907 and 9173 are added to maintain the structure of the dataframe
df_DCA_2_aminophenol <- df_feature_req_filtered %>%
  dplyr::filter(.[[1]] %in% c("18901", "22097", "9173"))

t_df_DCA_2_aminophenol <- transpose_and_format_df(df_DCA_2_aminophenol)
t_df_DCA_2_aminophenol <- t_df_DCA_2_aminophenol[,-c(2,4)]

#merge with metadata
df_DCA_2_aminophenol_wholegrain <- merge (t_df_DCA_2_aminophenol, df_metadata_T2T3, by.x = "filename", by.y = "UniqueSubjectID", all = FALSE)
df_DCA_2_aminophenol_wholegrain_all <- merge (t_df_DCA_2_aminophenol, df_metadata, by.x = "filename", by.y = "UniqueSubjectID", all = FALSE)

#calculate the whlgrain diet score
df_DCA_2_aminophenol_wholegrain <- df_DCA_2_aminophenol_wholegrain %>%
  mutate(grain_score = whole_grains + 2 * Rye_bread + Corn)

df_DCA_2_aminophenol_wholegrain_all <- df_DCA_2_aminophenol_wholegrain_all %>%
  mutate(grain_score = whole_grains + 2 * Rye_bread + Corn)

#remove the samples for which it is all zeros FD002, FD005, FD018, FD019, FD020, FD026
df_2_AP_filtered <- df_DCA_2_aminophenol_wholegrain_all %>%
  filter(!str_detect(filename, "FD002|FD005|FD007|FD010|FD018|FD019|FD020|FD026|PD002|PD005|PD018|PD019|PD020|PD026"))

ggplot(df_2_AP_filtered, aes(x = Sampletype, y = `log_18901`)) +
  geom_boxplot() +
  labs(title = "Box Plot of 18901 by Sampletype",
       x = "Sample Type",
       y = "18901") +
  theme_minimal()

#######Show for T2 and T3 for non-zero people by box plots or variation in each patient
df_2_AP_filtered_T2T3 <- df_2_AP_filtered %>%
  filter(timepoint %in% c("T2", "T3"))

df_2_AP_filtered_T2T3 <- df_2_AP_filtered_T2T3 %>%
  mutate(person_id = sub("_T[23]", "", filename))

df_2_AP_filtered_T2T3_fecal <- df_2_AP_filtered_T2T3 %>% filter(Sampletype %in% c("Fecal"))

# Reshape the data to long format
df_long_T2T3 <- df_2_AP_filtered_T2T3_fecal %>%
  dplyr::select(timepoint, `18901`, person_id) %>%  # Use dplyr::select to avoid conflicts
  pivot_longer(cols = `18901`, 
               names_to = "measurement", 
               values_to = "peak_area")
df_long_T2T3_noFD003 <- df_long_T2T3 %>%
  filter(!person_id %in% c("FD003"))

df_long_T2T3 <- df_long_T2T3 %>%
  mutate(log_peak_area = log2(`peak_area` + 1))

##########FIANL BOX PLOT + LINE PLOT###############
###color the line plots based on the person_id in increasing order of delta grain score
# Define a color palette for the lines (one color per person)

person_colors <- c(
  "FD003" = "#FFFFBA",  # Light Yellow
  "FD028" = "#FFF176",  # Soft Yellow
  "FD017" = "#FFE066",  # Yellow
  "FD013" = "#FFD54F",  # Golden Yellow
  "FD024" = "#FFB74D",  # Light Orange
  "FD006" = "#FF8A65",  # Orange-Red
  "FD022" = "#FF7043",  # Red-Orange
  "FD004" = "#FF5252",  # Bright Red
  "FD001" = "#FF1744"   # Dark Red
)

# Define custom thickness (linewidth) for each person
person_linewidth <- c(
  "FD003" = 1.5,
  "FD028" = 1.5,
  "FD017" = 1.5,
  "FD013" = 1.5,
  "FD024" = 1.5,
  "FD006" = 1.5,
  "FD022" = 1.5,
  "FD004" = 1.5,
  "FD001" = 1.5
)


df_long_T2T3$person_id <- factor(df_long_T2T3$person_id, levels = c("FD003", "FD028", "FD017", "FD013", "FD024", "FD006", "FD022", "FD004", "FD001"))

plot <- ggplot(df_long_T2T3, aes(x = timepoint, y = log_peak_area)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +  # Create box plots side by side and remove outliers
  geom_jitter(aes(color = timepoint), position = position_dodge(width = 0.8), size = 5, alpha = 1) +  # Add jittered points with jitter dodge for better spacing
  geom_line(aes(group = person_id, color = person_id, linewidth = person_id), alpha = 0.5) +  # Color lines by person_id and adjust linewidth
  scale_color_manual(values = person_colors) +  # Use the predefined color palette for person_id
  scale_linewidth_manual(values = person_linewidth) +  # Use the predefined linewidths
  scale_x_discrete(labels = c("T2" = "Before_diet_intervention", "T3" = "After_diet_intervention")) +  # Rename x-axis labels
  theme_minimal() +  # Minimal theme for a clean look
  labs(
    x = "Timepoint",
    y = "Peak area of DCA-2-aminophenol",
    color = "Person ID",  # Label for the color legend
    linewidth = "Person ID"  # Label for the linewidth legend
  ) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text at 45 degrees for better readability
    axis.text.y = element_text(size = 12),  # Adjust the size of y-axis text for better readability
    axis.title = element_text(size = 14),  # Adjust the size of axis titles
    legend.position = "right"  # Keep the legend on the right side for clarity
  )

plot 

ggsave("Guma_boxplot_plus_lineplot_T2T3.svg", plot = plot, width = 8, height = 6, device = "svg")

ggsave("Guma_boxplot_plus_lineplot_T2T3.png", plot = plot, width = 8, height = 6, device = "png")

####Tests for normality of data
# Histogram of peak_area
ggplot(df_long_T2T3, aes(x = peak_area)) +
  geom_histogram(binwidth = 500, color = "black", fill = "skyblue") +
  labs(title = "Histogram of Peak Area",
       x = "Peak Area",
       y = "Frequency") +
  theme_minimal()
# Shapiro-Wilk test for normality
shapiro_test_result <- shapiro.test(df_long_T2T3$peak_area)
# View the p-value
shapiro_test_result$p.value #my value is 2.5e-08 so that data is not normally distributed - will use wilcox_test for stats

####Statistical analysis
# Perform a Wilcoxon test
wilcox_test_result <- wilcox.test(log_peak_area ~ timepoint, data = df_long_T2T3)

# View the result of the Wilcoxon test
wilcox_test_result$p.value