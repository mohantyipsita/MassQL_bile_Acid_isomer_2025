######Script for Figure 2f - HIV-HNRC study#########
#GNPS2 job link: https://gnps2.org/status?task=c68fc90125424231bece7655426f9e17

setwd("C:/Users/imoha/OneDrive - University of California, San Diego Health/Bile_Acid_Fragmentation_tree/CASE_HNRC")

#function to transpose
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

#function to barplot
create_bar_plot <- function(data, x_column, fill_column = NULL, gradient = 2) {
  # Check if the specified columns exist in the dataframe
  if (!x_column %in% names(data)) {
    stop("The specified x_column does not exist in the dataframe.")
  }
  if (!is.null(fill_column) && !fill_column %in% names(data)) {
    stop("The specified fill_column does not exist in the dataframe.")
  }
  
  # Create the bar plot
  p <- ggplot(data, aes_string(x = x_column, fill = fill_column)) +
    geom_bar(position = "dodge") +
    labs(title = paste("Bar Plot of", x_column, ifelse(is.null(fill_column), "", paste("by", fill_column))),
         x = x_column,
         y = "Count") +
    scale_y_continuous(breaks = seq(0, max(table(data[[x_column]])), by = gradient)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}


###loading libraries
library(dplyr)
library(pheatmap)
library(stringr)
library(mixOmics)
library(ggpubr)
library(CoDaSeq)
library(vegan)
library(tibble)
library(svglite)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(impute)
library(data.table)

df_feature <- read.csv("HNRC_U19_test2_quant.csv", check.names = FALSE) #read the feature table
df_feature_req <- df_feature[,-c(2:6,9:13)] #subset for required columns

df_feature_req2 <- df_feature_req %>%
  dplyr::select(-which(colnames(df_feature_req) == "")) #filter out blank columns

df_feature_req_filtered <- df_feature_req2 %>%
  dplyr::select(-dplyr::matches("QC")) #filter out QC files

df_feature_req_filtered <- df_feature_req_filtered %>%
  rename_with(~ ifelse(.x == names(df_feature_req_filtered)[1], .x, sub("\\.mzML Peak area$", "", .x)), everything())


###############################################
###Testing for how many adduct forms of BAs####
df_libhit <-  fread("HNRC_FBMN_libhit.csv")
df_libhit_Req <- df_libhit [,-c(1,3:14,16:31,33:46)]

df_merge <- merge(df_feature_req_filtered, df_libhit, by.x = "row ID", by.y = "Scan", all.x = TRUE)

#subset for only those in GNPS-BILE-ACID-MODIFICATIONS library
df_merge_BA <- df_merge %>% filter(Organism == "GNPS-BILE-ACID-MODIFICATIONS" | Organism == "BILELIB19") #subset for Bile acid matches in the GNPS library
df_merge_BA_sub <- df_merge_BA %>% dplyr::select(c("row ID","correlation group ID"))
df_merge_BA_sub <- merge(df_merge_BA_sub, df_libhit, by.x = "row ID", by.y = "Scan", all = FALSE) 
df_merge_BA_sub <- df_merge_BA_sub[,-c(4)]

#editing the full feature table
df_feature <- df_feature[,-c(2:13)]
names(df_feature)[1] <- "Scan"
df_feature <- df_feature[, -ncol(df_feature)]
df_feature <- df_feature %>%
  rename_with(~ ifelse(.x == names(df_feature)[1], .x, sub("\\.mzML Peak area$", "", .x)), everything())

#transpose the feature table
t_feature <- transpose_and_format_df(df_feature)

#subset for the feature belonging to DCA-2-aminophenol
new_df <- t_feature[, c("filename", "25844")]

#filter for only biological samples
new_df_filtered <- new_df %>% dplyr::filter((str_detect(filename, "Sample"))) %>% mutate(SampleID = gsub(".*_(Sample_\\w+).*", "\\1", filename))

#rename the feature column to DCA_aminophenol
colnames(new_df_filtered)[2] <- "DCA_aminophenol"

write.csv(new_df_filtered, "DCA_2_aminophenol_peak_area.csv", row.names = FALSE)

#read metadata from metabolomics
df_metadata <- fread("Helena_metadata.csv") %>% mutate(SampleID = gsub(".*_(Sample_\\w+).*", "\\1", filename))

#read metadata from HNRC with more variables
df_large_metadata <- fread("hnrp_mibi_master_metadata_072424.csv") 
df_large_metadata_sub <- df_large_metadata %>%
  dplyr::select(`sample name`, host_age, host_sex, host_body_mass_index, cd4cd8_ratio, glucose_serum, mwgrain)

#subsetting the metadata for those samples that have mwgrain value (keep 0,1,2,3,4)
df_large_metadata_sub_mwgrain <- df_large_metadata_sub %>%
  filter(mwgrain %in% c(0, 1, 2, 3, 4))
df_large_metadata_sub_mwgrain <- df_large_metadata_sub_mwgrain %>%
  mutate(`sample name` = gsub("\\.", "_", `sample name`))

#merge both 
df_final_metadata <- merge(df_metadata, df_large_metadata_sub_mwgrain, by.x = "SampleID", by.y = 'sample name', all.y = TRUE)


#merge the the feature 25844 with the metadata
df_25844 <- merge(new_df_filtered, df_final_metadata, by = "SampleID", all.x = TRUE)
df_25844 <- df_25844[,-c(4)]

df_25844 <- df_25844 %>%
  mutate(Vial = sub("_Sample_.*", "", filename.x))

#df_25844_RT <- df_25844 %>%
filter(DCA_aminophenol!= 0) %>% filter(grepl("P3|P4", Vial)) 

######Plotting the number of samples where DCA-2-aminophenol was detected 
df_25844_sub <- df_25844[,c(3,51)]

df_25844_sub <- df_25844_sub %>%
  filter(!is.na(mwgrain))

#preparing datframe
df_25844_summary <- df_25844_sub %>%
  group_by(mwgrain) %>%  # Group by mwgrain only
  summarise(
    grain_group = first(mwgrain),  # mwgrain value
    number = sum(DCA_aminophenol != 0),  # Count of non-zero DCA_aminophenol values
    total = n(),  # Total number of rows for that grain group
    .groups = 'drop'  # Ensure that grouping is removed
  )

df_25844_summary <- df_25844_summary[,-c(1)]

df_25844_summary <- df_25844_summary %>%
  mutate(ratio = number / total)


plot <- ggplot(df_25844_summary, aes(x = as.numeric(grain_group), y = ratio)) +
  geom_point(size = 3, color = "#F95454", alpha = 0.7) +  # Scatter points
  geom_line(color = "#F95454", alpha = 0.7) +  # Line connecting the dots
  geom_smooth(method = "lm", se = FALSE, color = "#0D92F4", size = 0.8, linetype = "dashed") +  # Less bold, dashed regression line
  labs(x = "Grain Group", y = "Ratio", title = "Scatter Plot with Connected Dots and Dashed Regression Line") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black")  # Add x and y axes
  )

# Fit a linear model
lm_model <- lm(ratio ~ grain_group, data = df_25844_summary)

# View the summary of the regression model
summary(lm_model)


ggsave("HNRC_percentage_detected_whlgrain_lineplot.svg", plot = plot, width = 8, height = 6, device = "svg")

ggsave("HNRC_percentage_detected_whlgrain_lineplot.png", plot = plot, width = 8, height = 6, device = "png")

