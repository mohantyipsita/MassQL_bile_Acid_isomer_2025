####Script for Figure2e - MIND diet study#########
#GNPS2 job final link: https://gnps2.org/status?task=8bfb999be1ee4747bb78d9c7015e34b5

setwd("C:/Users/imoha/OneDrive - University of California, San Diego Health/Bile_Acid_Fragmentation_tree/CASE_U19_MIND")

library(stringr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(data.table)
library(ggplot2)
library(tibble)
library(gridExtra)

#function to transpose a feature table
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
df_libhit <- fread("C:/Users/imoha/OneDrive - University of California, San Diego Health/Bile_Acid_Fragmentation_tree/CASE_U19_MIND/MIND_FBMN_libhit.tsv") 
df_libhit_req <- df_libhit[,c(2,15,32)]

#read the FBMN table 
df_feature <- read.csv("MIND_featuretable.csv", check.names = FALSE)
df_feature_req <- df_feature[,-c(2:3)]

#remove columns that are not needed
df_feature_req_filtered <- df_feature_req %>%
  dplyr::select(1, dplyr::matches("_BL|_Y1|_Y2|_Y3"))

#substitute the name
df_feature_req_filtered <- df_feature_req_filtered %>%
  rename_with(~ ifelse(.x == names(df_feature_req_filtered)[1], .x, sub("\\.mzML Peak area$", "", .x)), everything())

#metadata
df_metadata_long <- fread("Metadata/Metadata_longitudinal.csv")
df_lucas_metadata <- fread("Metadata/Lucas_MINS_q2-ready-md.tsv")

df_lucas_metadata <- df_lucas_metadata %>%
  slice(-1) %>%
  mutate(mid = str_extract(`#SampleID`, "(?<=\\.)\\d{6}\\.\\w{2}(?=\\.)")) %>%
  mutate(mid = str_replace(mid, "\\.", "_"))

#read the food frequency questionnaire 
df_ffq_Y2 <- fread("Metadata/DUA_M24minddsserv.csv")
df_ffq_Y2_whlgrain <- df_ffq_Y2 %>%
  dplyr::select(1, dplyr::matches("wholegrain|whole_grains"))
# First, create the 'id' column
df_ffq_Y2_whlgrain <- df_ffq_Y2_whlgrain %>%
  mutate(sample_id = paste0(mid, "_Y2"))
# Then, subset to keep only the columns you need
df_ffq_Y2_whlgrain <- df_ffq_Y2_whlgrain[, c(4, 2)]


df_ffq_Y3 <- fread("Metadata/DUA_M36minddsserv.csv")
df_ffq_Y3_whlgrain <- df_ffq_Y3 %>%
  dplyr::select(1, dplyr::matches("wholegrain|whole_grains"))
# First, create the 'id' column
df_ffq_Y3_whlgrain <- df_ffq_Y3_whlgrain %>%
  mutate(sample_id = paste0(mid, "_Y3"))
# Then, subset to keep only the columns you need
df_ffq_Y3_whlgrain <- df_ffq_Y3_whlgrain[, c(4, 2)]


df_metadata_long_whlgrain <-  df_metadata_long %>%
  dplyr::select(1, dplyr::matches("wholegrain|whole_grains"))
df_metadata_long_whlgrain <- df_metadata_long_whlgrain[,c(1,3)]

#merge BL,Y1,Y2 and Y3
df_all_whlgrain <- df_metadata_long_whlgrain %>%
  left_join(df_ffq_Y2_whlgrain, by = "sample_id") %>%
  left_join(df_ffq_Y3_whlgrain , by = "sample_id")


#df_all_whlgrain_merged <- df_all_whlgrain %>%
mutate(whlgrain_wk_all = coalesce(wholegrain_wk, m24wholegrain_wk, m36wholegrain_wk)) %>%
  select(sample_id, whlgrain_wk_all)

df_all_whlgrain_merged <- df_all_whlgrain %>%
  mutate(whlgrain_wk_all = coalesce(wholegrain_wk, m24wholegrain_wk, m36wholegrain_wk))

df_all_whlgrain_merged <- df_all_whlgrain_merged[,c(1,5)]

####################################################################################
####I will plot box plots for those with labels#############
###################################################################################

#subset for feature ID #58941 which is the DCA-2-aminophenol
#df_DCA_2_aminophenol <- df_feature_req_filtered %>%
dplyr::filter(.[[1]] == 58941)

df_DCA_2_aminophenol <- df_feature_req_filtered %>%
  dplyr::filter(.[[1]] %in% c(58941, "355", "2825"))

t_df_DCA_2_aminophenol <- transpose_and_format_df(df_DCA_2_aminophenol)
t_df_DCA_2_aminophenol <- t_df_DCA_2_aminophenol[,-c(2,3)]


#merge with metadata
df_DCA_2_aminophenol_wholegrain <- merge (t_df_DCA_2_aminophenol, df_all_whlgrain_merged, by.x = "filename", by.y = "sample_id", all = FALSE)
#write.csv(df_DCA_2_aminophenol_wholegrain, "Foodomics/Whlgrain_metadata.csv", row.names = FALSE)

#read the intervention status
df_interveton <- fread("Metadata/DUA_intervention_status.csv")
df_other_vars_onlyBL <- fread("Metadata/METADATA_BASELINE_KEY_VARS_SUBSET.csv")

##Baseline
df_DCA_2_aminophenol_wholegrain_BL <-  df_DCA_2_aminophenol_wholegrain %>%
  dplyr::filter(grepl("_BL", filename))
#df_DCA_2_aminophenol_wholegrain_BL_2 <- df_DCA_2_aminophenol_wholegrain_BL[,c(1,2,4)]
df_DCA_2_aminophenol_wholegrain_BL_2  <- df_DCA_2_aminophenol_wholegrain_BL  %>%
  mutate(ID = sub("_.*", "", filename))
#merge with intervention status
df_DCA_2_aminophenol_wholegrain_BL_2 <- merge(df_DCA_2_aminophenol_wholegrain_BL_2, df_interveton, by.x = "ID", by.y = "mid", all.x = TRUE)
#baseline merge with other vars
df_DCA_2_aminophenol_wholegrain_BL_2 <- merge(df_DCA_2_aminophenol_wholegrain_BL_2, df_other_vars_onlyBL, by.x = "ID", by.y = "mid", all.x = TRUE)


##Year1
df_DCA_2_aminophenol_wholegrain_Y1 <-  df_DCA_2_aminophenol_wholegrain %>%
  dplyr::filter(grepl("_Y1", filename))
#df_DCA_2_aminophenol_wholegrain_Y1_2 <- df_DCA_2_aminophenol_wholegrain_Y1[,c(1,2,4)]
df_DCA_2_aminophenol_wholegrain_Y1_2  <- df_DCA_2_aminophenol_wholegrain_Y1  %>%
  mutate(ID = sub("_.*", "", filename))
#merge with intervention status
df_DCA_2_aminophenol_wholegrain_Y1_2 <- merge(df_DCA_2_aminophenol_wholegrain_Y1_2, df_interveton, by.x = "ID", by.y = "mid", all.x = TRUE)


##Year2
df_DCA_2_aminophenol_wholegrain_Y2 <-  df_DCA_2_aminophenol_wholegrain %>%
  dplyr::filter(grepl("_Y2", filename))
#df_DCA_2_aminophenol_wholegrain_Y1_2 <- df_DCA_2_aminophenol_wholegrain_Y1[,c(1,2,4)]
df_DCA_2_aminophenol_wholegrain_Y2_2  <- df_DCA_2_aminophenol_wholegrain_Y2  %>%
  mutate(ID = sub("_.*", "", filename))
#merge with intervention status
df_DCA_2_aminophenol_wholegrain_Y2_2 <- merge(df_DCA_2_aminophenol_wholegrain_Y2_2, df_interveton, by.x = "ID", by.y = "mid", all.x = TRUE)



##Year3
df_DCA_2_aminophenol_wholegrain_Y3 <-  df_DCA_2_aminophenol_wholegrain %>%
  dplyr::filter(grepl("_Y3", filename))
#df_DCA_2_aminophenol_wholegrain_Y1_2 <- df_DCA_2_aminophenol_wholegrain_Y1[,c(1,2,4)]
df_DCA_2_aminophenol_wholegrain_Y3_2  <- df_DCA_2_aminophenol_wholegrain_Y3  %>%
  mutate(ID = sub("_.*", "", filename))
#merge with intervention status
df_DCA_2_aminophenol_wholegrain_Y3_2 <- merge(df_DCA_2_aminophenol_wholegrain_Y3_2, df_interveton, by.x = "ID", by.y = "mid", all.x = TRUE)


#merging BL and T1
df_BL_sub <- df_DCA_2_aminophenol_wholegrain_BL_2[,c(1,3,4,7)]
df_Y1_sub <- df_DCA_2_aminophenol_wholegrain_Y1_2[,c(1,3,4)]

df_BL_sub <- df_BL_sub %>%
  rename(
    BL_58941 = `58941`,
    BL_whlgrain_wk_all = whlgrain_wk_all
  )

df_Y1_sub <- df_Y1_sub %>%
  rename(
    T1_58941 = `58941`,
    T1_whlgrain_wk_all = whlgrain_wk_all
  )

df_BL_Y1 <- merge(df_BL_sub, df_Y1_sub, by = "ID", all = FALSE)

#calculate delta area and delta whlgrain
df_BL_Y1  <- df_BL_Y1 %>%
  mutate(delta_whlgrain = T1_whlgrain_wk_all - BL_whlgrain_wk_all,  # Calculate the difference in grain_score
         delta_58941 = T1_58941 - BL_58941)

df_BL_Y1_nozero <- df_BL_Y1 %>%
  filter(delta_58941 != 0)

########Figure between Mind group BL and Y1 and ctrl group BL and Y1
df_mind_BLY1 <- df_BL_Y1 %>% filter(assignment == "MIND Diet Group")
df_ctrl_BLY1 <- df_BL_Y1 %>% filter(assignment == "Control Diet Group")

# Reshape the dataframe to long format for plotting
df_long_mind <- df_mind_BLY1 %>%
  dplyr::select(ID, BL_58941, T1_58941, delta_58941) %>%
  pivot_longer(cols = c(BL_58941, T1_58941), 
               names_to = "Timepoint", 
               values_to = "Peak_Area")

# Modify the Timepoint labels for clarity
df_long_mind <- df_long_mind %>%
  mutate(Timepoint = ifelse(Timepoint == "BL_58941", "BL", "Y1"))

df_long_mind <- df_long_mind %>%
  mutate(delta_class = ifelse(delta_58941 > 0, "Positive", "Negative"))


######remove samples where delta_peak_area is o
df_long_mind_nonzero <- df_long_mind %>% filter(delta_58941 != 0)

# Add a new column with log10 of Peak_Area
df_long_mind_nonzero$log_peak_area <- log10(df_long_mind_nonzero$Peak_Area + 1)

plot <- ggplot(df_long_mind_nonzero, aes(x = Timepoint, y = Peak_Area)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, size = 1) +  # Create box plots side by side and remove outliers
  geom_jitter(aes(color = Timepoint), position = position_dodge(width = 0.8), size = 4.5, alpha = 0.8) +  # Add jittered points
  geom_line(aes(group = ID, color = delta_class), alpha = 0.6, size = 1.2) +  # Color lines based on delta_class
  scale_x_discrete(labels = c("BL" = "Before MIND diet", "Y1" = "After MIND diet")) +  # Rename x-axis labels
  scale_color_manual(values = c("Positive" = "#987D9A", "Negative" = "#EECEB9")) +  # Color mapping for positive/negative changes
  theme_minimal() +  # Minimal theme for a clean look
  labs(
    x = "Timepoint",
    y = "Peak area of DCA-2-aminophenol",
    color = "Change Type"  # Legend label for the line color
  ) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text at 45 degrees for readability
    axis.text.y = element_text(size = 12),  # Adjust the size of y-axis text for readability
    axis.title = element_text(size = 14),  # Adjust the size of axis titles
    legend.position = "right",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

plot

ggsave("U19_boxplot_plus_lineplot_BLY1.svg", plot = plot, width = 8, height = 6, device = "svg")

ggsave("U19_boxplot_plus_lineplot_BLY1.png", plot = plot, width = 8, height = 6, device = "png")

# Reshape the data into wide format
df_wide_mind_nonzero <- df_long_mind_nonzero %>%
  dplyr::select(ID, Timepoint, Peak_Area) %>%
  pivot_wider(names_from = Timepoint, values_from = Peak_Area) %>%
  filter(!is.na(BL), !is.na(Y1))  # Ensure no missing values for BL or T3

# Perform a Wilcoxon signed-rank test on the paired data (BL vs Y1)  
wilcox_test_result <- wilcox.test(df_wide_mind_nonzero$BL, df_wide_mind_nonzero$Y1, paired = TRUE)

# View the result of the Wilcoxon test
wilcox_test_result$p.value
