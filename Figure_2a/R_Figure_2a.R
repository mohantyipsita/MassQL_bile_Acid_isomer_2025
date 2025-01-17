####Script for Figure2a - Paired intestinal and feces samples#########
#GNPS2 job final link: https://gnps2.org/status?task=4e5f76ebc4c6481aba4461356f20bc35

setwd("C:/Users/imoha/OneDrive - University of California, San Diego Health/Collaborations/Audrey")

library(stringr)
library(dplyr)
library(pheatmap)

#function to read and transpose a feature table
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

#read the original Bile acid libhit
df_libhit <- read.csv("FBMN_libhit_iimn_70peakheight_2.csv", check.names = FALSE) 
df_libhit_req <- df_libhit[,c(2,15,32)]

#read the FBMN feaure table
df_feature <- read.csv("MZMine/MSV000094551_iimn_70perpeakheight_2_quant.csv", check.names = FALSE)
df_feature_req <- df_feature[,-c(2:6,9:13)] #subset for required columns

#subset for only those files with the prefix "Gamme" or "QC"
df_feature_req_filtered <- df_feature_req %>%
  dplyr::select(-dplyr::matches("Gamme|QC"))

df_feature_req_filtered <- df_feature_req_filtered %>%
  rename_with(~ ifelse(.x == names(df_feature_req_filtered)[1], .x, sub("\\.mzXML Peak area$", "", .x)), everything())

#blank subtraction 
blank_substraction_edited <- function(df) {
  # Identify columns that contain "Blank" in their names
  blank_cols <- grep("Blank", names(df))
  
  # Ensure that blank columns are present
  if (length(blank_cols) == 0) {
    stop("No columns containing 'Blank' found.")
  }
  
  # Identify the columns to be averaged (excluding the first column)
  value_cols <- setdiff(seq_along(df), c(1, blank_cols))
  
  # Loop through each row to perform the required calculations
  df <- df %>%
    rowwise() %>%
    filter(mean(c_across(all_of(value_cols)), na.rm = TRUE) > 
             5 * mean(c_across(all_of(blank_cols)), na.rm = TRUE)) %>%
    ungroup()  # Ungroup to return to the original structure
  
  return(df)
}

#create the blank substracted feature table
df_feature_req_filtered_blksubs <- blank_substraction_edited(df_feature_req_filtered)
#write.csv(df_feature_req_filtered_blksubs, "Feature_table_blanksubstracted_IIMN_70perpeakheight_2.csv", row.names = FALSE)

#read the tsv from the MassQL queries run on python
df_massql <- read.table('BA_isomers/MSV000094551_iimn_70perpeakheight_2_massql_ALL_NEW_water_loses.tsv', header = TRUE, sep = '\t')

#manually read the stage 2 query results
df_stage2 <- read.table('BA_isomers/MSV000094551_iimn_70perpeakheight_2_massql_stage2_NEW_water_loses.tsv', header = TRUE, sep = '\t')

#combine both dataframes
df_comb <- rbind(df_stage2, df_massql)
df_comb_edit <- df_comb #duplicated to make edits on this dataframe
#write.csv(df_comb, "BA_isomers/All_Scans_massql.csv", row.names = FALSE)


########BILE ACID ISOMER LABEL ASSIGNMENT#######
# Convert the scan_list column from JSON-like string to an actual list
df_comb_edit$scan_list <- lapply(df_comb_edit$scan_list, function(x) {
  scans <- gsub("[\\[\\]]", "", x)  # Remove brackets
  unlist(strsplit(scans, ", "))     # Split by comma and space
})

# Extract the scan lists for "Monohydroxy" 
monohydroxy_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Monohydroxy_all_stage2"])
mono_3_oh_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Mono-3-OH"])
mono_3a_oh_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Mono-3a-OH"])
mono_7a_oh_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Mono-7a-OH"])
mono_7b_oh_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Mono-7b-OH"])

# Extract the scan lists for "Dihydroxy" 
dihydroxy_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Dihydroxy_all_stage2"])
di_3or7_12a_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Di-3,12a-OH; 7,12a-OH"])
di_non12a_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Di-3,7-OH; 3,6-OH; 3,12b-OH; 7,12b-OH"])
di_non_3_6_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Di-3,7-OH; 3,12b-OH; 7,12b-OH"])
di_3_12a_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Di-3,12a-OH"])
di_7_12a_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Di-7,12a-OH"])
di_3_6_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Di-3,6-OH"])
di_3_7_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Di-3,7-OH"])
di_3or7_12b_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Di-3,12b-OH; 7,12b-OH"])

# Extract the scan lists for "Trihydroxy" 
trihydroxy_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "Trihydroxy_all_stage2"])
tri_all_ketones_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "All_ketone_NEW"])
tri_3k_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "3k_7or12_OH"])
tri_12k_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "3a_12k"])
tri_3_6_7_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "3,6,7_NEW"])
tri_3_7a_12a_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "3,7a,12a-OH_NEW"])
tri_3_6b_7b_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "3,6b,7b-OH"])
tri_3_6a_7a_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "3,6a,7a-OH"])
tri_3_6b_7a_or_3_6a_7b_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "3,6b,7a-OH_3,6a,7b-OH"])
tri_3_6b_7a_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "3,6b,7a-OH"])
tri_3_6a_7b_scans <- unlist(df_comb_edit$scan_list[df_comb_edit$query == "3,6a,7b-OH"])

# Find the common scanns tracing the branch of the MassQL tree
common_scans_mono_3_oh <- intersect(monohydroxy_scans, mono_3_oh_scans)
common_scans_mono_3a <- intersect(monohydroxy_scans, mono_3a_oh_scans)
common_scans_mono_7a <- intersect(monohydroxy_scans, mono_7a_oh_scans)
common_scans_mono_7b <- intersect(monohydroxy_scans, mono_7b_oh_scans)

common_scans_di_3_12a <- Reduce(intersect, list(dihydroxy_scans, di_3or7_12a_scans, di_3_12a_scans))
common_scans_di_7_12a <- Reduce(intersect, list(dihydroxy_scans, di_3or7_12a_scans, di_7_12a_scans))
common_scans_di_3_6 <- Reduce(intersect, list(dihydroxy_scans, di_non12a_scans, di_3_6_scans))
common_scans_di_3_7 <- Reduce(intersect, list(dihydroxy_scans, di_non12a_scans, di_non_3_6_scans, di_3_7_scans))
common_scans_di_3or7_12b <- Reduce(intersect, list(dihydroxy_scans, di_non12a_scans, di_3or7_12b_scans))

common_scans_tri_ketone <- Reduce(intersect, list(trihydroxy_scans, tri_all_ketones_scans))
common_scans_tri_3k <- Reduce(intersect, list(trihydroxy_scans, tri_all_ketones_scans, tri_3k_scans))
common_scans_non_3k <- setdiff(Reduce(intersect, list(trihydroxy_scans, tri_all_ketones_scans)), tri_3k_scans)
common_scans_tri_12k <- Reduce(intersect, list(trihydroxy_scans, tri_all_ketones_scans, tri_12k_scans))
common_scans_tri_3_6_7 <- Reduce(intersect, list(trihydroxy_scans, tri_3_6_7_scans))
common_scans_tri_3_7a_12a <- Reduce(intersect, list(trihydroxy_scans, tri_3_7a_12a_scans))
common_scans_tri_3_6b_7b <- Reduce(intersect, list(trihydroxy_scans, tri_3_6_7_scans, tri_3_6b_7b_scans))
common_scans_tri_3_6a_7a <- Reduce(intersect, list(trihydroxy_scans, tri_3_6_7_scans, tri_3_6a_7a_scans))
common_scans_tri_3_6b_7a_or_3_6a_7b <- Reduce(intersect, list(trihydroxy_scans, tri_3_6_7_scans, tri_3_6b_7a_or_3_6a_7b_scans))
common_scans_tri_3_6b_7a <- Reduce(intersect, list(trihydroxy_scans, tri_3_6_7_scans, tri_3_6b_7a_or_3_6a_7b_scans, tri_3_6b_7a_scans))
common_scans_tri_3_6a_7b <- Reduce(intersect, list(trihydroxy_scans, tri_3_6_7_scans, tri_3_6b_7a_or_3_6a_7b_scans, tri_3_6a_7b_scans))

##the isomer bins that are not detected in this dataset are commented out
df_mono_3a <- data.frame(query = "Mono_3a", scan_list = common_scans_mono_3a)
#df_mono_7b <- data.frame(query = "Mono_7b", scan_list = common_scans_mono_7b)

df_di_3_12a <- data.frame(query = "Di_3_12a", scan_list = common_scans_di_3_12a)
#df_di_7_12a <- data.frame(query = "Di_7_12a", scan_list = common_scans_di_7_12a)
df_di_3_6 <- data.frame(query = "Di_3_6", scan_list = common_scans_di_3_6)
df_di_3_7 <- data.frame(query = "Di_3_7", scan_list = common_scans_di_3_7)
df_di_3or7_12b <- data.frame(query = "Di_3or7_12b", scan_list = common_scans_di_3or7_12b)

df_tri_ketone <- data.frame(query = "Tri_ketone", scan_list = common_scans_tri_ketone)
df_tri_3k <- data.frame(query = "Tri_3k", scan_list = common_scans_tri_3k)
df_tri_non_3k <- data.frame(query = "Tri_non_3k", scan_list = common_scans_non_3k)
df_tri_12k <- data.frame(query = "Tri_12k", scan_list = common_scans_tri_12k)
df_tri_3_6_7 <- data.frame(query = "Tri_3_6_7", scan_list = common_scans_tri_3_6_7)
df_tri_3_7a_12a <- data.frame(query = "Tri_3_7a_12a", scan_list = common_scans_tri_3_7a_12a)
df_tri_3_6b_7b <- data.frame(query = "Tri_3_6b_7b", scan_list = common_scans_tri_3_6b_7b)
#df_tri_3_6a_7a <- data.frame(query = "Tri_3_6a_7a", scan_list = common_scans_tri_3_6a_7a)
#df_tri_3_6b_7a_or_3_6a_7b <- data.frame(query = "Tri_3_6b_7a_or_3_6a_7b", scan_list = common_scans_tri_3_6b_7a_or_3_6a_7b) 
#df_tri_3_6b_7a <- data.frame(query = "Tri_3_6b_7a", scan_list = common_scans_tri_3_6b_7a)
#df_tri_3_6a_7b <- data.frame(query = "Tri_3_6a_7b", scan_list = common_scans_tri_3_6a_7b)

# Combine all dataframes into one df_new
df_new <- bind_rows(df_mono_3a, 
                    df_di_3_12a, df_di_3_6, df_di_3_7,df_di_3or7_12b, 
                    df_tri_non_3k, df_tri_3k, df_tri_12k, df_tri_3_6_7, df_tri_3_7a_12a, df_tri_3_6b_7b)

df_new2 <- bind_rows(df_mono_3a, 
                     df_di_3_12a, df_di_3_6, df_di_3_7,df_di_3or7_12b, 
                     df_tri_ketone, df_tri_3_6_7, df_tri_3_7a_12a, df_tri_3_6b_7b)


#combine all the scans into a one row with all scan lists
df_new <- df_new %>%
  group_by(query) %>%
  summarize(scan_list = paste("[", paste(scan_list, collapse = ", "), "]", sep = ""))

# Initialize an empty dataframe to store the results
MassQL_scan <- data.frame(Scan = integer(), class = character(), stringsAsFactors = FALSE)

# Iterate over each row of the dataframe
library(jsonlite)
for (i in 1:nrow(df_new)) {
  # Extract the query and scan list for the current row
  query <- df_new$query[i]
  
  # Check if scan_list is NA and skip to the next iteration if it is
  if (is.na(df_new$scan_list[i])) {
    next
  }
  
  # Convert the scan_list string to a vector using fromJSON
  scan_list <- fromJSON(df_new$scan_list[i])
  
  # Create a temporary dataframe with Scan and class columns
  temp_df <- data.frame(Scan = scan_list, class = query, stringsAsFactors = FALSE)
  
  # Append the temporary dataframe to the result dataframe
  MassQL_scan <- bind_rows(MassQL_scan, temp_df)
}

#remove the scans that have more than one class label
MassQL_scan_unique <- MassQL_scan[!duplicated(MassQL_scan$Scan) & !duplicated(MassQL_scan$Scan, fromLast = TRUE), ]

#merge feature table with libhit
df_merge <- merge(df_feature_req_filtered_blksubs, df_libhit, by.x = "row_ID", by.y = "Scan", all.x = TRUE)
#subset for the matches to bile acid library
df_merge_BA <- df_merge %>% filter(Organism == "GNPS-BILE-ACID-MODIFICATIONS"| Organism == "BILELIB19")

# Extract the Scan column values from MassQL_scan_unique
scan_values <- MassQL_scan_unique$Scan

# Identify rows from df_feature_req_filtered_blksubs that are in scan_values but absent in df_merge_BA - this is to add those BA annotations from libraries other the GNPS modification and BILELIB19
additional_rows <- df_feature_req_filtered_blksubs %>%
  filter(row_ID %in% scan_values & !(row_ID %in% df_merge_BA$row_ID))

# Append the additional rows to df_merge -- total 1499 BA libhits
df_merge_BA_new <- bind_rows(df_merge_BA, additional_rows)

#subset for required columns
df_merge_BA_new <- df_merge_BA_new[,-c(2,3,17,33:77)]

t_df_merge_BA <- transpose_and_format_df(df_merge_BA_new)

#read the metadata
df_metadata <- read.csv("Audrey_metadata_updated.csv", check.names = FALSE)

#merge with metadata
df_merge_BA_meta <- merge(t_df_merge_BA, df_metadata, by = "filename", all = FALSE)
df_merge_BA_meta <- df_merge_BA_meta[,-c(1501:1504,1506:1507)]
df_merge_BA_meta_2 <- df_merge_BA_meta[,c(1, 1501, 2:1500)]
names(df_merge_BA_meta_2)[2] <- "Label"

#write.csv(df_merge_BA_meta_2, "MetaboanalystR/Input_final_metaboanalystR_new_waterlosses.csv", row.names = FALSE)

###############################################################
##############VOLCANO PLOT WITH MetaboanalystR#################
############################################################### used for paper final figure

#library(EnhancedVolcano)
library(MetaboAnalystR)

mSet<-InitDataObjects("pktable", "stat", FALSE)

mSet<-Read.TextData(mSet, "MetaboanalystR/Input_final_metaboanalystR_new.csv", "rowu", "disc")
mSet<-SanityCheckData(mSet)

mSet<-PreparePrenormData(mSet) #You have manually change the name from data_original.qs to data_proc.qs 
mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_nonscaled", format ="png", dpi=72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_nonscaled", format = "png", dpi=72, width=NA)

# Perform fold-change analysis on uploaded data, unpaired -- # this analysis is removing 157 features
mSet<-FC.Anal(mSet, 1, 0, FALSE)

# Plot fold-change analysis
mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)

# To view fold-change 
mSet$analSet$fc$fc.log

# Perform T-test (parametric)
mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, "fdr", TRUE)

# Plot of the T-test results
mSet<-PlotTT(mSet, imgName = "ttest_0_", format = "png", dpi = 72, width=NA)

# Perform the volcano analysis
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.05, TRUE, "fdr")

# Create the volcano plot
mSet<-PlotVolcano(mSet, "volcano_0_", 1, 0, format ="png", dpi=72, width=NA)


##merge the foldchange and the p-value tables
df_foldchange <- read.csv("MetaboanalystR/FINAL/fold_change.csv", check.names = FALSE)
df_ttest <- read.csv("MetaboanalystR/FINAL/t_test_all.csv", check.names = FALSE)

df_volcano_all <- merge(df_foldchange, df_ttest, by = "Scan", all = FALSE)
df_volcano_all_sub <- df_volcano_all[,c(1,3,6,7)]

#merge with libhits
df_libhit_sub <- df_libhit [,-c(1,3:14,16:31,33:46)]

#subset for only those in GNPS-BILE-ACID-MODIFICATIONS library
df_libhit_BA <- df_libhit_sub %>% filter(Organism == "GNPS-BILE-ACID-MODIFICATIONS"| Organism == "BILELIB19")
df_volcano_all_sub_annot <- merge(df_volcano_all_sub, df_libhit_BA, by = "Scan", all.x = TRUE)
df_volcano_all_sub_annot <- df_volcano_all_sub_annot[,-c(6)]


df_volcano_all_sub_annot_class <- merge(df_volcano_all_sub_annot, MassQL_scan_unique, by = "Scan", all.x = TRUE)
df_volcano_all_sub_annot_class[is.na(df_volcano_all_sub_annot_class)] <- "NA"

#write.csv(df_volcano_all_sub_annot_class, "MassQL_Scan_list_tocheck.csv", row.names = FALSE)

# Create a named vector of colors for each class
class_colors_mono <- c("Mono_3a" = "#FAA828", "NA" = "#758694")

class_colors_di <- c("Di_3_12a" = "#4D4B9F", "Di_3_6" = "#FAA828", "Di_3_7" = "#DC0C83", "Di_3or7_12b" = "#4DBBEB", "NA" = "#758694")

class_colors_tri <- c("Tri_3_7a_12a" = "#8ABFA3", "Tri_3_6_7" = "#C62E2E", "Tri_3k" = "#9D0191", "Tri_non_3k" ="#091057", "NA" = "#758694") 


# Create a named vector of shapes for each class
class_shapes_mono <- c("Mono_3a" = 18, "NA" = 1)

class_shapes_di <- c("Di_3_12a" = 15, "Di_3_6" = 15, "Di_3_7" = 17, "Di_3or7_12b" = 18, "NA" =1)

class_shapes_tri <- c("Tri_3_7a_12a" = 15, "Tri_3_6_7" = 15, "Tri_3k" = 17, "Tri_non_3k" = 17, "NA" = 1)


# Apply the colors based on the class column
df_volcano_all_sub_annot_class$color_mono <- class_colors_mono[df_volcano_all_sub_annot_class$class]
df_volcano_all_sub_annot_class$color_di <- class_colors_di[df_volcano_all_sub_annot_class$class]
df_volcano_all_sub_annot_class$color_tri <- class_colors_tri[df_volcano_all_sub_annot_class$class]

df_volcano_all_sub_annot_class$shape_mono <- class_shapes_mono[df_volcano_all_sub_annot_class$class]
df_volcano_all_sub_annot_class$shape_di <- class_shapes_di[df_volcano_all_sub_annot_class$class]
df_volcano_all_sub_annot_class$shape_tri <- class_shapes_tri[df_volcano_all_sub_annot_class$class]


# Check if lengths match
stopifnot(length(df_volcano_all_sub_annot_class$color) == nrow(df_volcano_all_sub_annot_class))
stopifnot(length(df_volcano_all_sub_annot_class$shape) == nrow(df_volcano_all_sub_annot_class))


# Creating a keyvals vector similar to the documentation example
keyvals <- df_volcano_all_sub_annot_class$color_di
names(keyvals) <- df_volcano_all_sub_annot_class$class

# Creating a named vector for shapes similar to the documentation example
shape_vals <- df_volcano_all_sub_annot_class$shape_di
names(shape_vals) <- df_volcano_all_sub_annot_class$class

# Create a custom size vector
point_size <- ifelse(df_volcano_all_sub_annot_class$class == "NA", 0.8, 3.5)

# Check the custom size vector
table(point_size)

# Create a vector of labels to display (those that are not "NA" in class)
labels_to_display <- df_volcano_all_sub_annot_class$Scan[df_volcano_all_sub_annot_class$class != "NA"]


#Plotting the volcano plot with EnhancedVolcano plot
library(EnhancedVolcano)

volcano_plot_di <- EnhancedVolcano(df_volcano_all_sub_annot_class,
                                   lab = df_volcano_all_sub_annot$Scan,
                                   x = 'log2FC',
                                   y = 'FDR',
                                   title = 'Fecal vs Instestine',
                                   pCutoff = 0.1,
                                   FCcutoff = 1,
                                   pointSize = point_size,
                                   labSize = 4,
                                   ylim = c(0, 8),
                                   #col=c('black', 'black', 'black', 'red3'),
                                   colCustom = keyvals,
                                   shapeCustom = shape_vals,
                                   colAlpha = 1,
                                   drawConnectors = FALSE,
                                   widthConnectors = 0.75,
                                   selectLab = labels_to_display,
                                   legendPosition = 'left',
                                   legendLabSize = 15,
                                   legendIconSize = 5.0,
                                   arrowheads = FALSE,
                                   gridlines.major = FALSE,
                                   gridlines.minor = FALSE,
                                   border = 'full',
                                   borderWidth = 1,
                                   borderColour = 'black')

volcano_plot_di

# Save the plot as SVG
ggsave("FINAL_volcano_plot_di.svg", plot = volcano_plot_di, width = 10, height = 8, units = "in", device = "svg")
ggsave("FINAL_volcano_plot_di.png", plot = volcano_plot_di, width = 10, height = 8, units = "in", device = "png")
