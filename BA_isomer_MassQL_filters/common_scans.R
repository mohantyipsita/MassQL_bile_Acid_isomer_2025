####Use this script to obtain MS/MS scans for each terminal bin in the filtering tree#####

#read the tsv from the MassQL queries run on the python script "massql_main.py" deposited in GitHub
df_massql <- read.table('BA_isomers/MSV000094551_iimn_70perpeakheight_2_massql_ALL_NEW_water_loses.tsv', header = TRUE, sep = '\t')

#read the stage 2 query results - if the Stage 2 queries give an error, then this will need to be run separately outside the python script using the MassQL workflow on GNPS2 (https://www.nature.com/articles/s41592-025-02660-z)
#this step can be avoided if the Stage 2 queries run successfully on the python script
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
common_scans_mono_3a <- intersect(monohydroxy_scans, mono_3a_oh_scans)
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
common_scans_tri_3_6a_7a <- Reduce(intersect, list(trihydroxy_scans, tri_3_6_7_scans, tri_3_6a_7a_scans))
common_scans_tri_3_6b_7a_or_3_6a_7b <- Reduce(intersect, list(trihydroxy_scans, tri_3_6_7_scans, tri_3_6b_7a_or_3_6a_7b_scans))
common_scans_tri_3_6b_7a <- Reduce(intersect, list(trihydroxy_scans, tri_3_6_7_scans, tri_3_6b_7a_or_3_6a_7b_scans, tri_3_6b_7a_scans))


##the isomer bins that are not detected in this dataset are commented out
df_mono_3a <- data.frame(query = "Mono_3a", scan_list = common_scans_mono_3a)
df_mono_7b <- data.frame(query = "Mono_7b", scan_list = common_scans_mono_7b)

df_di_3_12a <- data.frame(query = "Di_3_12a", scan_list = common_scans_di_3_12a)
df_di_7_12a <- data.frame(query = "Di_7_12a", scan_list = common_scans_di_7_12a)
df_di_3_6 <- data.frame(query = "Di_3_6", scan_list = common_scans_di_3_6)
df_di_3_7 <- data.frame(query = "Di_3_7", scan_list = common_scans_di_3_7)
df_di_3or7_12b <- data.frame(query = "Di_3or7_12b", scan_list = common_scans_di_3or7_12b)

df_tri_ketone <- data.frame(query = "Tri_ketone", scan_list = common_scans_tri_ketone)
df_tri_3k <- data.frame(query = "Tri_3k", scan_list = common_scans_tri_3k)
df_tri_non_3k <- data.frame(query = "Tri_non_3k", scan_list = common_scans_non_3k)
df_tri_12k <- data.frame(query = "Tri_12k", scan_list = common_scans_tri_12k)
df_tri_3_6_7 <- data.frame(query = "Tri_3_6_7", scan_list = common_scans_tri_3_6_7)
df_tri_3_7a_12a <- data.frame(query = "Tri_3_7a_12a", scan_list = common_scans_tri_3_7a_12a)
df_tri_3_6a_7a <- data.frame(query = "Tri_3_6a_7a", scan_list = common_scans_tri_3_6a_7a)
df_tri_3_6b_7a_or_3_6a_7b <- data.frame(query = "Tri_3_6b_7a_or_3_6a_7b", scan_list = common_scans_tri_3_6b_7a_or_3_6a_7b) 
df_tri_3_6b_7a <- data.frame(query = "Tri_3_6b_7a", scan_list = common_scans_tri_3_6b_7a)

# Combine all dataframes into one df_new
df_new <- bind_rows(df_mono_3a, df_mono_7b, 
                    df_di_3_12a, df_di_7_12a, df_di_3_6, df_di_3_7, df_di_3or7_12b, 
                    df_tri_non_3k, df_tri_3k, df_tri_12k, df_tri_3_6_7, df_tri_3_7a_12a, df_tri_3_6a_7a, df_tri_3_6b_7a)


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
