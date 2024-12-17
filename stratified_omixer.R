library(tidyverse)
library(omixerRpm)


################################################
####   ARGUMENTS GIVEN BY THE MAIN SCRIPT   ####
input <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]
db_string <- commandArgs(trailingOnly = TRUE)[3]

################################################



# Determine the set of modules used for the analysis
db <- loadDB(listDB()[grepl(pattern = db_string, x = listDB(), ignore.case = TRUE)])



###########################################
####   POSSIBLE MODIFIABLE ARGUMENTS   ####

minimum.coverage <-  0.01
score.estimator <-  "sum"
###########################################






# READING WOLTKA'S INPUT
print("IMPORTING PICRUSt2 RESULTS AND RE-FORMATTING")
input_file <- paste0(input, "/pred_metagenome_contrib.tsv.gz")

omixer_input <- read_tsv(file = input_file) %>%
  select(sample, taxon, `function`, taxon_function_abun) %>% 
  rename(taxa = taxon, ko = `function`) %>% 
  pivot_wider(id_cols = c(taxa, ko), 
              names_from = "sample",
              values_from = "taxon_function_abun", 
              values_fill = 0)



# Running omixer one sample at a time
final_abundance_table <- data.frame(taxa = character(), module_number = character(), module_name = character())
final_coverage_table <- data.frame(taxa = character(), module_number = character(), module_name = character())


# starts from the column 3 because column 1 and 2 has taxonomy and gene, respectively
print("start")
for (col in 3:ncol(omixer_input)) {
  
  sample_df <- omixer_input[, c(1,2,col)]
  
  sample_rpm <- rpm(x = sample_df, 
                    minimum.coverage = minimum.coverage,
                    module.db = db,
                    annotation = 2,
                    score.estimator = score.estimator)
  
  # Calculating the abundance of each module for this sample
  sample_module_abundance <- asDataFrame(sample_rpm, type = "abundance")
  # the third column will always be the column with the sample name
  colnames(sample_module_abundance) <- c("taxa", "module_number", "module_name", colnames(sample_df)[3])
  

  # force the columns to be character to avoid the clash of classes
  sample_module_abundance <- sample_module_abundance %>%
	mutate(taxa = as.character(taxa),
	       module_number = as.character(module_number),
	       module_name = as.character(module_name))


  
  # Merging the abundance of the modules in this sample with the previous analaysed
  final_abundance_table <- full_join(x = final_abundance_table, 
                                     y = sample_module_abundance, 
                                     by = c("taxa", "module_number", "module_name"))
  
  
  
  # Calculating the coverage of each module for this sample
  sample_module_coverage <- asDataFrame(sample_rpm, type = "coverage")
  # the third column will always be the column with the sample name
  colnames(sample_module_coverage) <- c("taxa", "module_number", "module_name", colnames(sample_df)[3])
  

  # force the columns to be character to avoid the clash of classes
  sample_module_coverage <- sample_module_coverage %>%
	mutate(taxa = as.character(taxa),
	       module_number = as.character(module_number),
	       module_name = as.character(module_name))  


  # Merging the coverage of the modules in this sample with the previous analaysed
  final_coverage_table <- full_join(x = final_coverage_table, 
                                    y = sample_module_coverage, 
                                    by = c("taxa", "module_number", "module_name"))
  
}

# Export abundance table
abudance_file <- paste0(output, "/stratified_", db_string, "_abundance", "_threshold", minimum.coverage, ".tsv")
write.table(x = final_abundance_table,
            file = abudance_file, 
            row.names = FALSE,
            sep = "\t")


# Export coverage table
coverage_file <- paste0(output, "/stratified_", db_string, "_coverage", "_threshold", minimum.coverage, ".tsv")
write.table(x = final_coverage_table,
            file = coverage_file, 
            row.names = FALSE,
            sep = "\t")
