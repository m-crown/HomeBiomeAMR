suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(pheatmap))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(ontologyIndex))
suppressMessages(library(readr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(optparse))
suppressMessages(library(pander))
suppressMessages(library(tidyverse))
suppressMessages(library(readr))

option_list = list(
  make_option(c("-m", "--metadata"), type="character", default="/raid/HomeBiome/meta/metadata/homebiome_metagenomics_metadata.csv",
              help="metadata file", metavar="character"),
  make_option(c("-s", "--shortbred_dir"), type="character", default="/raid/HomeBiome/meta/snakemake_version/analysis/shortBRED",
              help="directory for shortbred analysis files to be found", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default='/raid/HomeBiome/meta/snakemake_version/analysis/shortBRED_analysis/',
              help="output_directory, must exist currently", metavar="character"),
  make_option(c("--aro_index"), type="character", default='/raid/HomeBiome/meta/snakemake_version/analysis/shortBRED_analysis/',
              help="aro_index", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Import the CARD ontology
aro_index = read_delim(file=opt$aro_index) %>%
  select(-c("CVTERM ID", "Model Sequence ID", "Model ID")) %>%
  distinct()
  
#set wd for shortBRED analysis and import metadata file
metadata <- read.csv(file = opt$metadata)

#getting the rpkm values for each aro from each sample and forming dataframe

shortbred_files=list.files(path=opt$shortbred_dir, recursive = TRUE, pattern="*_shortbred_results.tsv", full.names=TRUE)
cols <- sub(".*/([^/]+)_shortbred_results\\.tsv", "\\1", shortbred_files)

shortbred_rpkm <- list.files(path=opt$shortbred_dir, recursive = TRUE, pattern="*_shortbred_results.tsv", full.names=TRUE) %>% 
  lapply(function(x){read.table(pipe(paste("cut -f1,2 ", x)), row.names = 1, header=T)}) %>% 
  bind_cols

colnames(shortbred_rpkm) <- cols

#filter DF to only include rows (AROs) which have an RPKM value in one of the samples

shortbred_rpkm <- subset(shortbred_rpkm, select = -c(REANEG))
shortbred_rpkm <- shortbred_rpkm[rowSums(shortbred_rpkm) >0,] 

#Extracting the ARO from the rowname which includes full FASTA header.  

aros <- as.data.frame(str_split_fixed(rownames(shortbred_rpkm), "\\|", Inf))
aros <- aros$V3
aros <- gsub("_", ":", aros)
rownames(shortbred_rpkm) <- aros

#Removing ARO:3003361 which appears to be a contaminant - is highly abundant in all samples but is not a well known ARO - see card website
#shortbred_rpkm <- shortbred_rpkm[!(rownames(shortbred_rpkm) %in% c("ARO:3003361")),]

shortbred_rpkm_mean <- as.data.frame(rowMeans(shortbred_rpkm))

#making a sample column for plotting on heatmap 
sample_col <- data.frame(sample = head(metadata$Continent, -1))
row.names(sample_col) <- head(metadata$DSUK_ID, -1)

#merging the correct aro information with the rpkm df 
shortbred_rpkm <- merge(shortbred_rpkm, aro_index[c("ARO Accession", "ARO Name", "AMR Gene Family", "Drug Class", "Resistance Mechanism", "CARD Short Name")], by.x = 0, by.y = "ARO Accession")
rownames(shortbred_rpkm) <- shortbred_rpkm$Row.names

#aggregate the individual AROs into their family BY THEIR MEAN

agg <- shortbred_rpkm %>% 
  group_by(`AMR Gene Family`) %>% 
  #mutate(`AMR Gene Family` = ifelse(n() == 1, `ARO Name`, `AMR Gene Family`)) %>%
  summarise(across(intersect(names(shortbred_rpkm), metadata$DSUK_ID), ~ mean(.x, na.rm = TRUE))) %>% 
  mutate(rowsums = rowSums(select(., -`AMR Gene Family`))) %>%
  arrange(desc(rowsums)) %>%
  mutate(family = ifelse(row_number() > 10, "Other Arg", `AMR Gene Family`)) %>%
  mutate(numbered_family = paste0("family", row_number()))


agg_top10_plot <- agg %>%
  filter(family != "Other Arg") %>%
  select(-c("family", "rowsums", "numbered_family")) %>%
  mutate(`AMR Gene Family` = str_wrap(`AMR Gene Family`, width = 100)) %>%
  column_to_rownames(var = "AMR Gene Family") %>%
  pheatmap(., annotation_col = sample_col, fontsize_row = 12, angle_col = 45, fontsize_col = 12)

ggsave(paste0(opt$out_dir, "top10_family_aro.png"), agg_top10_plot, height = 16, width = 16, dpi = 600)

agg_top10_plot2 <- agg %>%
  filter(family != "Other Arg") %>%
  select(-c("family", "rowsums", "AMR Gene Family")) %>%
  column_to_rownames(var = "numbered_family") %>%
  pheatmap(., annotation_col = sample_col, fontsize_row = 12, fontsize_col = 12, angle_col = 45)

ggsave(paste0(opt$out_dir, "top10_family_aro2.png"), agg_top10_plot2, height = 16, width = 16, dpi = 600)

top10 <- shortbred_rpkm %>% 
  mutate(rowsums = rowSums(select(., intersect(names(shortbred_rpkm), metadata$DSUK_ID)))) %>%
  arrange(desc(rowsums)) %>%
  mutate(`CARD Short Name` = ifelse(row_number() > 10, "Other ARO", `CARD Short Name`))

top10_plot <- top10 %>%
  filter(`CARD Short Name` != "Other ARO") %>%
  remove_rownames() %>% 
  column_to_rownames(var = "CARD Short Name") %>%
  select(intersect(names(shortbred_rpkm), metadata$DSUK_ID)) %>%
  pheatmap(., annotation_col = sample_col, fontsize_row = 12, angle_col = 45, fontsize_col = 12)

ggsave(paste0(opt$out_dir, "top10_aro.png"), top10_plot, height = 16, width = 16, dpi = 600)

#Total antibiotic resistance gene families (in RPKMs) - so basically a sum of each of the gene families identified above per sample. 
print(length(rownames(shortbred_rpkm)))
print(length(unique(shortbred_rpkm$`AMR Gene Family`)))

per_sample_RPKM <- shortbred_rpkm %>%
  select(-c("Row.names", "CARD Short Name", "ARO Name", "AMR Gene Family", "Resistance Mechanism", "Drug Class")) %>%
  summarise(across(
    everything(),
    list(sum_rpkms = ~sum(., na.rm = TRUE)),
    .names = "{.col}"
  )) %>%
  pivot_longer(cols = starts_with("DSUK"), names_to = "DSUK_ID", values_to = "sum_rpkms") %>%
  left_join(metadata[c('DSUK_ID', "Continent")], by = "DSUK_ID") %>%
  rename(Location = Continent)

total_rpkms_country <- ggplot(per_sample_RPKM, aes(x=Location, y=sum_rpkms, fill=Location)) + 
  geom_violin() + 
  geom_boxplot(width=0.1) + geom_point() + 
  labs(title = "Total antimicrobial resistance by location", y = "Total Antimicrobial Resistance (RPKM)") +
  theme(legend.position = "none")

ggsave(paste0(opt$out_dir, "total_rpkms_country.png"), total_rpkms_country, height = 16, width = 12, dpi = 600)

save(shortbred_rpkm, file = paste0(opt$out_dir, "shortbred_rpkm.Rdata"))

shortbred_rpkm_krona_resistance <- shortbred_rpkm %>%
  mutate(rowsums = rowSums(select(., intersect(names(.), metadata$DSUK_ID)))) %>%
  group_by(`Resistance Mechanism`) %>%
  summarize(group_total = sum(rowsums)) %>%
  ungroup() %>%
  mutate(total_rowsums = sum(group_total)) %>%
  mutate(group_perc_total = group_total/total_rowsums) 

shortbred_rpkm_krona_resistance <- shortbred_rpkm_krona_resistance %>%
  left_join(shortbred_rpkm, by = "Resistance Mechanism") %>%
  mutate(`Resistance Mechanism` = if_else(group_perc_total <= 0.01, "other", as.character(`Resistance Mechanism`))) %>%
  mutate(rowsums = rowSums(select(., intersect(names(.), metadata$DSUK_ID)))) %>%
  select(c("rowsums","Resistance Mechanism","CARD Short Name"))

shortbred_rpkm_krona_drug <- shortbred_rpkm %>%
  mutate(rowsums = rowSums(select(., intersect(names(.), metadata$DSUK_ID)))) %>%
  group_by(`Drug Class`) %>%
  summarize(group_total = sum(rowsums)) %>%
  ungroup() %>%
  mutate(total_rowsums = sum(group_total)) %>%
  mutate(group_perc_total = group_total/total_rowsums) 

shortbred_rpkm_krona_drug <- shortbred_rpkm_krona_drug %>%
  left_join(shortbred_rpkm, by = "Drug Class") %>%
  mutate(`Drug Class` = if_else(group_perc_total <= 0.01, "other", as.character(`Drug Class`))) %>%
  mutate(rowsums = rowSums(select(., intersect(names(.), metadata$DSUK_ID)))) %>%
  select(c("rowsums","Drug Class","CARD Short Name"))

#writing tables for processing in krona
write.table(shortbred_rpkm_krona_resistance, paste0(opt$out_dir, "rpkms_sum_res_combined.tsv"), sep = "\t" ,row.names = FALSE, col.names = FALSE , quote = FALSE)
write.table(shortbred_rpkm_krona_drug, paste0(opt$out_dir, "rpkms_sum_drug_combined.tsv"), sep = "\t" ,row.names = FALSE, col.names = FALSE , quote = FALSE)

#adding the resistance mechanisms to the per sample files to do individual krona plots for each location
#also removing any zero AROs before writing

for (ci in intersect(names(shortbred_rpkm), metadata$DSUK_ID)) {
  subset_res = shortbred_rpkm %>%
    mutate(rowsums = rowSums(select(., intersect(names(.), metadata$DSUK_ID)))) %>%
    group_by(`Drug Class`) %>%
    summarize(group_total = sum(rowsums)) %>%
    ungroup() %>%
    mutate(total_rowsums = sum(group_total)) %>%
    mutate(group_perc_total = group_total/total_rowsums) 
  
  subset_res <- subset_res %>%
    left_join(shortbred_rpkm, by = "Drug Class") %>%
    mutate(`Resistance Mechanism` = if_else(group_perc_total <= 0.01, "other", as.character(`Resistance Mechanism`))) %>%
    mutate(rowsums = rowSums(select(., intersect(names(.), metadata$DSUK_ID)))) %>%
    select(c( all_of(ci), "Resistance Mechanism", "CARD Short Name")) %>%
    filter(!!sym(ci) > 0)
  
  subset_drug <- shortbred_rpkm %>%
    mutate(rowsums = rowSums(select(., intersect(names(.), metadata$DSUK_ID)))) %>%
    group_by(`Drug Class`) %>%
    summarize(group_total = sum(rowsums)) %>%
    ungroup() %>%
    mutate(total_rowsums = sum(group_total)) %>%
    mutate(group_perc_total = group_total/total_rowsums) 
  
  subset_drug <- subset_drug %>%
    left_join(shortbred_rpkm, by = "Drug Class") %>%
    mutate(`Resistance Mechanism` = if_else(group_perc_total <= 0.01, "other", as.character(`Resistance Mechanism`))) %>%
    mutate(rowsums = rowSums(select(., intersect(names(.), metadata$DSUK_ID)))) %>%
    select(c(all_of(ci), "Drug Class","CARD Short Name")) %>%
    filter(!!sym(ci) > 0)
  
  write.table(subset_drug, paste0(opt$out_dir, ci, "_drug.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE , quote = FALSE)
  write.table(subset_res, paste0(opt$out_dir, ci, "_res.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE , quote = FALSE)
}

file.create(paste0(opt$out_dir, "REANEG_drug.tsv"))
file.create(paste0(opt$out_dir, "REANEG_res.tsv"))



#produce a table of the amr gene families and percentage occurences - essentially a flat version of the krona plots

resistance_mechanism_table <- shortbred_rpkm %>%
  mutate(rowsums = rowSums(select(., intersect(names(.), metadata$DSUK_ID)))) %>%
  group_by(`Resistance Mechanism`) %>%
  summarise(total_count = sum(rowsums)) %>%
  mutate(total_sum = sum(total_count)) %>%
  mutate(percentage = (total_count/total_sum) * 100) %>%
  arrange(desc(percentage))

write.csv(resistance_mechanism_table, file = paste0(opt$out_dir, "resistance_mechanism_table.csv"), row.names = FALSE)

drug_class_table <- shortbred_rpkm %>%
  mutate(rowsums = rowSums(select(., intersect(names(.), metadata$DSUK_ID)))) %>%
  group_by(`Drug Class`) %>%
  summarise(total_count = sum(rowsums)) %>%
  mutate(total_sum = sum(total_count)) %>%
  mutate(percentage = (total_count/total_sum) * 100) %>%
  arrange(desc(percentage))

write.csv(drug_class_table, file = paste0(opt$out_dir, "drug_class_table.csv"), row.names = FALSE)

#table of aros, with name, amr gene family, total rpkms, percentage of total, #samples present in

aro_table <- shortbred_rpkm %>%
  mutate(sample_count = rowSums(select(., intersect(names(shortbred_rpkm), metadata$DSUK_ID)) != 0)) %>%
  mutate(sample_count_percentage = (sample_count / length(intersect(names(shortbred_rpkm), metadata$DSUK_ID))) * 100) %>%
  mutate(rpkm_count_sum = rowSums(select(., intersect(names(shortbred_rpkm), metadata$DSUK_ID)))) %>%
  mutate(rpkm_count_sum_percentage = (rpkm_count_sum / sum(rpkm_count_sum)) * 100) %>%
  mutate(rpkm_count_mean = rowMeans(select(., intersect(names(shortbred_rpkm), metadata$DSUK_ID)))) %>%
  select(c("ARO Name", "AMR Gene Family", "sample_count", "sample_count_percentage", "rpkm_count_sum", "rpkm_count_sum_percentage", "rpkm_count_mean"))

write.csv(aro_table, file = paste0(opt$out_dir, "aro_table.csv"), row.names = FALSE)

# Renaming ARO names to 1-x
aro_table_fd <- aro_table %>%
  arrange(desc(rpkm_count_sum)) %>%
  mutate(ARO_Index = row_number())

# Plotting
aro_fd_plot <- ggplot(aro_table_fd, aes(x = ARO_Index, y = rpkm_count_sum)) +
  geom_point() +
  labs(x = "AMR genes (1-x)", y = "RPKM Count Sum") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(size = 12),  
        axis.title.y = element_text(size = 12))

ggsave(paste0(opt$out_dir, "aro_fd_plot.png"), aro_fd_plot, height = 12, width = 16, dpi = 600)

aro_table_gene_family <- shortbred_rpkm %>% 
  group_by(`AMR Gene Family`) %>% 
  summarise(across(intersect(names(shortbred_rpkm), metadata$DSUK_ID), ~ mean(.x, na.rm = TRUE))) %>% 
  mutate(sample_count = rowSums(select(., intersect(names(shortbred_rpkm), metadata$DSUK_ID)) != 0)) %>%
  mutate(sample_count_percentage = (sample_count / length(intersect(names(shortbred_rpkm), metadata$DSUK_ID))) * 100) %>%
  mutate(rpkm_count_sum = rowSums(select(., intersect(names(shortbred_rpkm), metadata$DSUK_ID)))) %>%
  mutate(rpkm_count_sum_percentage = (rpkm_count_sum / sum(rpkm_count_sum)) * 100) %>%
  mutate(rpkm_count_mean = rowMeans(select(., intersect(names(shortbred_rpkm), metadata$DSUK_ID))))
  
write.csv(aro_table_gene_family, file = paste0(opt$out_dir, "aro_table_gene_family.csv"), row.names = FALSE)  