library(ggplot2)
library(ggpubr)
library(reshape2)
library(scales)
library(stringr)
library(dplyr)
library(tidyr)
setwd("/raid/HomeBiome/meta")

raw_stats <- read.table("snakemake_version/analysis_2023/read_stats/raw_read_stats.tsv", header = TRUE)
raw_sample_names <- str_extract(raw_stats$file, "(?<=fastq/)([A-Za-z0-9]+)_ME_L001_R([1,2])_001.fastq.gz", group = 1) 
raw_sample_read_num <- str_extract(raw_stats$file, "(?<=fastq/)([A-Za-z0-9]+)_ME_L001_R([1,2])_001.fastq.gz", group = 2)
raw_stats$file <- raw_sample_names
raw_stats$read_num <- raw_sample_read_num
raw_stats <- raw_stats %>%
  select(-c(format,type,sum_len)) %>%
  rename(raw_seqs = num_seqs, raw_min_len = min_len,raw_avg_len = avg_len, raw_max_len = max_len)

qcd_stats <- read.table("snakemake_version/analysis_2023/read_stats/bbmap_read_stats.tsv", header = TRUE)
qcd_sample_read_num <- str_extract(qcd_stats$file, "(?<=analysis_2023/human_removed/)([A-Za-z0-9]+)_clean_R(1|2)", group = 2) 
qcd_sample_names <- str_extract(qcd_stats$file, "(?<=analysis_2023/human_removed/)([A-Za-z0-9]+)_clean_R(1|2)", group = 1)
qcd_stats$file <- qcd_sample_names
qcd_stats$read_num <- qcd_sample_read_num

qcd_stats <- qcd_stats %>%
  select(-c(format,type,sum_len)) %>%
  rename(qcd_seqs = num_seqs, qcd_min_len = min_len,qcd_avg_len = avg_len, qcd_max_len = max_len)

qcd_for_summary <- qcd_stats %>% filter(read_num == 1)

kraken_stats <- read.table("snakemake_version/analysis_2023/read_stats/kraken_read_stats.tsv", header = TRUE)
kraken_sample_names <- str_extract(kraken_stats$file, "(?<=analysis_2023/kraken/)([A-Za-z0-9]+)_([a-z]+)_(1|2)", group = 1)
kraken_sample_read_num <- str_extract(kraken_stats$file, "(?<=analysis_2023/kraken/)([A-Za-z0-9]+)_([a-z]+)_(1|2)", group = 3)
kraken_sample_state <- str_extract(kraken_stats$file, "(?<=analysis_2023/kraken/)([A-Za-z0-9]+)_([a-z]+)_(1|2)", group = 2)
kraken_stats$file <- kraken_sample_names
kraken_stats$read_num <- kraken_sample_read_num
kraken_stats$kraken_state <- kraken_sample_state

kraken_stats <- kraken_stats %>%
  select(-c(format,type,sum_len)) %>%
  rename(kraken_seqs = num_seqs, kraken_min_len = min_len,kraken_avg_len = avg_len, kraken_max_len = max_len)

summary_read_stats <- raw_stats %>%
  left_join(qcd_stats, by = c("file", "read_num")) %>%
  left_join(kraken_stats, by = c("file", "read_num")) %>%
  filter(read_num == 1)

write.csv(file = "snakemake_version/analysis_2023/read_stats/summary_read_stats_df.csv", summary_read_stats)

order <- c("human_reads", "unclassified", "classified (bacterial)")

summary_read_stats <- summary_read_stats %>%
  mutate(human_reads = raw_seqs-qcd_seqs) %>%
  mutate(unclassified = qcd_seqs - kraken_seqs) %>%
  rename(`classified (bacterial)` = kraken_seqs) %>%
  select(c(human_reads,unclassified,`classified (bacterial)`,file)) %>%
  pivot_longer(
    cols = c(human_reads, unclassified, `classified (bacterial)`), 
    names_to = "type", 
    values_to = "count")

summary_read_stats$type <- factor(summary_read_stats$type, levels = order)

summary_read_plot <- summary_read_stats %>%
  ggplot(aes(x = file, y = count, fill = type)) +
  geom_bar(stat = "identity") +
  labs(title = "Reads Breakdown", y = "Number of Reads", x = "Sample", fill = "Read Status") +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("snakemake_version/analysis_2023/read_stats/read_stats_plot.png", width = 13.2, height= 9, summary_read_plot)



