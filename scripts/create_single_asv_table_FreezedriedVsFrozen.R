#!/usr/bin/env Rscript

# create_single_asv_table_FreezedriedVsFrozen.R
#
# combine reada samples sequences twuce and add relative abundances 
#
# Author: daniel.lundin@lnu.se hanna.bensch@lnu.se

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(vegan))

# read in asvtable with frozen and freeze-died samples and negative controls
asvsFD <- read_tsv("asvs_FreezedriedVsFrozen.tsv", 
                   col_types = cols(.default = col_character(),
                                    count = col_integer()))


  
# Filter away negative control ASVs and calculate relative abundances
asvsFDfilt <- asvsFD %>%
    anti_join(
      asvs %>% filter(grepl('neg', sample)) %>% distinct(asv),
      by = c('asv')
     )
  
# Write out table summed over sample and asv (use when confident that all is OK)
asvsFDfilt %>% group_by(asv, sample) %>% summarise(count = sum(count), .groups = 'drop') %>%
    group_by(sample) %>% mutate(relab = count/sum(count)) %>% ungroup() %>%
    write_tsv("asv_summed_FreezedriedVsFrozen.tsv")
  