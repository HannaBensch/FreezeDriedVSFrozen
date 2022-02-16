# create permanova table freeze-dried vs frozen study
#Author: Hanna Bensch

library(tidyverse)
library(phyloseq)
library(vegan)

# load phyloseq object created in Freezedried_vs_min80_Unifraq_betadiversity.Rmd used for calculating alpha diversity and unifrac distances
pseq <- readRDS( "../data/pseq.rds")

wuni_dist2 <- phyloseq::distance(subset_samples(pseq, SampleNumber != "4"), "wunifrac")
uni_dist2 <- phyloseq::distance(subset_samples(pseq, SampleNumber != "4"), "unifrac")

#######################
# *On un-weighted unifrac distance excluding sample 4*
# with plate number as strata argument
perm <- how(nperm = 9999)
set.seed(3000)
setBlocks(perm) <- with(as.tibble(pseq@sam_data) %>% filter( SampleNumber != "4"), Plate_No)
set.seed(3000)
permanova <- adonis2(uni_dist2 ~ Treatment + NbReads + SampleNumber,data = as(sample_data(subset_samples(pseq, SampleNumber != "4")), "data.frame"), permutations = perm,  by = "margin")
print(as.data.frame(permanova)) # Treatment significant but just small proportion
t1<-print(as.data.frame(permanova)) %>% tibble::rownames_to_column("Factor") %>% mutate(Dissimilarity_matrix ="Unweighted UniFrac", Permanova = 1) %>% rename(p = "Pr(>F)") %>% 
  filter(Factor %in% c("Treatment", "NbReads", "SampleNumber")) %>%
  select(Permanova, Dissimilarity_matrix, Factor, `F`, R2, p)


# treatment betadisp
set.seed(3000)  
perm <- how(nperm = 9999)
betadisp <- betadisper(uni_dist2, as(sample_data(subset_samples(pseq, SampleNumber != "4")), "data.frame")$Treatment, 
                       type = "centroid", bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
t1$p_betadisp[t1$Factor == "Treatment"]<- print(as.data.frame(anova(betadisp))) %>% rename(p = "Pr(>F)")  %>% select(p) %>% slice(1) %>% pull()

# plate no
perm <- how(nperm = 9999)
set.seed(3000)  
permanova <- adonis2(uni_dist2 ~ Plate_No, data = as(sample_data(subset_samples(pseq, SampleNumber != "4")), "data.frame"), permutations = perm,  by = "margin")
t2 <- print(as.data.frame(permanova)) %>% tibble::rownames_to_column("Factor") %>% mutate(Dissimilarity_matrix ="Unweighted UniFrac", Permanova = 2)  %>% rename(p = "Pr(>F)") %>% slice(1) %>%
  select(Permanova, Dissimilarity_matrix, Factor, `F`, R2, p)
set.seed(3000)  
betadisp <- betadisper(uni_dist2, as(sample_data(subset_samples(pseq, SampleNumber != "4")), "data.frame")$Plate_No, 
                       type = "centroid", bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
t2$p_betadisp <- print(as.data.frame(anova(betadisp))) %>% rename(p = "Pr(>F)") %>% select(p) %>% slice(1) %>% pull()

##################
# *On weighted unifrac distance excluding sample 4*
# with plate number as strata argument
perm <- how(nperm = 9999)
set.seed(3000)
setBlocks(perm) <- with(as.tibble(pseq@sam_data) %>% filter( SampleNumber != "4"), Plate_No)
set.seed(3000)
permanova <- adonis2(wuni_dist2 ~ Treatment + NbReads + SampleNumber,data = as(sample_data(subset_samples(pseq, SampleNumber != "4")), "data.frame"), permutations = perm,  by = "margin")
print(as.data.frame(permanova)) # Treatment significant but just small proportion
t3<-print(as.data.frame(permanova)) %>% tibble::rownames_to_column("Factor") %>% mutate(Dissimilarity_matrix ="Weighted UniFrac", Permanova = 3) %>% rename(p = "Pr(>F)") %>% 
  filter(Factor %in% c("Treatment", "NbReads", "SampleNumber")) %>%
  select(Permanova, Dissimilarity_matrix, Factor, `F`, R2, p)


# treatment betadisp
set.seed(3000)  
perm <- how(nperm = 9999)
betadisp <- betadisper(wuni_dist2, as(sample_data(subset_samples(pseq, SampleNumber != "4")), "data.frame")$Treatment, 
                       type = "centroid", bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
t3$p_betadisp[t3$Factor == "Treatment"]<- print(as.data.frame(anova(betadisp))) %>% rename(p = "Pr(>F)")  %>% select(p) %>% slice(1) %>% pull()


# plate no
perm <- how(nperm = 9999)
set.seed(3000)  
permanova <- adonis2(wuni_dist2 ~ Plate_No, data = as(sample_data(subset_samples(pseq, SampleNumber != "4")), "data.frame"), permutations = perm,  by = "margin")
t4 <- print(as.data.frame(permanova)) %>% tibble::rownames_to_column("Factor") %>% mutate(Dissimilarity_matrix ="Weighted UniFrac", Permanova = 4)  %>% rename(p = "Pr(>F)") %>% slice(1) %>%
  select(Permanova, Dissimilarity_matrix, Factor, `F`, R2, p)
set.seed(3000)  
betadisp <- betadisper(uni_dist2, as(sample_data(subset_samples(pseq, SampleNumber != "4")), "data.frame")$Plate_No, 
                       type = "centroid", bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
#"boxplot"(betadisp, ylab = "Distance to centroid", xlab = "Treatment")  
t4$p_betadisp <- print(as.data.frame(anova(betadisp))) %>% rename(p = "Pr(>F)") %>% select(p) %>% slice(1) %>% pull()


############### 
# Hellinger 
asvs <- data.frame(otu_table(pseq)) %>% rownames_to_column("asv") %>%
  pivot_longer(starts_with("A_"), names_to = "Asample", values_to = "count") %>% filter(count > 0)
hell <- asvs %>%
  select(asv, Asample, count) %>% arrange(Asample) %>%
  pivot_wider(names_from = 'asv', values_from = 'count', values_fill = 0) %>%
  as.data.frame() %>% tibble::column_to_rownames('Asample') %>%
  decostand(method = 'hellinger') %>%
  data.frame() %>% dist(method = "euclidian")
metadf<- as(sample_data(pseq), "data.frame") %>% rownames_to_column("Asample") %>%  arrange(Asample)

# with plate number as strata argument
perm <- how(nperm = 9999)
set.seed(3000)
setBlocks(perm) <- with(metadf, Plate_No)
set.seed(3000)
permanova <- adonis2(hell ~ Treatment +NbReads + SampleNumber, data = metadf, permutations = perm,  by = "margin") 
t5<-print(as.data.frame(permanova)) %>% tibble::rownames_to_column("Factor") %>% mutate(Dissimilarity_matrix ="Hellinger", Permanova = 5) %>% rename(p = "Pr(>F)") %>% 
  filter(Factor %in% c("Treatment", "NbReads", "SampleNumber")) %>%
  select(Permanova, Dissimilarity_matrix, Factor, `F`, R2, p)

# treatment  betadisp
set.seed(3000)  
perm <- how(nperm = 9999)
betadisp <- betadisper(hell, metadf$Treatment, 
                       type = "centroid", bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
#"boxplot"(betadisp, ylab = "Distance to centroid", xlab = "Treatment")  
t5$p_betadisp[t5$Factor == "Treatment"]<- print(as.data.frame(anova(betadisp))) %>% rename(p = "Pr(>F)")  %>% select(p) %>% slice(1) %>% pull()

# plate no
perm <- how(nperm = 9999)
set.seed(3000)  
permanova <- adonis2(hell ~ Plate_No, data = metadf, permutations = perm,  by = "margin")
t6 <- print(as.data.frame(permanova)) %>% tibble::rownames_to_column("Factor") %>% mutate(Dissimilarity_matrix ="Hellinger", Permanova = 6)  %>% rename(p = "Pr(>F)") %>% slice(1) %>%
  select(Permanova, Dissimilarity_matrix, Factor, `F`, R2, p)
set.seed(3000)  
betadisp <- betadisper(hell, metdf$Plate_No, 
                       type = "centroid", bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
#"boxplot"(betadisp, ylab = "Distance to centroid", xlab = "Treatment")  
t6$p_betadisp <- print(as.data.frame(anova(betadisp))) %>% rename(p = "Pr(>F)") %>% select(p) %>% slice(1) %>% pull()


###########
# CLR

m1 <-  asvs  %>% select(asv,Asample) %>%
  complete(asv, Asample) %>% left_join(asvs %>% select(asv, count, Asample)) %>%
  mutate( count = replace_na(count, 0),
          count = count +0.65) %>% 
  pivot_wider(names_from = 'asv', values_from='count') %>% 
  tibble::column_to_rownames('Asample') %>% as.matrix()
clr <- t(apply(t(m1), 2, compositions::clr))
clr<- clr %>% dist( method = "euclidian") 
metadf <- metadf %>% tibble::column_to_rownames('Asample')
str(metadf)



# with plate number as strata argument
perm <- how(nperm = 9999)
set.seed(3000)
setBlocks(perm) <- with(metadf, Plate_No)
permanova <- adonis2(clr ~ Treatment +NbReads +SampleNumber, data = metadf, permutations = perm,  by = "margin") 
t7<-print(as.data.frame(permanova)) %>% tibble::rownames_to_column("Factor") %>% mutate(Dissimilarity_matrix ="CLR", Permanova = 7) %>% rename(p = "Pr(>F)") %>% 
  filter(Factor %in% c("Treatment", "NbReads", "SampleNumber")) %>%
  select(Permanova, Dissimilarity_matrix, Factor, `F`, R2, p)

# treatment betadsip
set.seed(3000)  
perm <- how(nperm = 9999)
betadisp <- betadisper(clr, metadf$Treatment, 
                       type = "centroid", bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
t7$p_betadisp[t7$Factor == "Treatment"]<- print(as.data.frame(anova(betadisp))) %>% rename(p = "Pr(>F)")  %>% select(p) %>% slice(1) %>% pull()

# plate no
perm <- how(nperm = 9999)
set.seed(3000)  
permanova <- adonis2(clr ~ Plate_No, data = metadf, permutations = perm,  by = "margin")
t8 <- print(as.data.frame(permanova)) %>% tibble::rownames_to_column("Factor") %>% mutate(Dissimilarity_matrix ="CLR", Permanova = 8)  %>% rename(p = "Pr(>F)") %>% slice(1) %>%
  select(Permanova, Dissimilarity_matrix, Factor, `F`, R2, p)
set.seed(3000)  
betadisp <- betadisper(clr, metadf$Plate_No, 
                       type = "centroid", bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
t8$p_betadisp <- print(as.data.frame(anova(betadisp))) %>% rename(p = "Pr(>F)") %>% select(p) %>% slice(1) %>% pull()



# combine tables
combinedt <- rbind(t1, t2, t3, t4, t5, t6, t7, t8)
combinedt <- combinedt  %>% mutate_if(is.numeric,
                     round,
                     digits = 3)
combinedt %>% mutate(Factor = case_when(Factor == "NbReads" ~ "Library size",
                                        Factor == "SampleNumber" ~ "Sample identity",
                                        Factor == "Plate_No" ~ "Plate number",
                                        TRUE ~ Factor),
                     p_betadisp = as.character(p_betadisp),
                     p_betadisp = case_when(p_betadisp == 0 ~ "<0.001",
                                            is.na(p_betadisp) ~ "", 
                                              TRUE ~ p_betadisp),
                     p = as.character(p),
                     p = case_when(p == "0" ~ "<0.001",
                                            TRUE ~ p)) %>%
  arrange(Permanova, Factor) #%>% write_csv("data/permanova_table.csv")

