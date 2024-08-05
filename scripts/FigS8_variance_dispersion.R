
# Load libraries 
library(ampvis2)
library(tidyverse)
library(vegan)

# set root dir
root.dir = rprojroot::find_rstudio_root_file()

#load data
d <- amp_load(otutable = 'data/ASV_table_combined.tsv',
              taxonomy = 'data/ASVs_combined_MiDAS5.2.R1.sintax',
              metadata = 'data/metadata.xlsx',
              pruneSingletons = T) # remove singletons

# merge replicates
d_rep <- amp_merge_replicates(d, 
                              merge_var = "SampleID_replicates", 
                              round = 'up')
d_rep$metadata <- d_rep$metadata %>%
  rename(SampleID = SampleID_replicates)

d_rep_filt_as <- amp_subset_samples(d_rep,
                                    SampleContent %in% 'AS',
                                    minreads = 10000, removeAbsentOTUs = F)


set.seed(22)
# rarefy
d_rare_as <- amp_rarefy(d_rep_filt_as, rarefy = 10000) # 10 000
dasn <- normaliseTo100(d_rare_as)
dasn <- amp_subset_samples(dasn, removeAbsentOTUs = T)

# put on groups for industrial load and size
dasn$metadata <- dasn$metadata %>% 
  mutate(Ind_text = case_when(
    `IndustrialLoad[%]` == 0 ~ "None",
    `IndustrialLoad[%]` > 0 & `IndustrialLoad[%]` <= 10 ~ "Very Low",
    `IndustrialLoad[%]` > 10 & `IndustrialLoad[%]` <= 30 ~ "Low",
    `IndustrialLoad[%]` > 30 & `IndustrialLoad[%]` <= 50 ~ "Medium",
    `IndustrialLoad[%]` > 50 & `IndustrialLoad[%]` < 100 ~ "High",
    TRUE ~ NA_character_  # For any values that don't match the above conditions
  ),
  PEdesign_text = case_when(
    PEdesign <= 25000 ~ "Small",
    PEdesign > 25000 & PEdesign <= 100000 ~ "Medium",
    PEdesign > 100000 ~ "Large",
    TRUE ~ NA_character_  # For any values that don't match the above conditions
  ))



##############################################################################
##############################################################################
################################ Functions ###################################
##############################################################################
##############################################################################

# getting distances from betadisper() object
betadisper_distances <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(distances = unlist(model$distances))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}
# getting eigenvalues out of betadisper() object
betadisper_eigenvalue <- function(model){
  temp <- data.frame(eig = unlist(model$eig))
  temp$PCoA <- row.names(temp)
  row.names(temp) <- NULL
  return(temp)
}
# getting the eigenvectors out of a betadisper() object
betadisper_eigenvector <- function(model){
  temp <- data.frame(group = model$group)
  temp2 <- data.frame(unlist(model$vectors))
  temp2$sample <- row.names(temp2)
  temp <- cbind(temp, temp2)
  temp <- dplyr::select(temp, group, sample, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}
# get centroids
betadisper_centroids <- function(model){
  temp <- data.frame(unlist(model$centroids))
  temp$group <- row.names(temp)
  temp <- dplyr::select(temp, group, dplyr::everything())
  row.names(temp) <- NULL
  return(temp)
}



get_betadisper_data <- function(model){
  temp <- list(distances = betadisper_distances(model),
               eigenvalue = betadisper_eigenvalue(model),
               eigenvector = betadisper_eigenvector(model),
               centroids = betadisper_centroids(model))
  return(temp)
}


##############################################################################
##############################################################################
##############################################################################







# make distance matrix
dist_matrix_bray <- vegdist(t(dasn$abund), method = "bray")

dasn$metadata$Region <- factor(dasn$metadata$Region, levels = c("Northern", "Central", "Southern", "Funen", "Zealand"))
dasn$metadata$Design <- factor(dasn$metadata$Design, levels = c("BioP", "ChemP", "Unknown"))
dasn$metadata$Ind_text <- factor(dasn$metadata$Ind_text, levels = c("High", "Medium", "Low",'Very Low','None'))
col_design <- c('#FDAE61','#ABD9E9','grey95')



##### plots
permdisp <- betadisper(dist_matrix_bray, dasn$metadata$PEdesign_text)
permdisp_ind <- betadisper(dist_matrix_bray, dasn$metadata$Ind_text)
#permutest(permdisp)

# get betadisper data ####
betadisper_dat <- get_betadisper_data(permdisp)
betadisper_dat_ind <- get_betadisper_data(permdisp_ind)

# add convex hull points ####
# this could be put in a function
betadisper_dat$chull <- group_by(betadisper_dat$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()

betadisper_dat_ind$chull <- group_by(betadisper_dat_ind$eigenvector, group) %>%
  do(data.frame(PCoA1 = .$PCoA1[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])],
                PCoA2 = .$PCoA2[c(chull(.$PCoA1, .$PCoA2), chull(.$PCoA1, .$PCoA2)[1])])) %>%
  data.frame()
# combine centroid and eigenvector dataframes for plotting
betadisper_lines <- merge(select(betadisper_dat$centroids, group, PCoA1, PCoA2), select(betadisper_dat$eigenvector, group, PCoA1, PCoA2), by = c('group'))
betadisper_lines_ind <- merge(select(betadisper_dat_ind$centroids, group, PCoA1, PCoA2), select(betadisper_dat_ind$eigenvector, group, PCoA1, PCoA2), by = c('group'))



##### Region
reg_pcoa <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, fill = group), shape = 21, 
             betadisper_dat$eigenvector, size = 1.8, alpha =0.8)  +
  geom_path(aes(PCoA1, PCoA2, col = group, group = group), betadisper_dat$chull,
            linewidth=0.8) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, 
                   group = row.names(betadisper_lines), col = group), betadisper_lines,
               linewidth=0.2)+
  geom_point(aes(PCoA1, PCoA2, fill = group), shape=21, betadisper_dat$centroids, size = 5)  +
  theme_bw(base_size = 12) +
  scale_fill_brewer('Region', palette="RdYlBu") +
  scale_color_brewer('Region', palette="RdYlBu") 


reg_box <- ggplot(betadisper_dat$distances, aes(group, distances, fill = group)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.5, alpha=1,
               position = position_dodge(width = 0.55)) +
  geom_point(aes(group, distances, fill = group),alpha=0.7 ,shape = 21, position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.3)) +
  theme_bw(base_size = 12)+
  ylab('Distance to centroid') +
  theme(legend.position = 'none') +
  xlab('')+
  scale_fill_brewer('Region', palette="RdYlBu") 



####### Design
des_pcoa <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, fill = group), shape = 21, 
             betadisper_dat$eigenvector, size = 1.8, alpha =0.8)  +
  geom_path(aes(PCoA1, PCoA2, col = group, group = group), betadisper_dat$chull,
            linewidth=0.8) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, 
                   group = row.names(betadisper_lines), col = group), betadisper_lines,
               linewidth=0.2)+
  geom_point(aes(PCoA1, PCoA2, fill = group), shape=21, betadisper_dat$centroids, size = 5)  +
  theme_bw(base_size = 12) +
  scale_fill_manual('Design', values = col_design) +
  scale_color_manual('Design', values = col_design ) 

des_box <- ggplot(betadisper_dat$distances, aes(group, distances, fill = group)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.5, alpha=1,
               position = position_dodge(width = 0.55)) +
  geom_point(aes(group, distances, fill = group),alpha=0.7 ,shape = 21, position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.3)) +
  theme_bw(base_size = 12)+
  ylab('Distance to centroid') +
  theme(legend.position = 'none') +
  xlab('')+
  scale_fill_manual('Design', values = col_design)


col_size <- c('#274e13','#9CA986','#C9DABF')
####### Size
size_pcoa <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, fill = group), shape = 21, 
             betadisper_dat$eigenvector, size = 1.8, alpha =0.8)  +
  geom_path(aes(PCoA1, PCoA2, col = group, group = group), betadisper_dat$chull,
            linewidth=0.8) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, 
                   group = row.names(betadisper_lines), col = group), betadisper_lines,
               linewidth=0.2)+
  geom_point(aes(PCoA1, PCoA2, fill = group), shape=21, betadisper_dat$centroids, size = 5)  +
  theme_bw(base_size = 12) +
  scale_fill_manual('Size (PE)', values = col_size) +
  scale_color_manual('Size (PE)', values = col_size ) 


size_box <- ggplot(betadisper_dat$distances, aes(group, distances, fill = group)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.5, alpha=1,
               position = position_dodge(width = 0.55)) +
  geom_point(aes(group, distances, fill = group),alpha=0.7 ,shape = 21, position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.3)) +
  theme_bw(base_size = 12)+
  ylab('Distance to centroid') +
  theme(legend.position = 'none') +
  xlab('')+
  scale_fill_manual('Size (PE)', values = col_size)




col_indu <- c('#542b09','darkred','#b48662','#faa663','#fddbc0')
####### Size
indu_pcoa <- ggplot() +
  geom_point(aes(PCoA1, PCoA2, fill = group), shape = 21, 
             betadisper_dat_ind$eigenvector, size = 1.8, alpha =0.8)  +
  geom_path(aes(PCoA1, PCoA2, col = group, group = group), betadisper_dat_ind$chull,
            linewidth=0.8) +
  geom_segment(aes(x = PCoA1.x, y = PCoA2.x, yend = PCoA2.y, xend = PCoA1.y, 
                   group = row.names(betadisper_lines_ind), col = group), betadisper_lines_ind,
               linewidth=0.2)+
  geom_point(aes(PCoA1, PCoA2, fill = group), shape=21, betadisper_dat_ind$centroids, size = 5)  +
  theme_bw(base_size = 12) +
  scale_fill_manual('Industrial\nload (%)', values = col_indu) +
  scale_color_manual('Industrial\nload (%)', values = col_indu ) 


indu_box <- ggplot(betadisper_dat_ind$distances, aes(group, distances, fill = group)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.5, alpha=1,
               position = position_dodge(width = 0.55)) +
  geom_point(aes(group, distances, fill = group),alpha=0.7 ,shape = 21, position = position_jitterdodge(dodge.width = 0.55, jitter.width = 0.3)) +
  theme_bw(base_size = 12)+
  ylab('Distance to centroid') +
  theme(legend.position = 'none') +
  xlab('')+
  scale_fill_manual('Industrial load (%)', values = col_indu)


# save plot
p_tot <- (reg_pcoa +reg_box)/
          (des_pcoa + des_box) /
          (size_pcoa +size_box)/
        (indu_pcoa +indu_box)

ggsave('output/FigS8_dispersion.png',
       p_tot, width=9, height=13,
       dpi = 320)

ggsave(filename="output/FigS8_dispersion.pdf", 
       p_tot, width=9, height=13, useDingbats=FALSE, limitsize=FALSE)

