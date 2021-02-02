library(tidyverse)

setwd('/Users/alicia/martin_lab/projects/neurogap/data/high_coverage/pca')

ref <- read.table(gzfile('hgdp_tgp_neurogap_pca_scores.txt.bgz'), header=T)
proj <- read.table(gzfile('hgdp_tgp_neurogap_pca_projected_scores.txt.bgz'), header=T)
hgdp_tgp <- read.delim('/Users/alicia/martin_lab/projects/hgdp_tgp/gnomad_meta_v1.txt', header=T, sep='\t') %>%
  dplyr::select(s, project_meta.title, starts_with('hgdp_tgp_meta')) 
pca <- bind_rows(ref, proj) %>%
  left_join(hgdp_tgp)
pca[is.na(pca[,'hgdp_tgp_meta.Genetic.region']),'hgdp_tgp_meta.Genetic.region'] <- 'AFR'
pca[is.na(pca[,'project_meta.title']),'project_meta.title'] <- 'NeuroGAP'

ggplot(pca, aes(x=PC1, y=PC2, color=hgdp_tgp_meta.Genetic.region, shape=project_meta.title)) +
  geom_point() +
  theme_classic()
