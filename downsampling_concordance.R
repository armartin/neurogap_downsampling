library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(gtable)
library(maptools)
library(mapdata)
library(raster)
library(ggmap)

setwd('/Users/alicia/martin_lab/projects/neurogap/data/high_coverage/stats')

sample_id_map <- read.table('../sample_id_map.txt', header=T)
depth <- read.table('../downsample-bam.txt') %>%
  separate(V1, into=c(NA, NA, NA, NA, NA, 'sample', NA, NA), sep='/') %>%
  left_join(sample_id_map, by=c('sample'='high_cov_id')) %>%
  filter(!sample %in% c('NGE0018', 'NGE0130'))

depth %>%
  group_by(site) %>%
  summarize(mean_cov=mean(V2), sd_cov=sd(V2), min_cov=min(V2), max_cov=max(V2))

read_file <- function(coverage, refined=FALSE, imputed=FALSE, array=FALSE, topmed=FALSE, gencove=FALSE, gencove_compare=FALSE, imputed_compare=FALSE) {
  if(refined) {
    filename <- paste0('neurogap_', coverage, 'X.refined_concordance_variants.tsv')
    process <- 'Refined'
  } else if (imputed) {
    filename <- paste0('neurogap_', coverage, 'X.imputed_concordance_variants.tsv')
    process <- 'Imputed'
  } else if (array) {
    filename <- paste0('NeuroGap_30x_', coverage, '.imputed_concordance_variants.tsv')
    process <- 'Imputed'
  } else if (topmed) {
    filename <- paste0('NeuroGap_30x_', coverage, '_topmed.dose_concordance_variants.tsv')
    process <- 'TOPMed'
  } else if (gencove) {
    filename <- paste0('../gencove/stats/merge_', coverage, '_concordance_variants.tsv')
    process <- 'Gencove'
  } else if (gencove_compare) {
    filename <- paste0('../gencove/stats/side_by_side/merge_intersect_', coverage, '_concordance_variants.tsv')
    process <- 'Gencove'
  } else if (imputed_compare) {
    filename <- paste0('../gencove/stats/side_by_side/neurogap_intersect_', coverage, 'X.imputed_concordance_variants.tsv')
    process <- 'BEAGLE'
  } else {
    filename <- paste0('neurogap_', coverage, 'X_concordance_variants.tsv')
    process <- 'Raw'
  }
  concordance <- read.table(filename, header=T)
  concordance <- concordance %>% 
    mutate(cov=coverage,
           is_snp=case_when(snp=='true' ~ 'SNP',
                            TRUE ~ 'INDEL'),
           process=process)
  return(concordance)
}

plot_concordance <- function(concordance) {
  p1 <- ggplot(concordance, aes(x=freq, y=non_ref_concordance, color=cov, size=n_variants)) +
    facet_grid(~is_snp) +
    geom_point() +
    #geom_line() +
    labs(x='Frequency', y='Non-reference concordance') +
    scale_color_brewer(palette='Dark2', name='Depth') +
    ylim(0,1) +
    guides(size=FALSE) +
    theme_classic() +
    theme(text = element_text(size=16))
  return(p1)
}

style_sheet <- read.delim('style_sheet.txt', header=T)
covs <- c('0.5', '1.0', '2.0', '4.0', '6.0', '10.0', '20.0')
lowcovs <- c('0.5', '1.0', '2.0', '4.0', '6.0')
arrays <- c('H3Africa', 'Omni2.5', 'MEGA', 'PsychChip', 'GSA')
#covs <- c('0.5', '1.0', '2.0')

raw_concordance <- map_dfr(covs, read_file)
refined_concordance <- map_dfr(lowcovs, function(x) {read_file(x, refined=TRUE)} )
imputed_concordance <- map_dfr(lowcovs, function(x) {read_file(x, imputed=TRUE)} )
gencove_concordance <- map_dfr(lowcovs, function(x) {read_file(x, gencove=TRUE)})
array_concordance <- map_dfr(arrays, function(x) {read_file(x, array=TRUE)} )
topmed_concordance <- map_dfr(arrays, function(x) {read_file(x, topmed=TRUE)} )
gencove_compare <- map_dfr(lowcovs, function(x) {read_file(x, gencove_compare=TRUE)})
imputed_compare <- map_dfr(lowcovs, function(x) {read_file(x, imputed_compare=TRUE)})
raw_concordance$cov <- factor(raw_concordance$cov, levels=c(rev(covs)))

all_concordance <- bind_rows(raw_concordance, refined_concordance, imputed_concordance, array_concordance) %>%
  left_join(style_sheet, by=c('cov'='Technology'))
all_concordance$process <- factor(all_concordance$process, levels=c('Raw', 'Refined', 'Imputed'))
all_concordance$cov <- factor(all_concordance$cov, levels=c(rev(covs), arrays))
snp_concordance <- subset(all_concordance, is_snp=='SNP')

side_by_side <- bind_rows(gencove_compare, imputed_compare) %>%
  left_join(style_sheet, by=c('cov'='Technology')) %>%
  filter(is_snp=='SNP')
side_by_side$process <- factor(side_by_side$process, levels=c('Gencove', 'BEAGLE'))
side_by_side$cov <- factor(side_by_side$cov, levels=rev(covs))

gencove_chip <- bind_rows(gencove_compare, array_concordance) %>%
  left_join(style_sheet, by=c('cov'='Technology')) %>%
  filter(is_snp=='SNP')
gencove_chip$process <- factor(gencove_chip$process, levels=c('Gencove', 'Imputed'))
gencove_chip$cov <- factor(gencove_chip$cov, levels=c(rev(covs), arrays))


snp_concordance <- subset(all_concordance, is_snp=='SNP')

color_vec <- as.character(style_sheet$Colors)
names(color_vec) <- style_sheet$Technology
shape_vec <- c(rep(21, 7), rep(4, 5))
names(shape_vec) <- style_sheet$Technology

p1 <- ggplot(subset(raw_concordance, freq!=0), aes(x=freq, y=non_ref_concordance, color=cov, fill=cov, size=n_variants, shape=cov)) +
  facet_grid(~is_snp) +
  geom_point() +
  labs(x='Frequency', y='Non-reference concordance') +
  scale_color_manual(name='Depth', values=color_vec) +
  scale_fill_manual(name='Depth', values=color_vec) +
  scale_shape_manual(name='Depth', values=shape_vec) +
  scale_y_continuous(minor_breaks = seq(0 , 1, length.out=11), breaks = seq(0, 1, length.out=6), limits=c(0,1)) +
  guides(size=FALSE) +
  theme_bw() +
  theme(text = element_text(size=16))

ggsave('raw_downsampling.pdf', p1, width=10, height=5)

p2 <- ggplot(snp_concordance, aes(x=freq, y=non_ref_concordance, color=cov, fill=cov, size=n_variants, shape=cov)) +
  facet_grid(~process) +
  geom_point() +
  labs(x='Frequency', y='Non-reference concordance') +
  scale_color_manual(name='Depth/Array', values=color_vec) +
  scale_fill_manual(name='Depth/Array', values=color_vec) +
  scale_shape_manual(name='Depth/Array', values=shape_vec) +
  scale_y_continuous(minor_breaks = seq(0 , 1, length.out=11), breaks = seq(0, 1, length.out=6), limits=c(0,1)) +
  guides(size=FALSE) +
  theme_bw() +
  theme(text = element_text(size=14))

ggsave('snp_downsampling.pdf', p2, width=14, height=5)

p2_1 <- ggplot(side_by_side, aes(x=freq, y=non_ref_concordance, color=cov)) +
  geom_line(aes(linetype=process)) +
  labs(x='Frequency', y='Non-reference concordance', title='Imputation methods comparison') +
  scale_color_manual(name='Depth', values=color_vec) +
  scale_linetype_manual(name='Imputation\nmethod', values=c("solid", "dashed")) +
  scale_y_continuous(breaks = seq(0, 1, by = .2), limits=c(0,1)) +
  #ylim(0,1) +
  theme_bw() +
  theme(text = element_text(size=14))

p2_3 <- ggplot(gencove_chip, aes(x=freq, y=non_ref_concordance, color=cov, fill=cov)) +
  geom_point(data=subset(gencove_chip, process=='Imputed'), aes(size=n_variants, shape=cov)) +
  geom_line(data=subset(gencove_chip, process=='Gencove'), aes(linetype=process)) +
  labs(x='Frequency', y='Non-reference concordance', title='Gencove vs chip comparison') +
  scale_color_manual(name='Depth/Array', values=color_vec) +
  scale_shape_manual(name='Depth/Array', values=shape_vec) +
  scale_y_continuous(breaks = seq(0, 1, by = .2), limits=c(0,1)) +
  guides(size=FALSE, linetype=FALSE, shape=FALSE, fill=FALSE) +
  theme_bw() +
  theme(text = element_text(size=14))
#fill=cov, size=n_variants, linetype=process

bottom_row <- plot_grid(p2_1, p2_3, labels=c('B', 'C'))
p2_complete <- plot_grid(p2, bottom_row, labels = c('A', ''), ncol = 1)
ggsave('snp_downsampling_compare.eps', p2_complete, width=10, height=8)

p2_5 <- ggplot(subset(snp_concordance, process=='Imputed'), aes(x=freq, y=non_ref_concordance, color=cov, fill=cov, size=n_variants, shape=cov)) +
  geom_point() +
  labs(x='Frequency', y='Non-reference concordance') +
  scale_color_manual(name='Depth/Array', values=color_vec) +
  scale_fill_manual(name='Depth/Array', values=color_vec) +
  scale_shape_manual(name='Depth/Array', values=shape_vec) +
  scale_y_continuous(minor_breaks = seq(0 , 1, length.out=11), breaks = seq(0, 1, length.out=6), limits=c(0,1)) +
  guides(size=FALSE) +
  theme_bw() +
  theme(text = element_text(size=16))

ggsave('snp_downsampling_imputed_only.pdf', p2_5, width=6, height=5)


p2_7 <- ggplot(subset(snp_concordance, process=='Imputed' & cov %in% c('2.0', '4.0', '6.0', 'GSA')), aes(x=freq, y=non_ref_concordance, color=cov, fill=cov, size=n_variants, shape=cov)) +
  geom_point() +
  labs(x='Frequency', y='Non-reference concordance') +
  scale_color_brewer(name='Depth/Array', palette='Set1') +
  scale_fill_brewer(name='Depth/Array', palette='Set1') +
  scale_shape_manual(name='Depth/Array', values=shape_vec) +
  scale_y_continuous(minor_breaks = seq(0 , 1, length.out=11), breaks = seq(0, 1, length.out=6), limits=c(0,1)) +
  guides(size=FALSE) +
  theme_bw() +
  theme(text = element_text(size=16),
        legend.position = c(.95,.05),
        legend.justification = c(1,0))

ggsave('snp_downsampling_imputed_only_grant.pdf', p2_7, width=5, height=5)

p2_8 <- ggplot(subset(snp_concordance, process=='Imputed' & cov %in% c('2.0', '4.0', '6.0', 'GSA', 'Omni2.5')), aes(x=freq, y=non_ref_concordance, color=cov, fill=cov, size=n_variants, shape=cov)) +
  geom_point() +
  labs(x='Frequency', y='Non-reference concordance') +
  scale_color_brewer(name='Depth/Array', palette='Set1') +
  scale_fill_brewer(name='Depth/Array', palette='Set1') +
  scale_shape_manual(name='Depth/Array', values=shape_vec) +
  scale_y_continuous(minor_breaks = seq(0 , 1, length.out=11), breaks = seq(0, 1, length.out=6), limits=c(0,1)) +
  guides(size=FALSE) +
  theme_bw() +
  theme(text = element_text(size=16),
        legend.position = c(.95,.05),
        legend.justification = c(1,0))

ggsave('snp_downsampling_imputed_only_hailiang.pdf', p2_8, width=4, height=4)

# Concordance by population -----------------------------------------------

read_file_site <- function(coverage, site, imputed=FALSE, gencove=FALSE, array=FALSE) {
  if (imputed) {
    filename <- paste0('neurogap_', coverage, 'X.imputed_', site, '_concordance_variants.tsv')
    process <- 'Imputed'
  } else if (gencove) {
    filename <- paste0('../gencove/stats/merge_', coverage, '_', site, '_concordance_variants.tsv')
  } else { # array
    filename <- paste0('NeuroGap_30x_', coverage, '.imputed_', site, '_concordance_variants.tsv')
    process <- 'Imputed'
  }
  concordance <- read.table(filename, header=T)
  concordance <- concordance %>% 
    mutate(cov=coverage,
           is_snp=case_when(snp=='true' ~ 'SNP',
                            TRUE ~ 'INDEL'),
           process=process,
           site=site)
  concordance <- concordance %>%
    group_by(cov, process, site) %>%
    mutate(snp_count=sum(total_concordant, total_discordant)) %>%
    mutate(snp_concordance=sum(total_concordant) / snp_count)
  return(concordance)
}

sites <- c('AddisEthiopia','CapeTownSouthAfrica','KEMRIKenya','MakerereUganda','MoiKenya')
sites_gencove <- c('AddisEthiopia','CapeTownSouthAfrica','MakerereUganda','MoiKenya')

imputed_concordance_site <- map_dfr(lowcovs, function(x) {map_dfr(sites, function(y) {read_file_site(x, y, imputed=TRUE)} ) } )
gencove_concordance_site <- map_dfr(lowcovs, function(x) {map_dfr(sites_gencove, function(y) {read_file_site(x, y, gencove=TRUE)} ) } )
array_concordance_site <- map_dfr(arrays, function(x) {map_dfr(sites, function(y) {read_file_site(x, y, array=TRUE)} ) } )

gencove_concordance_site %>%
  ungroup() %>%
  dplyr::select(cov, site, snp_concordance) %>%
  unique() %>%
  spread(site, snp_concordance)

imputed_concordance_site$cov <- factor(imputed_concordance_site$cov, levels=c(rev(covs)))
snp_concordance_site <- subset(imputed_concordance_site, is_snp=='SNP')

all_concordance_sites <- bind_rows(imputed_concordance_site, array_concordance_site) %>%
  left_join(style_sheet, by=c('cov'='Technology'))
all_concordance_sites$cov <- factor(all_concordance_sites$cov, levels=c(rev(covs), arrays))
snp_concordance_sites <- subset(all_concordance_sites, is_snp=='SNP')

p3 <- ggplot(snp_concordance_sites, aes(x=freq, y=non_ref_concordance, color=cov, fill=cov, size=n_variants, shape=cov)) +
  facet_grid(~site) +
  geom_point() +
  labs(x='Frequency', y='Non-reference concordance') +
  scale_color_manual(name='Depth/Array', values=color_vec) +
  scale_fill_manual(name='Depth/Array', values=color_vec) +
  scale_shape_manual(name='Depth/Array', values=shape_vec) +
  scale_y_continuous(minor_breaks = seq(0 , 1, length.out=11), breaks = seq(0, 1, length.out=6), limits=c(0,1)) +
  guides(size=FALSE) +
  theme_bw() +
  theme(text = element_text(size=16))

ggsave('site_concordance.pdf', p3, width=12, height=5)


# Downsampling ancestry ethnicity language --------------------------------

gencove_anc <- read.delim('../gencove/ancestry/gencove_metadata_anc.txt', header=T, sep='\t') 

master_phenos <- read.csv('/Users/alicia/daly_lab/neurogap/data/phenotypes/ngpsych_FROZEN_08apr2019.csv', header=T)

gencove_anc %>%
  left_join(master_phenos, by=c('sample'='collaborator_participant_id'))
  
gencove_anc_grouped <- gencove_anc %>%
  gather('ancestry', 'ancestry_fraction', CAFRICA:WAFRICA)

p_anc <- ggplot(gencove_anc_grouped, aes(x=ancestry_fraction, fill=ancestry)) +
  facet_grid(~site) +
  geom_histogram(alpha=0.8, position='identity') +
  labs(x='Gencove ancestry fraction', y='Count') +
  theme_bw()

ggsave('gencove_anc.pdf', p_anc, width=12, height=4)

primary_eth <- count(gencove_anc %>% filter(!sample %in% c('NGE0018', 'NGE0130')), site, ethnicity_1)
write.table(primary_eth, 'primary_eth.txt', quote=F, row.names=F, sep='\t')


# Microbiome analyses -----------------------------------------------------

read_kraken <- function(kraken_sample, kraken_path) {
  kraken_file <- read.delim(gzfile(as.character(kraken_path)), header=F) %>%
    mutate(sample_id=kraken_sample)
  return(kraken_file)
}


kraken <- read.table('../gencove/microbiome/kraken_reports.txt', header=T)

kraken_outputs <- map2_dfr(kraken$sample_id, kraken$path, read_kraken) %>%
  left_join(sample_id_map, by=c('sample_id'='high_cov_id'))
colnames(kraken_outputs)[1:6] <- c('read_percentage', 'num_reads_covered', 'num_reads_assigned', 'rank', 'ncbi_taxonomy_id', 'scientific_name')
kraken_outputs <- kraken_outputs %>%
  group_by(sample_id, rank) %>%
  mutate(relative_abundance=read_percentage/sum(read_percentage),
         sci_name = str_trim(scientific_name))

kraken_summary <- kraken_outputs %>%
  group_by(rank, sci_name) %>%
  summarize(mean_vals = mean(relative_abundance)) %>%
  arrange(desc(mean_vals))

phyla <- subset(kraken_outputs, rank=='P'&relative_abundance>0)
phyla$sci_name <- factor(phyla$sci_name, levels=rev(subset(kraken_summary, rank=='P')$sci_name))
ind_order <- subset(phyla, sci_name==subset(kraken_summary, rank=='P')$sci_name[1]) %>%
  arrange(relative_abundance)
phyla$sample_id <- factor(phyla$sample_id, levels=ind_order$sample_id)

color_vec <- colorRampPalette(brewer.pal(9, 'Set1'))(length(unique(phyla$scientific_name)))
p_phyla <- ggplot(phyla, aes(x=sample_id, y=relative_abundance, fill=sci_name)) +
  facet_wrap(~site, scales='free_x', nrow=1) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=rev(color_vec), name='Phylum') +
  theme_classic() +
  labs(y='Relative abundance', x='Sample') +
  theme(axis.text.x = element_text(angle=90))

order <- subset(kraken_outputs, rank=='O'&relative_abundance>0)
order$sci_name <- factor(order$sci_name, levels=rev(subset(kraken_summary, rank=='O')$sci_name))
ind_order <- subset(order, sci_name==subset(kraken_summary, rank=='O')$sci_name[1]) %>%
  arrange(relative_abundance)
order$sample_id <- factor(order$sample_id, levels=ind_order$sample_id)

color_vec <- colorRampPalette(brewer.pal(9, 'Set1'))(length(unique(order$scientific_name)))
p_order <- ggplot(order, aes(x=sample_id, y=relative_abundance, fill=sci_name)) +
  facet_wrap(~site, scales='free_x', nrow=1) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=rev(color_vec), name='Order') +
  theme_classic() +
  labs(y='Relative abundance', x='Sample') +
  guides(fill=guide_legend(ncol=4)) +
  theme(axis.text.x = element_text(angle=90))

ggsave('microbiome_order.pdf', p_order, width=12, height=4)


# PCA and Africa map ------------------------------------------------------

data(wrld_simpl)
world <- fortify(wrld_simpl)

# latitude/longitude bounds for africa
xlim = c(-20,50)
ylim = c(-35,35)

lims = SpatialPoints(coords = data_frame(x = xlim, y = ylim), proj4string = CRS("+proj=longlat +datum=WGS84"))%>%
  spTransform(CRS("+init=epsg:3857"))

# Read the latitude/longitdue data for populations in AGVP/TGP
pop_pos <- read.csv('../african_pop_neurogap.csv', header=T)#, stringsAsFactors = F)

color_vec <- as.character(pop_pos$Color)
names(color_vec) <- pop_pos$Population
shape_vec <- pop_pos$Shape
names(shape_vec) <- pop_pos$Population

neurogap_country_codes <- c('ZAF', 'ETH', 'KEN', 'UGA')
neurogap_colors <- data.frame(id=unique(world$id),
                              color='lightyellow')
neurogap_colors$color <- as.character(neurogap_colors$color)
neurogap_colors$color[neurogap_colors$id %in% neurogap_country_codes] <- 'snow2'
world <- world %>% left_join(neurogap_colors, by='id')

# plot the map of Africa with data points labeled
p_map <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group=group, fill=color), color='lightgrey') +
  scale_fill_identity() +
  geom_point(data = pop_pos, aes(Longitude, Latitude, color=Population), shape=19, size=3) +
  coord_fixed(xlim = c(-20,50), ylim = c(-35,35)) +
  labs(x='Longitude', y='Latitude') +
  theme_classic() +
  scale_color_brewer(palette='Set1', name = "Site") +
  #scale_shape_manual(name = "Population", values = shape_vec) +
  theme(panel.background = element_rect(fill = "lightblue"),
        plot.background = element_rect(fill = "transparent", color = NA),
        #legend.position='bottom',
        text = element_text(size=16),
        axis.text = element_text(color='black'))
p_map <- p_map + guides(fill=F, color=F)

pca <- read.table(gzfile('../pca/NeuroGap_30x_Pilot_Callset_lightfilt_maf05scores.txt.bgz'), header=T)
sample_map <- read.table('../sample_id_map.txt', header=T)
pca <- pca %>%
  left_join(sample_map, by=c('s'='high_cov_id')) %>%
  filter(!s %in% c('NGE0018', 'NGE0130'))

p_pca <- ggplot(pca, aes(x=PC1, y=PC2, color=site)) +
  geom_point(size=3) +
  theme_classic() +
  scale_color_brewer(palette='Set1', name='Site') +
  #guides(color=F) +
  theme(text = element_text(size=16),
        axis.text = element_text(color='black'),
        legend.position = 'left')

p_map_pca <- plot_grid(p_map, p_pca, nrow=1, rel_widths = c(0.4,0.6))
save_plot('map_pca.pdf', p_map_pca, nrow=1, base_width=12, base_heigh=4)


# Overall concordance -----------------------------------------------------

setwd('/Users/alicia/daly_lab/neurogap/data/high_coverage/gencove/stats/side_by_side')
setwd('/Users/alicia/daly_lab/neurogap/data/high_coverage/stats')

depths <- c('0.5', '1.0', '2.0', '4.0', '6.0')
chips <- c('H3Africa', 'Omni2.5', 'MEGA', 'PsychChip', 'GSA')

depth_concordance <- function(depth, prefix='neurogap_intersect_', suffix='X.imputed_') {
  concordance <- read.table(paste0(prefix, depth, suffix, 'concordance_variants.tsv'), header=T) %>%
    filter(snp=='true')
  snp_count <- sum(concordance$total_concordant, concordance$total_discordant)
  snp_concordance <- sum(concordance$total_concordant) / (sum(concordance$total_concordant) + sum(concordance$total_discordant))
  return(list(depth=depth, snp_count=snp_count, snp_concordance=snp_concordance))
}

imputed_concordance_overall <- map_dfr(depths, depth_concordance)
#imputed_concordance_overall <- map_dfr(depths, function(x) {depth_concordance(x, prefix='neurogap_')})
gencove_concordance_overall <- map_dfr(depths, function(x) {depth_concordance(x, prefix='merge_intersect_', suffix='_')})
chip_concordance_overall <- map_dfr(chips, function(x) {depth_concordance(x, prefix='NeuroGap_30x_', suffix='.imputed_')})


raw_concordance <- map_dfr(covs, read_file)
refined_concordance <- map_dfr(lowcovs, function(x) {read_file(x, refined=TRUE)} )
