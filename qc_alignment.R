library(tidyverse)
library(cowplot)

setwd('/Users/alicia/daly_lab/neurogap/data/high_coverage')

pdo <- read.csv('PDO-18489_Seq Metrics_FINAL_071519.csv', header=T) %>%
  filter(!External.ID %in% c('NGE0018', 'NGE0130'))

mean(pdo[,"PF.Reads.Aligned.."])
sd(pdo[,"PF.Reads.Aligned.."])

pdo_1X <- read.csv('../low_coverage/PDO-17659_Seq_Metrics.csv', header=T)
mean(pdo_1X[,"PF.Reads.Aligned.."])
sd(pdo_1X[,"PF.Reads.Aligned.."])
mean(pdo_1X[,"Mean.Coverage..Raw..Full"])
sd(pdo_1X[,"Mean.Coverage..Raw..Full"])
p1 <- ggplot(pdo_1X, aes(x=Mean.Coverage..Raw..Full)) + 
  geom_density(fill='steelblue') +
  labs(x='Coverage', y='Density') +
  geom_vline(xintercept=1, linetype='dashed') +
  theme_classic() +
  theme(text = element_text(size=16, color='black'),
        axis.text = element_text(color='black'))

p2 <- ggplot(pdo_1X, aes(x=Mean.Coverage..Raw..Full, y=PF.Reads.Aligned..)) +
  geom_point(color='steelblue') +
  labs(x='Coverage', y='% reads aligned') +
  geom_vline(xintercept=1, linetype='dashed') +
  theme_classic() +
  theme(text = element_text(size=16, color='black'),
        axis.text = element_text(color='black'))

cor.test(pdo_1X$Mean.Coverage..Raw..Full, pdo_1X$PF.Reads.Aligned..)
p <- plot_grid(p1, p2, labels = c('A', 'B'), align='h')
ggsave('../low_coverage/cov_alignment.pdf', p, width=12, height=6)
ggsave('../low_coverage/cov_alignment.png', p, width=12, height=6)


ggplot(pdo_1X, aes(x=Coverage_Adjustment)) + geom_density()
