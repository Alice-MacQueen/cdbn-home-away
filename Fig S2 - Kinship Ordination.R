# Figure S2 - Ordinations

library(here)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(wesanderson)

popn = 'Race'
colors = wes_palette('FantasticFox1')

PC_scree = here('results', 'figure data', 'Fig S2 - ordination scree.csv') %>% 
  read.csv(stringsAsFactors=FALSE)
ordplots= here('results', 'figure data', 'Fig S2 - ordination ax1-2 scores.csv') %>% 
  read.csv(stringsAsFactors=FALSE)


screeplot = ggplot(PC_scree,
                   aes_string(x='PC',
                              y='pVar',
                              fill=popn)) +
  scale_x_continuous(breaks=seq(1, 10, by=2)) +
  facet_wrap(paste('~', popn) %>% formula) +
  scale_fill_manual(values=colors, guide=NULL) +
  geom_col(position='dodge') +
  labs(x = 'Principal Coordinate Axis',
       y = 'Variance [%]',
       color=NULL,
       fill=NULL) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        panel.grid.major.y=element_line(size=0.2, color='white'),
        panel.ontop=TRUE,
        axis.ticks=element_line(color='gray',
                                size=0.2),
        strip.background=element_blank(),
        strip.text=element_blank(),
        axis.title=element_text(size=9))


ordination = ggplot(ordplots,
                    aes_string(x='PC1',
                               y='PC2',
                               color=popn)) +
  scale_color_manual(values=colors, guide=NULL) +
  coord_equal() +
  facet_wrap(paste('~', popn) %>% formula) +
  geom_hline(aes(yintercept=0),
             size=0.2,
             color='gray') +
  geom_vline(aes(xintercept=0),
             size=0.2,
             color='gray') +
  geom_point(size=0.5) +
  labs(x='PC 1',
       y='PC 2') +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.ticks=element_line(color='gray',
                                size=0.2),
        strip.text=element_text(size=10,
                                face='bold'),
        axis.title=element_text(size=9))

# Combine and export
p1 = ordination + ggtitle('A)')
p2 = screeplot + ggtitle('B)')

out = "Figure S2 - Kinship Ordination and Scree.jpg" %>% 
  here('results', .)
jpeg(out, width=6.5, height=5, units='in', res=300)
grid.arrange(p1, p2, nrow=2)
dev.off()