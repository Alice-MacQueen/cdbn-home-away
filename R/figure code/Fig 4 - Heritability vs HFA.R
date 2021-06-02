# HFA vs heritability, annotated into quadrants.

library(here)
library(ggplot2)
library(gridExtra)
library(magrittr)
library(wesanderson)
library(sf)

popn = 'Race'
colors = wes_palette('FantasticFox1')

site_H_HFA = here('results', 'figure data', 'Fig 4 - HFA vs H2.csv') %>% 
  read.csv(stringsAsFactors=FALSE)
bssp = here('results', 'figure data', 'Fig 4 - best breeding sites.gpkg') %>% 
  st_read
basemap =   here('results', 'figure data', 'CDBN Basemap.gpkg') %>% 
  st_read

annotater = data.frame(label=c('II', 'I', 'III', 'IV') %>% 
                         paste('bold(', ., ')', sep=''),
                       x=c(0.05, 0.95, 0.05, 0.95),
                       y=c(1400, 1400, 0, 0))

medians=c(H = median(site_H_HFA$H, na.rm=TRUE),
          HFA = median(site_H_HFA$HFA, na.rm=TRUE))

site_hfa_plot = 
  ggplot(site_H_HFA, 
         aes_string(x='H',
                    y='HFA',
                    color=popn)) +
  scale_color_manual(values=colors, guide=NULL) +
  facet_wrap(paste('~', popn) %>% formula,
             nrow=1) +
  geom_vline(aes(xintercept=medians['H']),
             color='gray',
             size=0.2) +
  geom_hline(aes(yintercept=medians['HFA']),
             color='gray',
             size=0.2) +
  geom_point(size=0.5) +
  geom_text(data=annotater,
            aes(label=label,
                x=x,
                y=y),
            hjust='inward',
            vjust='inward',
            size=3,
            inherit.aes=FALSE,
            parse=TRUE) +
  theme_minimal() +
  labs(x='Heritability',
       y='Home Field Advantage [kg/ha]') +
  theme(panel.grid=element_blank(),
        axis.ticks=element_line(color='gray',
                                size=0.2),
        axis.title=element_text(size=9))
site_hfa_plot

best_sites_plot = ggplot(basemap) +
  geom_sf(fill='white',
          size=0.2) +
  geom_sf_text(data=bssp,
               aes(label=label),
               size=3,
               alpha=0.5,
               color='black') +
  theme_minimal() +
  labs(fill=expression(paste('Mean Yield [kg ', 'ha'^-1, ']')),
       parse=TRUE) +
  theme(axis.ticks=element_line(color='gray',
                                size=0.2),
        axis.title=element_text(size=9),
        legend.title=element_text(size=9),
        panel.grid=element_line(color='gray',
                                size=0.2),
        legend.position='bottom')

out = "Figure 4 - HFA vs H2.jpg" %>% 
  here('results', .)
jpeg(out, 
     width=6.5, 
     height=6,
     units='in', res=300)
grid.arrange(grobs=list(site_hfa_plot, 
                        best_sites_plot),
             heights=c(1,2))
dev.off()