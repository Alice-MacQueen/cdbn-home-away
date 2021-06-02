library(here)
library(magrittr)
library(sf)
library(ggplot2)
library(gridExtra)
library(wesanderson)

bs_plot = here('results', 'figure data', 'Fig S3 - faceted best breeding sites.gpkg') %>% 
  st_read
basemap =   here('results', 'figure data', 'CDBN Basemap.gpkg') %>% 
  st_read
colors = wes_palette('FantasticFox1')

popn = 'Race'
site = 'Location_code'


bs_plot$label = apply(bs_plot, 1, function(x) ifelse(x['plotval']==1, x[site], ''))

best_sites_facet_plot = ggplot(basemap) +
  facet_wrap(paste('~', popn) %>% formula,
             nrow=1) +
  scale_color_manual(values=colors, guide=NULL) +
  geom_sf(fill='white',
          size=0.2) +
  geom_sf_text(data=bs_plot,
               aes_string(label='label',
                          color=popn),
               size=3) +
  theme_minimal() +
  labs(fill=expression(paste('Mean Yield [kg ', 'ha'^-1, ']')),
       parse=TRUE,
       x=NULL,
       y=NULL) +
  theme(axis.ticks=element_line(color='gray',
                                size=0.2),
        axis.title=element_text(size=9),
        legend.title=element_text(size=9),
        panel.grid=element_line(color='gray',
                                size=0.2),
        legend.position='bottom')

out = "Figure S3 - Best Breeding Sites by Race.jpg" %>% 
  here('results', .)
jpeg(out, 
     width=6.5, 
     height=2.5,
     units='in', res=300)
grid.arrange(best_sites_facet_plot)
dev.off()