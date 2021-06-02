# Figure 1: Yield by location and race
library(here)
library(magrittr)
library(ggplot2)
library(wesanderson)
library(sf)

popn = 'Race'
pheno = 'SY'
colors = wes_palette('FantasticFox1')


# load points data
site_map = here('results', 'figure data', 'Figure 1 - Yield by Site.gpkg') %>% 
  st_read

# rnaturalearth ne_countries basemap, cropped and projected to match site_map
basemap =   here('results', 'figure data', 'CDBN Basemap.gpkg') %>% 
  st_read

# map
site_plot = ggplot(basemap) +
  facet_wrap(paste('~', popn) %>% formula,
             nrow=1) +
  scale_fill_gradientn(colors=colors[c(3,2,5)], 
                       limits=range(data.frame(site_map)[, pheno])) + 
  geom_sf(fill='white',
          size=0.2) +
  geom_sf(data=site_map,
          aes_string(fill=pheno),
          shape=21,
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
        # axis.text=element_blank(),
        panel.grid=element_line(color='gray',
                                size=0.2),
        legend.position='bottom')
site_plot  

out = "Figure 1 - Yield by Site.jpg" %>% 
  here('results', .)
jpeg(out, width=6.5, height=3, units='in', res=300)
dev.off()
