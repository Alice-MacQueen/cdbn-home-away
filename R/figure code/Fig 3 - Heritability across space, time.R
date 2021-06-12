# Figure 3b - heritability vs year

library(here)
library(magrittr)
library(ggplot2)
library(ggpmisc)
library(gridExtra)
library(wesanderson)
library(stars)
library(sf)

popn = 'Race'
colors = wes_palette('FantasticFox1')

in_dir = file.path('results', 'figure data')

sites  = here(in_dir, 'Figure 1 - Yield by Site.gpkg') %>% 
  st_read
basemap =   here(in_dir, 'CDBN Basemap.gpkg') %>% 
  st_read
H_map = here(in_dir, 'kriged_bean_heritability.tif') %>%
  read_stars(proxy=FALSE)

## Figure 3a: Heritability across space
# buffer to crop
site_buff = sites %>% 
  st_union %>% 
  st_buffer(100000) %>% 
  st_convex_hull %>% 
  st_transform(st_crs(H_map))

H_map %<>%
  st_warp(cellsize=0.1,
          method='bilinear',
          use_gdal=TRUE) %>%
  st_crop(site_buff) %>% 
  st_transform(st_crs(basemap))
# basemap %<>% st_transform(st_crs(H_map))

# st_as_sf %>% 
# set_colnames(c('value', 'geometry'))

H_space = ggplot(basemap) +
  scale_fill_gradientn(colors=colors[c(3,2,5)]) +
  geom_sf(fill='white',
          size=0.2) +
  geom_stars(data=H_map,
             alpha=0.8) +
  geom_sf(data=sites,
          color='black',
          size=0.1) +
  labs(fill='Yield\nHeritability') +
  theme_minimal() +
  theme(axis.ticks=element_line(color='gray',
                                size=0.2),
        axis.title=element_text(size=9),
        legend.title=element_text(size=9),
        panel.grid=element_line(color='gray',
                                size=0.2)) +
  ggtitle('A)') +
  theme(plot.title=element_text(size=10,
                                hjust=-0.38))
H_space

## Figure 3b: Heritability across time

se = function(x) sd(x) / sqrt(length(x))
se_high = function(x) mean(x) + se(x)
se_low = function(x) mean(x) - se(x)

H_df =   here('results', 'figure data', 'Fig 3b - heritability vs year.csv') %>% 
  read.csv(stringsAsFactors=FALSE)

H_yr = ggplot(H_df,
              aes_string(x='Y0',
                         y='Heritability',
                         color=popn,
                         fill=popn)) +
  facet_wrap(paste('~', popn) %>% formula,
             nrow=1) +
  scale_x_continuous(limits=c(0,35),
                     breaks=seq(0, 35, by=10),
                     labels=seq(1980, 2015, by=10)) +
  scale_color_manual(values=colors, guide=NULL) +
  scale_fill_manual(values=colors, guide=NULL) +
  stat_summary(fun=mean,
               fun.max=se_high,
               fun.min=se_low,
               shape=20,
               size=0.2) +
  geom_smooth(method='lm',
              formula=y~x,
              alpha=0.1,
              size=0.2,
              fullrange=TRUE,
              level=0.95) +
  stat_poly_eq(
    aes(label=paste('list(',
                    stat(eq.label),
                    ',italic(p) < 0.001)',
                    sep="")),
    formula=y~x,
    parse=TRUE,
    label.y.npc=0.9,
    size=3,
    color='black'
  ) +
  labs(x = 'Year',
       y = 'Yield Heritability [au]',
       color=NULL,
       fill=NULL) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.ticks=element_line(color='gray',
                                size=0.2),
        axis.title=element_text(size=9)) +
  ggtitle('B)') + 
  theme(plot.title=element_text(size=10))
H_yr


"Figure 3 - Heritability across Space, Time.jpg" %>% 
  here('results', .) %>% 
  jpeg(width=6.5, height=5, units='in', res=300)
grid.arrange(grobs=list(H_space,
                        H_yr),
             heights=c(5, 4))
dev.off()
