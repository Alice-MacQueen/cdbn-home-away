# Various summary statistics about heritability and/or yield and time.

library(ggplot2)
library(ggpmisc)
library(magrittr)
library(here)
library(wesanderson)
library(gridExtra)

pheno = 'SY'
popn = 'Race'
colors = wes_palette('FantasticFox1')

H_stats = here('results', 'figure data', 'Fig S1 - Heritability Site Statistics.csv') %>% 
  read.csv(stringsAsFactors=FALSE)
yield_H = here('results', 'figure data', 'Fig S1 - Heritability vs Yield.csv') %>% 
  read.csv(stringsAsFactors=FALSE)

gg_themer = function(x, remove_strip=FALSE) {
  require(ggplot2)
  out = x + 
    theme_minimal() +
    theme(panel.grid=element_blank(),
          axis.ticks=element_line(color='gray',
                                  size=0.2),
          axis.title=element_text(size=9))
  
  if (remove_strip) {  # for the bottom panel of fully faceted figures
    out = out + theme(strip.text=element_blank())
  }
  return(out)
}

sd_H_plt = ggplot(H_stats,
                  aes_string(x='site_year',
                             y='sd_H')) +
  geom_point(size=0.5) +
  geom_smooth(method='lm',
              formula='y~x',
              alpha=0.2,
              size=0.2,
              color='darkgray',
              fill='lightgray') +
  stat_poly_eq(
    aes(label=paste('atop(',
                    stat(eq.label),
                    ',italic(p) == ', 0.014, ')',
                    sep="")),
    formula=y~x,
    parse=TRUE,
    label.y.npc=0.0,
    label.x.npc=0.9, 
    size=3,
    color='black'
  ) +
  labs(x='Year  at Location',
       y='Heritability Standard Deviation')

sd_H_plt %<>% 
  gg_themer

sd_H_plt



label_df = data.frame(x=0.35,
                      y=0.05,
                      label='n.s.')
mean_H_plt = ggplot(H_stats,
                    aes_string(x='sd_H',
                               y='mean_H')) +
  geom_point(size=0.5) +
  geom_text(data=label_df,
            aes(x=x,
                y=y,
                label=label),
            size=3,
            color='black',
            inherit.aes=FALSE) +
  labs(x='Heritability Standard Deviation',
       y='Mean Heritability')

mean_H_plt %<>% 
  gg_themer


yield_H_plt = ggplot(data=yield_H,
                     aes_string(x='Heritability',
                                y=pheno)) +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors, guide=NULL) +
  # facet_wrap(paste('~', popn) %>% formula) +
  geom_point(aes_string(color=popn),
             size=0.5) +
  geom_smooth(aes_string(x='Heritability',
                         y=pheno),
              method='lm',
              formula='y~x',
              color='darkgray',
              fill='lightgray',
              alpha=0.2,
              size=0.2) +
  stat_poly_eq(
    aes(label=paste('atop(',
                    stat(eq.label),
                    ',italic(p) == ', 0.012, ')',
                    sep="")),
    formula=y~x,
    parse=TRUE,
    label.y.npc=0,
    label.x.npc=0.95, 
    size=3,
    color='black'
  ) +
  labs(x='Heritability',
       y='Mean Yield [kg/ha]',
       color=NULL) 

yield_H_plt %<>% gg_themer


"Figure S1 - Heritability Statistics and Relationships.jpg" %>% 
  here('results', .) %>% 
  jpeg(width=4.5, height=5, units='in', res=300)
grid.arrange(grobs=list(sd_H_plt, mean_H_plt, yield_H_plt),
             widths=c(1,2,2,1),
             layout_matrix=rbind(c(1,1,2,2),
                                 c(3,3,3,NA))
)
dev.off()