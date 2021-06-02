# Figure 3b - heritability vs year

library(here)
library(magrittr)
library(ggplot2)
library(ggpmisc)
library(gridExtra)
library(wesanderson)

popn = 'Race'
colors = wes_palette('FantasticFox1')

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
        axis.title=element_text(size=9))

H_yr = H_yr + 
  ggtitle('B)') + 
  theme(plot.title=element_text(size=10))


"Figure 3b - Heritability vs Time.jpg" %>% 
  here('results', .) %>% 
  jpeg(width=6.5, height=2.5, units='in', res=300)
grid.arrange(H_yr)
dev.off()
