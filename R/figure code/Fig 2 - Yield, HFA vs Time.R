# Figure 2: Trends in yield and HFA across time.
library(here)
library(magrittr)
library(ggplot2)
library(ggpmisc)
library(gridExtra)
library(wesanderson)

popn = 'Race'
pheno = 'SY'
year = 'Year'
colors = wes_palette('FantasticFox1')

# Fig 2a Data
df_2a = here('results', 'figure data', 'Fig 2a - yield by year.csv') %>% 
  read.csv(stringsAsFactors=FALSE)

# Fig 2b Data
thfa_df = here('results', 'figure data', 'Fig 2b - HFA across time.csv') %>% 
  read.csv(stringsAsFactors=FALSE)
annotate = here('results', 'figure data', 'Fig 2b - HFA across time Annotations.csv') %>% 
  read.csv(stringsAsFactors=FALSE)
exp_hfa = here('results', 'figure data', 'Fig 2b - HFA across time Expectation.csv') %>% 
  read.csv(stringsAsFactors=FALSE)

# Fig 2a - Yield ~ Time
se = function(x) sd(x) / sqrt(length(x))
se_high = function(x) mean(x) + se(x)
se_low = function(x) mean(x) - se(x)

yield_yr = ggplot(df_2a,
                  aes_string(x='Year0', 
                             y=pheno,
                             color=popn,
                             fill=popn,
                             grp.label=popn)) +
  facet_wrap(paste('~', popn) %>% formula,
             nrow=1) +
  scale_x_continuous(breaks=seq(0, 35, by=10),
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
       y = 'Yield [kg/ha]',
       color=NULL,
       fill=NULL) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.ticks=element_line(color='gray',
                                size=0.2),
        axis.title=element_text(size=9),
        axis.text.x=element_text(angle=90,
                                 hjust=0.5,
                                 vjust=0.5))

# Fig 2b - HFA ~ Time
year_num = paste(year, 'num', sep='_')
xstep = 10

xrange = thfa_df[, year_num] %>% 
  range
xmin = xrange[1] %>% 
  divide_by(xstep) %>% 
  floor %>% 
  multiply_by(xstep)
xmax = xrange[2]
xbreaks = seq(xmin, xmax, by=xstep)
xbreaks_minor = xbreaks + (xstep/2)


thfa_plot = thfa_df %>% 
  ggplot(aes_string(x=paste0(year, '_num'),
                    y='Estimate',
                    color=popn,
                    fill=popn,
                    grp.label=popn)) +
  facet_wrap(paste('~', popn) %>% formula,
             nrow=1) +
  scale_color_manual(values=colors, guide=NULL) +
  scale_fill_manual(values=colors, guide=NULL) +
  scale_x_continuous(breaks=xbreaks, minor_breaks=xbreaks_minor) +
  geom_hline(aes(yintercept=0),
             size=0.2,
             color='gray') +
  geom_ribbon(data=exp_hfa,
              aes_string(ymin='p05',
                         ymax='p95',
                         x='x_val'),
              fill='#EEEEEE',
              alpha=1,
              inherit.aes=FALSE) +
  geom_pointrange(aes(ymin=Estimate - Std.Error,
                      ymax=Estimate + Std.Error),
                  shape=20,
                  size=0.2) +
  geom_smooth(method='lm',
              formula=y~x,
              alpha=0.1,
              size=0.2,
              level=0.95,
              fullrange=TRUE) +
  geom_text(data=annotate,
            aes(x=x,
                y=y,
                label=text),
            color='black',
            hjust='inward',
            vjust='inward',
            parse=TRUE,
            size=3) +
  labs(x='Year',
       y='Home Field Advantage [kg/ha]',
       color=NULL) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.title=element_text(size=9),
        axis.ticks=element_line(color='gray',
                                size=0.2))

# Combine
fig2a = yield_yr +
  ggtitle('A)') +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x=element_line(size=0.2, color='gray'),
        plot.title=element_text(size=10),
        strip.text=element_text(face='bold'))
fig2b = thfa_plot +
  ggtitle('B)') +
  theme(strip.text=element_blank(),
        plot.title=element_text(size=10))

grid.arrange(fig2a, fig2b, nrow=2)

# Export
out = "Figure 2 - Yield, HFA Trends.jpg" %>% 
  here('results', .)
jpeg(out, width=6.5, height=5, units='in', res=300)
grid.arrange(fig2a, fig2b, nrow=2)
dev.off()
