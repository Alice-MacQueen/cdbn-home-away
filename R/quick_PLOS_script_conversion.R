#! /bin/R
#
# PME Updated June 17, 2019

libraries = c('lme4', 'lattice', 'car', 'reshape2', 'magrittr', 'quantreg', 'ggplot2')
for (i in libraries) library(i, character.only=TRUE)

#### Housekeeping ####

dir_in = "/home/patrick/Cloud/Side Projects/CDBN HFA/Data"
tall_in = "Phenotypes_of_sequenced_individuals_home_away_analysis.csv"
factors = c('Location_code', 'Year', 'Taxa', 'CDBN_ID', 'Seq_ID', 'Race', 'Origin')

race_to_origin = c(Durango = 'Mesoamerican',
                   Jalisco = 'Mesoamerican',
                   Mesoamerican = 'Mesoamerican',
                   `Nueva Granada` = 'Andean')



df = file.path(dir_in, tall_in) %>% 
  read.csv(stringsAsFactors=FALSE)

df$Origin = df[,'Race'] %>% 
  as.character %>%
  race_to_origin[.]



df[, factors] %<>% lapply(as.factor)
df[, factors] %<>% lapply(tolower)  # semi-important!! A few names aren't consistent typecases.

names(df) %<>%
  gsub('Location_code', 'LOCATION', .) %>% 
  gsub('Year', 'YEAR', .) %>% 
  gsub('Taxa', 'GENOTYPE', .) %>% 
  gsub('SY', 'YIELD', .)
df[, 'LOCYR'] = paste(df[, 'LOCATION'], df[, 'YEAR'], sep="_") %>% 
  as.factor

# Remove rare sites (cannot be estimated for site or site-year effects)
tt = df[, c('LOCATION', 'YEAR')] %>% 
  unique %>% 
  aggregate(YEAR ~ LOCATION, ., length) %>% 
  subset(YEAR < 3, select=LOCATION)
df %<>% subset(!(LOCATION %in% tt[, 1]))

# Convert bu/ac to kg/ha (15.5%, 56 lb bu)
# df[, 'YIELD'] %<>% multiply_by(.06262) # Mg/ha

# factors
# facts = c('REGION', 'YEAR', 'COMPANY', 'HYBRID', 'LOCATION', 'GENOTYPE', 'LOCYR')
# for(i in facts) df[, i] %<>% as.factor

#### First Look ####
# densityplot(~YIELD | YEAR, group=REGION, data=df, auto.key=TRUE)


#### Model ####
# The linear model is:
#   Yield ~ Hybrid (genotype) + Location + Is_Home + error - for each year.
#   Is_Home ~ Time
#
# To do this:
#   1. Subset checks
#   2. ID the home location. This is the location where the hybrid does the best relative to others.
#     a. Center/scale by year and location
#     b. Aggregate mean yield by hybrid and location
#     c. ID max location for each hybrid.
#     d. Flag home location in checks df
#   3. Run model for each year
#   4. Extract home_advantage coefficient
#   5. Regress Is_Home ~ Time

#### ID CHECKS ####

# Varieties planted at each site-year
# id_checks = split(df, df$YEAR) %>% 
#   lapply(function(x) {
#     nsites = length(unique(x$LOCATION)) # sites planted in year X
#     at_all_sites = aggregate(YIELD ~ YEAR + GENOTYPE,
#                              data=x,
#                              function(x) length(x) == nsites) # is planted at all sites?
#     return(at_all_sites)
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame
# names(id_checks)[3] = 'IS_CHECK'

# checks = merge(df, id_checks, by=c('YEAR', 'GENOTYPE')) %>% 
#   subset(IS_CHECK)

#### Find home location ####
# Find home location - the highest average relative yield (standardized mean=0, sd=1)
# Relative yields within location-year

# Center yields by location and year
checks = df
centered = split(checks, checks[, 'LOCYR']) %>% 
  lapply(function(x) {
    x[, 'YIELD'] %<>% scale
    return(x) 
  }) %>% 
  do.call(rbind, .)

# Find highest relative yield for each genotype
centered = split(centered, centered[, 'GENOTYPE']) %>% 
  lapply(function(x) {
    # aggregate to mean (relative) yield within location across years. Preserve genotype column.
    aggd = aggregate(YIELD ~ LOCATION + GENOTYPE, data=x, function(y) mean(y, na.rm=T))
    names(aggd)[3] = 'YIELD'
    
    # find hightest mean yield
    tt = aggd[, 'YIELD'] == max(aggd[, 'YIELD'])
    aggd[, 'IS_HOME'] = tt
    aggd[, 'YIELD'] = NULL
    
    # Merge back to non-aggregated df
    x = merge(x, aggd, by=c('LOCATION', 'GENOTYPE'))
    return(x)
  }) %>% 
  do.call(rbind, .)
rownames(centered) = NULL

# Send home info back to original df
merge_by = c('GENOTYPE', 'YEAR', 'LOCATION')
checks1 = centered[, c(merge_by, 'IS_HOME')] %>% 
  merge(checks, ., by=merge_by)

checks = checks1; rm(checks1)

#### G*E, across years ####
# checks[, 'GENO_LOC'] = paste(checks[, 'GENOTYPE'], checks[, 'LOCATION'], sep="_") %>%
#   as.factor
# check_sub = duplicated(checks$GENO_LOC) %>%
#   which %>%
#   c(duplicated(checks$GENO_LOC, fromLast=TRUE) %>%
#       which) %>%
#   unique %>%
#   checks[., ]
# 
# # Years per genotype
# geno_years = split(check_sub, check_sub$GENOTYPE) %>%
#   sapply(function(x) length(unique(x$YEAR))) %>%
#   as.matrix %T>% 
#   print
# Only three were planted more than two years. That's not enough to do a G*E

facts = c('GENOTYPE', 'LOCATION', 'YEAR')
for(i in facts) checks[, i] %<>% as.factor

checks = subset(checks, Origin=='andean')


Var_Mod = lm(YIELD ~ GENOTYPE + LOCATION + YEAR + YEAR:LOCATION, data=checks)
a = Anova(Var_Mod, type="II") %T>% # G, E, Y, E*Y highly significant. G*Y marginal. G*E not significan.t
  print
qqPlot(residuals(Var_Mod))
influenceIndexPlot(Var_Mod)
ss_Var = data.frame(
  PREDICTOR = rownames(a),
  SUMSQ = a$`Sum Sq`,
  pVAR = round(a$`Sum Sq` / sum(a$`Sum Sq`) *100, 2)
) %T>% print

#### Local Adaptation Model ####
LA_mod = lm(YIELD ~ YEAR + GENOTYPE + LOCATION + 
              LOCATION:YEAR + 
              IS_HOME, 
            data=checks)

qqPlot(residuals(LA_mod))
influenceIndexPlot(LA_mod)
ss = summary(LA_mod)
ss$adj.r.squared
a = Anova(LA_mod, type="II") %T>% print

ss_LA = data.frame(
  PREDICTOR = rownames(a),
  SUMSQ = a$`Sum Sq`,
  pVAR = round(a$`Sum Sq` / sum(a$`Sum Sq`) *100, 2)
) %T>% print

HFA = summary(LA_mod)$coefficients['IS_HOMETRUE',]
mean(checks$YIELD)
sd(checks$YIELD)/sqrt(length(checks$YIELD))
HFA['Estimate']/mean(checks$YIELD)*100
HFA


# Compare
anova(Var_Mod, LA_mod)
AIC(LA_mod, Var_Mod)
BIC(LA_mod, Var_Mod)

# Local Adaptation Effect Size
print("Local Adaptation Effect Size")
ss$coefficients["IS_HOMETRUE", ]

# Local Adaptation Effect Relative
print("Local Adaptation Effect Relative")
ss$coefficients["IS_HOMETRUE", 'Estimate']/ss$coefficients['(Intercept)', 'Estimate']

anv = Anova(LA_mod, type="II")
anv['IS_HOME', 'Sum Sq'] / sum(anv[, 'Sum Sq'])


#### Donut Plot of Variance Partitioning! ####
sum(ss_LA$pVAR)
ss_LA$Source = c('Year', 'Genotype', 'Location', 'Home Site', 'Location\nby Year', 'Residuals')
source_order = c(3,5,2,4,1,6)
ss_LA$Source = factor(ss_LA$Source, levels = ss_LA$Source[source_order])
pltdf = ss_LA[source_order,]
mn = pltdf$pVAR
for (i in 2:length(mn)) {
  mn[i] = min(mn[i] + mn[i-1], 100)
}
pltdf$MX = mn
pltdf$MN = pltdf$MX - pltdf$pVAR

ww = FALSE
# ww = TRUE
if(ww) {
  fout = paste("PLOT DATA - Yield Variance Partitioning Donut ",
               Sys.Date(),
               ".csv",
               sep="") %>% 
    file.path(dir_in, "Figures", .)
  write.csv(pltdf, fout, row.names=FALSE)
}

donut = ggplot(pltdf, aes(fill=Source,
                          xmin=2,
                          xmax=5,
                          ymin=MN,
                          ymax=MX)) +
  scale_fill_manual(values=c('steelblue1', 'steelblue2', 'steelblue3', 'tomato2', 'steelblue4', 'gray'),
                    name='Source of Yield\nVariation') +
  geom_rect() +
  coord_polar(theta='y') +
  xlim(c(0, 5)) +
  theme_light() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid=element_blank(),
        panel.border=element_blank())
donut

png(file.path(dir_in, 'Donut Variance.pdf'), height=4, width=5, units='in', res=300)
donut
dev.off()

#### Local Adaptation - By Year ####
yy = checks[, c('YEAR', 'YIELD', 'GENOTYPE', 'LOCATION', 'IS_HOME')]
yy[, c('YEAR', 'GENOTYPE', 'LOCATION')] %<>% lapply(droplevels)
yy %<>%
  split(yy[, 'YEAR'])

# Year local adaptation
local_adapt = lapply(yy, function(x) {
  if (sum(x$IS_HOME) >= 3) {
    # print(x[1, 'YEAR'])
    mod = lm((YIELD) ~ GENOTYPE + LOCATION + IS_HOME, data=x) # robust regression (median quantile) each year
    vals = summary(mod)#, se='boot')                            # bootstrap standard errors
    vals = vals$coefficients[c('(Intercept)', 'IS_HOMETRUE'), c('Estimate', 'Std. Error')]  # make a table
    # mod = lm((YIELD) ~ GENOTYPE + LOCATION + IS_HOME, data=x) # robust regression (median quantile) each year
    # vals = summary(mod)#, se='boot')                            # bootstrap standard errors
    # vals = vals$coefficients[c('(Intercept)', 'IS_HOMETRUE'), c('Estimate', 'Std. Error')]  # make a table
    # mod = lmer((YIELD) ~ (1|GENOTYPE) + (1|LOCATION) + IS_HOME, data=x) # robust regression (median quantile) each year
    # vals = summary(mod)#, se='boot')                            # bootstrap standard errors
    # vals = vals$coefficients[c('(Intercept)', 'IS_HOMETRUE'), c('Estimate', 'Std. Error')]  # make a table
    return(vals)
  } else {
    return(NULL)
  }
}) 

is_null = sapply(local_adapt, is.null)
local_adapt %<>% .[!is_null]

# Make a nice dataframe
local_adapt = lapply(names(local_adapt), function(x) {
  out = local_adapt[[x]] %>% 
    data.frame
  out$MEASURE = rownames(out)
  out$YEAR = x
  return(out)
}) %>% 
  do.call(rbind, .) %>% 
  data.frame %>% 
  split(.$MEASURE) %>% 
  lapply(function(x) {
    names(x) = c('VALUE', 'ERROR', 'MEASURE', 'YEAR')
    x$YEAR %<>% 
      as.character %>% 
      as.numeric
    x$YEAR_NUMBER = x$YEAR - min(x$YEAR)
    return(x)
  })

home = local_adapt[['IS_HOMETRUE']]
intercept = local_adapt[['(Intercept)']]

mod = lm(VALUE ~ YEAR_NUMBER, data=home)
summary(mod)

print("#####    Mean Local Adaptation Across Years    #####")
mean(home$VALUE)

mn = mean(home$VALUE)
se = sd(home$VALUE)/sqrt(length(home$VALUE)-1)
# mn_vec = list(mn, se, 'IS_HOMETRUE', 'Overall\nMean', 20, TRUE)
# home$IS_MEAN = FALSE
# home = as.data.frame(
#   rbind(home, mn_vec))
# home$YEAR = as.factor(home$YEAR)


overwrite=FALSE
if(overwrite) {
  write.csv(home, 
            file.path(dir_in,
                      paste("PLOT DATA - Local Adaptation by Year ", 
                            Sys.Date(), 
                            ".csv",
                            sep="")),
            row.names=FALSE)
}

#

#### Plot Local Adaptation By Year ####

home$YEAR = as.factor(home$YEAR)
home$YEAR_NUM = as.numeric(home$YEAR)
plt = ggplot(home) +
  # scale_fill_manual(values=c('steelblue3', 'tomato2')) +
  # ylim(c(0, 1.5)) +
  xlab("Year") +
  ylab(expression(paste("Observed Home Field Advantage [Mg ", "ha"^"-1", "]"))) +
  # geom_ribbon(aes(x=YEAR,
  #                 ymin=mn-se,
  #                 ymax=mn+se)) +
  geom_hline(aes(yintercept=mn+se),
             color='#999999',
             linetype='dashed') +
  geom_hline(aes(yintercept=mn-se),
             color="#999999",
             linetype='dashed') +
  geom_errorbar(aes(x=YEAR,
                    ymin=VALUE-ERROR,
                    ymax=VALUE+ERROR),
                width=0.0,
                color="#666666") +
  geom_point(aes(x=YEAR,
                 y=VALUE),
             color='steelblue3',
             size=4) +
  theme_light() +
  theme(axis.text.x=element_text(angle=45,
                                 vjust=0.6))
plt

png(file.path(dir_in, 'home_field_advantage_yrs.png'), height=4, width=5, units='in', res=300)
plt
dev.off()

#### Yield across years ####
yields = aggregate(YIELD ~ YEAR, data=checks, mean) %>% 
  data.frame(ERROR = aggregate(YIELD ~ YEAR, 
                               data=checks,
                               function(x) sd(x)/sqrt(length(x-1)))[, 2])
yields$YEAR %<>% as.character %>% as.numeric


ff = formula('YIELD ~ YEAR')
yield_year_mod = lm(ff, data=yields)
summary(yield_year_mod)
Anova(yield_year_mod)

#### Plot Yield by Year ####
w = FALSE
# w = TRUE
if(w) {
  fout = paste("PLOT DATA - Yield vs Year ",
               Sys.Date(),
               ".csv", sep="") %>% 
    file.path(dir_in, .)
  write.csv(yields, 
            fout,
            row.names=FALSE)
}

plt = ggplot(yields) +
  xlab("Year") +
  ylab("Yield [Mg/ha]") +
  geom_errorbar(aes(x=YEAR,
                    ymin=YIELD-ERROR,
                    ymax=YIELD+ERROR),
                width=0.2,
                color="#666666") +
  geom_smooth(aes(x=YEAR,
                  y=YIELD),
              color='tomato4',
              method='lm',
              weight=0.1,
              linetype="dashed") +
  geom_point(aes(x=YEAR,
                 y=YIELD),
             color='steelblue3',
             size=4) +
  theme_light()
plt

png(file.path(dir_in, 'yield across years.png'), height=4, width=5, units='in', res=300)
plt
dev.off()
#### Local adaptation vs yield anomoly ####
# Yield anomaly = residuals of yield ~ year
# Local adaptation as caluclated in each year
LA = local_adapt[['IS_HOMETRUE']][, c('YEAR', 'VALUE', 'ERROR', 'YEAR_NUMBER')]
names(LA)[2] = 'LOCAL_ADAPTATION'
YA = residuals(yield_year_mod) %>% 
  as.data.frame
names(YA) = 'YIELD_ANOMOLY'
YA$YEAR_NUMBER = rownames(YA) %>% 
  as.numeric %>% 
  subtract(1)
# Add 1 to years 13+, because no data for year 13 (2012)
# YA$YEAR_NUMBER %<>% ifelse(.>= 13, .+1, .)

yield_anomoly = merge(LA, YA, by='YEAR_NUMBER', all=TRUE)

# Analyze and plot
yy = 'LOCAL_ADAPTATION'
xx = 'YIELD_ANOMOLY'
ee = 'ERROR'
dd = yield_anomoly

ff = paste(yy, '~', xx) %>% 
  as.formula

mod = lm(ff, data=dd)
qqPlot(residuals(mod))
Anova(mod)
summary(mod)

plt = ggplot() +
  xlab("Yield Anomoly [bu/ac]") +
  ylab("Home Advantage [bu/ac]") +
  geom_errorbar(data=dd,
                aes_string(x=xx,
                           ymin=paste(yy, ee, sep="-"),
                           ymax=paste(yy, ee, sep="+")),
                width=0.2,
                color="#666666") +
  geom_hline(aes(yintercept=mean(dd[, yy])),,
             linetype='dashed',
             color="#666666") +
  # geom_smooth(aes_string(x=xx,
  #                        y=yy),
  #             color='tomato4',
  #             method='lm',
  #             weight=0.1,
  #             linetype="dashed") +
  geom_point(data=dd,
             aes_string(x=xx,
                        y=yy),
             color='steelblue3',
             size=4) +
  theme_light()
plt

png(file.path(dir_in, 'home_field_advantage_vs_yield_anomoly.png'), height=4, width=5, units='in', res=300)
plt
dev.off()
#### Genotype Home Advantage ~ Yield Potential ####
mean_yld = aggregate(YIELD ~ GENOTYPE, 
                     data=checks, 
                     function(x) mean(x, na.rm=TRUE))
home_yld = aggregate(YIELD ~ GENOTYPE, 
                     data=subset(checks, IS_HOME==TRUE), 
                     function(x) mean(x, na.rm=TRUE))
names(home_yld) = c('GENOTYPE', 'HOME_YIELD')

# Genetic potential of each genotype, accounting for year.

genetic_potential = lm(YIELD ~ YEAR + GENOTYPE,# + LOCATION + YEAR:LOCATION, 
                       data=checks) %>% 
  coefficients %>% 
  sapply(function(x) ifelse(is.na(x), 0, x))
# save intercept as the expected genetic potential
mean_potential = genetic_potential['(Intercept)']
# select genotype coefficients
tt = grepl("GENOTYPE", names(genetic_potential))
genetic_potential = genetic_potential[tt] %>% 
  add(mean_potential)
names(genetic_potential) %<>% gsub("GENOTYPE", "", .)
genetic_potential %<>% data.frame(
  GENOTYPE = names(.),
  GENETIC_POTENTIAL = .
) 

home_yld = merge(home_yld, genetic_potential, by='GENOTYPE')
home_yld$YIELD_BUMP = home_yld$HOME_YIELD - home_yld$GENETIC_POTENTIAL

xx = 'GENETIC_POTENTIAL'
yy = 'YIELD_BUMP'
dd = home_yld

ff = paste(yy, xx, sep="~") %>% 
  as.formula

mod = lm(ff, data=home_yld)
qqPlot(residuals(mod))
Anova(mod)
summary(mod)

xyplot(ff, data=dd, pch=20, cex=2, alpha=0.5)

#### Home advantage vs variance by variety ####
# This is broad vs local adaptation. A small variance gives broad adaptation.
# Algorithm:
#   For each variety:
#     Center
#     Calculate SD/Var
#     Extract "Home" value
#     Save in DF

rdf = lm(YIELD ~ LOCATION*YEAR, data=checks) %>% 
  residuals %>% 
  cbind(checks[, c('IS_HOME', 'GENOTYPE')]) %>% 
  as.data.frame
names(rdf)[1] = 'RESID'

adapt_type = split(rdf, rdf$GENOTYPE) %>% 
  lapply(function(x) {
    x$RESID %<>% scale(scale=FALSE) # center
    vv = var(x$RESID, na.rm=TRUE) %>%
      round(3)
    hm = x[x$IS_HOME, 'RESID'] %>%
      mean(na.rm=TRUE) %>%
      round(3)
    out = c(hm, vv)
    return(out)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame
names(adapt_type) = c('HOME_ADVANTAGE', 'VARIANCE')


xyplot(HOME_ADVANTAGE ~ sqrt(VARIANCE), data=adapt_type, type=c('r','p'))

mod = lm(HOME_ADVANTAGE ~ sqrt(VARIANCE), data=adapt_type)
summary(mod)
Anova(mod)

w = FALSE
# w = TRUE
if (w) {
  fout = paste("PLOT DATA - Home Advantage vs Variance ",
               Sys.Date(),
               ".csv", sep="") %>% 
    file.path(dir_in,  .)
  write.csv(adapt_type, fout, row.names=FALSE)
}

plt = ggplot(adapt_type) +
  xlab("Variation Across Environments [bu/ac]") +
  ylab("Home Advantage [bu/ac]") +
  geom_smooth(aes(x=sqrt(VARIANCE),
                  y=HOME_ADVANTAGE),
              color='tomato4',
              method='lm',
              weight=0.1,
              linetype="dashed") +
  geom_point(aes(x=sqrt(VARIANCE),
                 y=HOME_ADVANTAGE),
             color='steelblue3',
             size=2,
             alpha=0.5) +
  theme_light()
plt


png(file.path(dir_in, 'home_field_advantage_vs_yield_variability.png'), height=4, width=5, units='in', res=300)
plt
dev.off()


#### The Fitness Set ####
checks$GROUP = ifelse(checks$REGION=='south', 'SOUTHERN', 'NORTHERN')

regions = aggregate(YIELD ~ GENOTYPE + GROUP + YEAR, data=checks, function(x) mean(x, na.rm=TRUE)) %>% 
  dcast(GENOTYPE + YEAR ~ GROUP, value.var='YIELD')
regions_se = aggregate(YIELD ~ GENOTYPE + GROUP + YEAR, data=checks, function(x) sd(x, na.rm=TRUE)/sqrt(length(x-1))) %>% 
  dcast(GENOTYPE + YEAR ~ GROUP, value.var='YIELD')
xyplot(NORTHERN ~ SOUTHERN, group=YEAR, data=regions)

fitsets = ggplot(regions) +
  scale_y_continuous(breaks=scales::pretty_breaks(3)) +
  scale_x_continuous(breaks=scales::pretty_breaks(3)) +
  facet_wrap(~YEAR,
             scales='free', 
             ncol=6) +
  geom_point(aes(x=SOUTHERN/1000,
                 y=NORTHERN/1000),
             alpha=0.5,
             color='steelblue4') +
  xlab(expression(paste("Yield in Southern Derived Region [Mg ", "ha"^"-1", "]"))) +
  ylab(expression(paste("Yield in Northern Derived Region [Mg ", "ha"^"-1", "]"))) +
  theme_light()
fitsets

w = FALSE
w = TRUE
if(w) {
  fout = paste("PLOT DATA - Fitness Sets ",
               Sys.Date(),
               ".csv",
               sep="") %>% 
    file.path(dir_in, "Figures", .)
  write.csv(regions, fout, row.names=FALSE)
}


#### Locally adapted is also max? ####
tt = split(checks, checks$LOCYR) %>% 
  lapply(function(x) {
    x$IS_MAX = ifelse(x$YIELD == max(x$YIELD), TRUE, FALSE)
    return(x)
  }) %>% 
  unsplit(checks$LOCYR)
tt$RANK_TYPE = 1*tt$IS_HOME + 2*tt$IS_MAX + 1
# is home, is max
# 1 = false/false
# 2 = true/false
# 3 = false/true
# 4 = true/true
RT = c('Neither', 'Home', 'Max', 'Both')
tt$RANK_TYPE %<>% sapply(function(x) return(RT[x])) %>% 
  factor(levels = RT)


ggplot(subset(tt, RANK_TYPE != 'Neither')) +
  ylim(c(0, 200)) +
  geom_bar(aes(x=RANK_TYPE)) # gives counts

RANK_COUNTS[4] / sum(RANK_COUNTS[c(2,4)])
RANK_COUNTS[4] / sum(RANK_COUNTS[c(3,4)])

