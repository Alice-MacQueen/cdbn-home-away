.scale_performance = function(x, performance) {
  # for applying
  colname = paste0('rel_', performance)
  x[, colname] = scale(x[, performance])
  return(x)
}

#' @title .id_best_performance
#' 
#' @description id the level of `x[, site]` where the mean relative value of 
#' `x[, performance]` is highest. By default, uses shrinkage estimates of mean 
#' relative `x[, performance]` (via `lme4`) if any site occurs in at least two 
#' years. Otherwise (or optionally) uses means.
#' 
#' @param x data.frame or matrix for a single unit, ex a genotype. Coerced to data.frame.
#' @param site character of grouping column, ex field site locations
#' @param performance character of response, ex yield, biomass
#' @param blup logical - use shrinkage estimates?
#' @param verbose logical
#' 
#' @returns x, with an additional logical column `"Home"`, which indicates 
#' whether the value of `x[, site]` is the home site.

.id_best_performance = function(x, site, performance, blup=TRUE, verbose=TRUE) {
  require(magrittr)
  
  x %<>% as.data.frame(stringsAsFactors=FALSE)
  
  x[, site] %<>%
    factor %>%
    droplevels
  
  # find max years at any one site. 
  if (blup) {
    site_years = x[, site] %>%
      tapply(., ., length) %>%
      max
    if (site_years < 3) { #(site_years < 2 & blup) {#
      if (verbose) message('Cannot blup home field. Using (unshrunk) mean values')
      blup = FALSE
    }
  }
  
  # calculate mean performance within site.
  if (blup) {
    require(lme4)
    
    if (verbose) message('Using lme4 to identify the home site')
    site_eff = paste(performance, '~ 0 + (1|', site, ')') %>% 
      formula %>% 
      lmer(data=x,
           control=lmerControl(calc.derivs = FALSE)) %>% # tip from lme4 performance vignette
      coef %>% 
      extract2(site)
    site_eff %<>% 
      .[,1] %>% 
      set_names(rownames(site_eff))
    
  } else {
    site_eff = tapply(x[, performance], x[, site], mean, simplify=TRUE)  # faster than aggregate
  }
  
  # identify site with max mean
  is_max = site_eff %>% 
    which.max
  max_site = names(site_eff)[is_max]
  
  x$is_home = x[, site] == max_site
  
  return(x)
}

#' @title id_home
#' 
#' @description Identify the home site based on best relative performance. 
#' 
#' @param df data.frame of performance data by site, year, and variety, for example. 
#' @param site character, column name indicating spatial units
#' @param year character, column name indicating temporal units
#' @param variety character, column name indicating genetic units
#' @param blup logical (`TRUE`). Use shrinkage estimates of mean relative performance within each site across years?
#' @param verbose logical (`TRUE`).
#' 
#' @details The home site for a given variety is defined as the location where 
#' variety performs best across years, relative to other varieties. It's calculated
#' by:
#' 
#' 1. Calculate relative performance within site-year by scaling (mean = 0) and 
#' scaling (st. dev. = 1) performance across varieties within site-years
#' 2. Find the expected relative performance for each variety within
#' within sites across years. This is either the raw mean or the shrunk/BLUP mean (via
#' a random intercept model using `lme4`). 
#' 3. For each variety, identify the site with the highest expected relative performance.
#' 
#' @returns a data.frame with two new columns:
#' 
#' - `'rel_<performance>'`: numeric, the relative performance of each variety within site-year
#' - `is_home`: logical of whether a site is a variety's home site.

id_home = function(df, site, year, variety, performance, blup=TRUE, verbose=TRUE) {
  require(magrittr)
  require(parallel)
  
  rel_colname = paste0('rel_', performance)
  
  # make site-year vector
  site_year = paste(df[, site], df[, year], sep='_')
  # center/scale performance within site-year
  df %<>% 
    split(site_year) %>% 
    lapply(.scale_performance, performance) %>% 
    do.call(rbind, .)
  # is_perf = colnames(df) == performance
  # colnames(df)[is_perf] = 'standard_perf'
  
  # find highest relative yield for each genotype
  df %<>%
    split(df[variety]) %>%
    lapply(.id_best_performance, site, rel_colname, blup, verbose) %>% 
    # mclapply(.id_best_performance, site, rel_colname, blup, verbose,
    #          mc.cores=detectCores()) %>%
    do.call(rbind, .)
  
  return(df)
}
