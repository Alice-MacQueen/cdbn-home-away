.scale_pheno = function(x, pheno) {
  # for applying
  colname = paste0('rel_', pheno)
  x[, colname] = scale(x[, pheno])
  return(x)
}

#' @title .id_top_pheno
#' 
#' @description id the level of `x[, site]` where the mean relative value of 
#' `x[, pheno]` is highest. By default, uses shrinkage estimates of mean 
#' relative `x[, pheno]` (via `lme4`) if any site occurs in at least two 
#' years. Otherwise (or optionally) uses means.
#' 
#' @param x data.frame or matrix for a single unit, ex a genotype. Coerced to data.frame.
#' @param site character of grouping column, ex field site locations
#' @param pheno character of response, ex yield, biomass
#' @param blup logical - use shrinkage estimates?
#' @param verbose logical
#' 
#' @returns x, with an additional logical column `"Home"`, which indicates 
#' whether the value of `x[, site]` is the home site.

.id_top_pheno = function(x, site, pheno, blup=TRUE, verbose=TRUE) {
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
  
  # calculate mean pheno within site.
  if (blup) {
    require(lme4)
    
    if (verbose) message('Using lme4 to identify the home site')
    site_eff = paste(pheno, '~ 0 + (1|', site, ')') %>% 
      formula %>% 
      lmer(data=x,
           control=lmerControl(calc.derivs = FALSE)) %>% # tip from lme4 pheno vignette
      coef %>% 
      extract2(site)
    site_eff %<>% 
      .[,1] %>% 
      set_names(rownames(site_eff))
    
  } else {
    site_eff = tapply(x[, pheno], x[, site], mean, simplify=TRUE)  # faster than aggregate
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
#' @description Identify the home site based on highest relative phenotype value. 
#' 
#' @param df data.frame of performance data by site, year, and variety, for example. 
#' @param site character, column name indicating spatial units
#' @param year character, column name indicating temporal units
#' @param geno character, column name indicating genetic units
#' @param pheno character, column name indicating the phenotype (ex. yield)
#' @param blup logical (`TRUE`). Use shrinkage estimates of mean relative performance within each site across years?
#' @param verbose logical (`TRUE`).
#' 
#' @details The home site for a given variety is defined as the location where 
#' variety performs best across years, relative to other varieties. It's calculated
#' by:
#' 
#' 1. Calculate relative phenotype within site-year by scaling (mean = 0) and 
#' scaling (st. dev. = 1) phneotype across varieties within site-years
#' 2. Find the expected relative phenotype for each variety within
#' within sites across years. This is either the raw mean or the shrunk/BLUP mean (via
#' a random intercept model using `lme4`). 
#' 3. For each variety, identify the site with the highest expected relative phenotype value. 
#' 
#' @returns a data.frame with two new columns:
#' 
#' - `'rel_<pheno>'`: numeric, the relative phenotype value of each variety within site-year
#' - `is_home`: logical of whether a site is a variety's home site.

id_home = function(df, site, year, geno, pheno, blup=TRUE, verbose=FALSE) {
  require(magrittr)
  require(parallel)
  
  rel_colname = paste0('rel_', pheno)
  
  # make site-year vector
  site_year = paste(df[, site], df[, year], sep='_')
  # center/scale performance within site-year
  df %<>% 
    split(site_year) %>% 
    lapply(.scale_pheno, pheno) %>% 
    do.call(rbind, .)
  
  # find highest relative yield for each genotype
  df %<>%
    split(df[geno]) %>%
    lapply(.id_top_pheno, site, rel_colname, blup, verbose) %>% 
    do.call(rbind, .)
  
  return(df)
}
