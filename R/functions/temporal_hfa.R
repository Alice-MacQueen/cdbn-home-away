#' @title calculate_temporal_hfa <internal>
#' 
#' @description calculate hfa for each year
#' 
#' @param x data.frame
#' @param ff formula
#' 
#' @details the formula should contain a term, <year>_num:is_home, and also specify a response (e.g. pheno)
#' 
#' @return a data.frame with the coefficients plus relevant 

.calculate_temporal_hfa = function(x, ff, year) {
  y = as.character(ff)[2]
  
  # get the coefficients
  home_coef = lm(ff, x) %>% 
    summary %>% 
    coef
  
  # format
  home_coef %<>% 
    rownames %>% 
    grepl(':is_homeTRUE', .) %>% 
    home_coef[., ]
  
  rownames(home_coef) %<>% 
    gsub(':is_homeTRUE', '', .) %>% 
    gsub(year, '', .)

  # return values
  home_coef %<>% 
    data.frame(year = rownames(.),
               year_num = type.convert(rownames(.)),
               .,
               p.adj = p.adjust(home_coef[, 'Pr(>|t|)']),
               stringsAsFactors=FALSE) %>% 
    set_rownames(NULL)
  colnames(home_coef) %<>% 
    gsub('year', year, .) %>% 
    gsub('Pr...t..', 'p.value', .) %>% 
    gsub('Std..Error', 'Std.Error', .)
  
  
  return(home_coef)
}

#' @title temporal_hfa
#' 
#' @description calculate trend in home field advantage across time. The formula is:
#' 
#' `pheno ~ year + site + geno + year:site + year:is_home`. 
#' 
#' The year:is_home coefficient is the home field advantage in each year. 
#' 
#' @param data data.frame
#' @param site column in `data` indicating spatial environment
#' @param year column in `data`indicating temporal environment
#' @param geno column in `data` indicating genotype
#' @param pheno column in `data` indicating phenotype
#' @param popn column in `data` indicating groups of genotypes
#' @param blup_home identify home site using blup?

temporal_hfa = function(data, site, year, geno, pheno, popn=NA, 
                        blup_home=TRUE, 
                        formula_override=NA,
                        parallel=TRUE) {
  require(magrittr)
  require(parallel)
  require(car)
  
  if (grepl('mingw', version$os)) parallel=FALSE
  ncpu = ifelse(parallel, detectCores(), 1)
  
  # set up dataframe
  
  if (is.na(popn)) {
    dd = list(data)
  } else {
    dd = split(data, data[, popn])
  }
  
  # formula
  if (is.na(formula_override)) {
    ff = "pheno ~ year + site + geno + year:site + year:is_home"  %<>%
      gsub('pheno', pheno, .) %>% 
      gsub('site', site, .) %>% 
      gsub('year', year, .) %>% 
      gsub('geno', geno, .) %>% 
      formula
  } else {
    ff = formula_override
  }
  
  # has home already been identified?
  id_home_site = colnames(data) %>% 
    grepl('is_home', .) %>% 
    any %>% 
    not
  if (id_home_site) {
    dd %<>% mclapply(function(x) {
      id_home(x, site, year, geno, pheno, blup=blup_home, verbose=FALSE)
    }, 
    mc.cores=ncpu)
  }
  
  # calculate annual hfa
  out = mclapply(dd, .calculate_temporal_hfa, ff, year,
                 mc.cores=ncpu) %>% 
    do.call(rbind, .)
  
  # format output
  # population, if specified
  if (!is.na(popn)) {
    popn_id = rownames(out) %>% 
      strsplit('\\.') %>% 
      sapply('[', 1)
    
    out = data.frame(
      popn = popn_id,
      out,
      stringsAsFactors=FALSE
    ) %>% 
      set_rownames(NULL)
    colnames(out) %<>% gsub('popn', popn, .)
  }
  
  # year/temporal
  if (is.factor(data[, year])) {
    out[, year] %<>% factor(levels=levels(data[, year]))
  }
  
  # run anova
  year_num = paste0(year, '_num')
  ff_aov = paste('Estimate', year_num, sep='~')
  if (!is.na(popn)) {
    ff_aov %<>% paste('*', popn)
  }
  
  out_lm = ff_aov %>% 
    formula %>% 
    lm(data=out)
  out_aov = Anova(out_lm)
  
  out = list(anova = out_aov,
             model = out_lm,
             temporal_hfa = out,
             formulas = c(hfa = ff,
                        anova = ff_aov))
  
  return(out)
}
