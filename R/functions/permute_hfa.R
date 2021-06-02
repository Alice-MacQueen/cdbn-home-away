#' @title .quick_resid <internal>
#' 
#' @description calculate residuals of a linear model as efficiently as possible (using base R).
#' Residuals can be used in subsequent regressions as partial effects.
#' 
#' @param ff formula
#' @param data data.frame containing formula
#' 
#' @return A vector of residuals, length = nrow(data)

.quick_resid = function(ff, data) {
  require(Matrix)

    X = model.matrix(ff, data) 
  X = qr(X)
  Y = model.frame(ff, data)[, 1]
  
  out = qr.resid(X, Y)

  return(out)
}


#' @title .generate_sets <internal>
#' 
#' @description generate permutation sets structured by site-year
#' 
#' @param x data.frame
#' @param site column identifying spatial environment
#' @param year column identifying temporal environment
#' @param times number of sets to produce. 
#' 
#' @return a data.frame of permutation sets - row numbers used to reorder x. The first column is the original data.

.generate_sets = function(x, site, year, times, seed=NULL) {
  require(magrittr)
  require(permute)
  
  control = paste(x[, site], x[, year], sep='_') %>% 
    as.factor
  N = nrow(x)
  
  set.seed(seed)
  ss = shuffleSet(N, times, control=how(blocks=control)) %>% # permutations
    rbind(seq_len(N),  # include original data
          .) %>% 
    t %>% 
    data.frame
  
  return(ss)
}

#' @title .calculate_hfa <internal>
#' 
#' @description calculate the home field advantage based on the formula, ff. Uses efficient
#' packages for fastest permutation results.
#' 
#' @param x data.frame
#' @param ff formula to use. part_pheno ~ geno */+ is_home
#' 
.calculate_hfa = function(x, ff, part_pheno, geno, year, LAPACK=TRUE) {
  require(Matrix)
  
  X = droplevels(x) %>%   # make sure all factors in ff are present for LAPACK
    model.matrix(ff, data=.) %>% 
    qr(LAPACK=LAPACK)
  
  Y = x[, part_pheno, drop=FALSE]
  
  home_coef = qr.coef(X, Y) %>% 
    as.matrix
  
  is_home = grepl('is_homeTRUE', rownames(home_coef))
  home_coef %<>% .[is_home, , drop=FALSE]
  
  rownames(home_coef) %<>% 
    gsub(':is_homeTRUE', '', .) %>% 
    gsub(geno, '', .) %>% 
    gsub(year, '', .) %>% 
    gsub('TRUE', '', .)
  
  return(home_coef)
}


#' @title .two_tailed <internal>
#' 
#' @description perform two-tailed test on permutation values
#' 
#' @param x matrix. Each column is the result of a different permutation. The first column is the original data.
#' 
#' @return a vector of p-values for each row in x.

.two_tailed = function(x) {
  require(magrittr)
  
  alpha = sweep(x, 1, x[, 1], '>=') %>% 
    rowSums %>% 
    divide_by(ncol(x))
  p_val = cbind(alpha, 1-alpha) %>% 
    apply(1, min) %>% 
    multiply_by(2)
  return(p_val)
}

#' @title .calculate_intervals <internal>
#' 
#' @description calculate median and 90% confidence intervals of difference of permutation and observed data
#' 
#' @param x permutation results. Rows are values, columns are permutations. The first column is observed.
#' 
#' @return a matrix of columns "median", "q05" (5th percentile), "q95" (95th percentile) of (observed - permutation)
.calculate_intervals = function(x) {
  require(magrittr)
  
  difference = sweep(x, 1, x[, 1], '-') %>%
    multiply_by(-1)
  
  out = cbind(
    median = apply(difference, 1, median, na.rm=TRUE),
    p05 = apply(difference, 1, quantile, 0.05, na.rm=TRUE),
    p95 = apply(difference, 1, quantile, 0.95, na.rm=TRUE)
  )
  
  return(out)
}

#' @title permute_hfa
#' 
#' @description test the magnitude and significance of the home field advantage
#' versus what would be expected by chance.
#' 
#' @param data data.frame
#' @param level parameter at which to apply the home field advantage. Default is population.
#' @param site column name containing spatial environmental information
#' @param year column name containing temporal environmental information
#' @param geno column name containing genotype or variety information
#' @param pheno column name containing phenotype or performance information
#' @param popn column name differnetiating sub-populations of genotypes within the dataset.
#' 
#' @return 
#' 
permute_hfa = function(data,
                       level=c('population', 'genotype', 'year', 'site'),
                       site=NA, 
                       year=NA, 
                       geno=NA, 
                       pheno=NA, 
                       popn=NA,
                       times=99, 
                       blup_home=TRUE,
                       parallel=TRUE,
                       seed=NULL) {
  require(magrittr)
  require(parallel)  # will need to transition to foreach for compatibility
  
  ncpu = ifelse(parallel, detectCores(), 1)
  level = match.arg(level)
  LAPACK = !(level %in% c('site', 'year'))  # cannot guarantee a home site in each year, or for each variety.
  # LAPACK=TRUE
  
  # new column names
  part_pheno = paste0('part_', pheno)
  rel_pheno = paste0('rel_', pheno)
  
  # formula for calculating HFA effect
  ff = switch(level,
              'population' = ' ~ geno + is_home',  # overall HFA
              'genotype' = ' ~ geno + geno:is_home',
              'year' = '~ geno + year:is_home',
              'site' = '~ geno + site:is_home')  # HFA for each genotype
  ff %<>%
    gsub('geno', geno, .) %>% 
    gsub('year', year, .) %>% 
    gsub('site', site, .) %>% 
    formula
  
  # select and format data into list of dataframes for each populations
  dd = 
    c(site, year, geno, pheno, popn) %>% 
    na.omit %>% 
    data[, .]
  
  if (! (is.na(popn))) {
    dd %<>% split(dd[, popn])
  } else {
    dd %<>% list
  }
  
  
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
  
  # calculate partial and relative yields
  dd %<>% mclapply(function(x) {
    
    # calc rel_pheno and is_home
    x %<>% id_home(site, year, geno, pheno, blup_home, verbose=FALSE)

    # Partial phenotype based on site-year. Relatively slow. 
    x[, part_pheno] = 
      paste(pheno, '~', year, '*', site) %>% 
      formula %>% 
      .quick_resid(data=x) %>% 
      scale(scale=FALSE)
    
    return(x)
  }, mc.cores=ncpu)
  
  # permute HFA within each population
  results = lapply(dd, function(x) {
    # Set up structured permutations
    sets = .generate_sets(x, site, year, times, seed)
    
    # Permute HFA
    coef_permute = mclapply(sets, function(ss) {
      # ID home site
      x[, c(rel_pheno, part_pheno)] %<>% .[ss, ]  # permute phenotypes within site-year
      x %<>%
        split(x[, geno]) %>%
        lapply(.id_top_pheno, site, rel_pheno, blup=blup_home, verbose=FALSE) %>%  # This is the rate-limiting step.
        do.call(rbind, .)
      
      # calculate HFA using qr decomposition and lapack
      home_coef = .calculate_hfa(x, ff, part_pheno, geno, year, LAPACK)
      
      return(home_coef)
    }, mc.cores=ncpu) %>% 
      do.call(cbind, .)
    
    colnames(coef_permute) = c('observed', 
                               paste0('perm', 1:(ncol(coef_permute)-1)))
    
    # Calculate p-values and effects
    test = cbind(
      intervals = .calculate_intervals(coef_permute),
      p_val = .two_tailed(coef_permute)
    )
    
    results = list(results = test,
                   perms = coef_permute)
    
    return(results)
  })
  
  if (level == 'population') {
    test_results = lapply(results, function(x) x$results) %>% 
      do.call(rbind, .) %>% 
      set_rownames(names(results))
    perms = lapply(results, function(x) x$perms) %>% 
      do.call(rbind, .) %>% 
      set_rownames(names(results))
  } else {
    lvl = switch(level, 
                 'genotype' = geno,
                 'year' = year,
                 'site' = site)
    test_results = names(results) %>% 
      lapply(function(x) {
      res = results[[x]]$results
      res = data.frame(popn = x,
                       lvl = rownames(res),
                       res,
                       stringsAsFactors=FALSE) %>% 
        set_rownames(NULL)
      names(res) %<>% 
        gsub('lvl', lvl, .) %>% 
        gsub('popn', popn, .)
      return(res)
    }) %>% 
      do.call(rbind, .)
    
    perms = lapply(results, function(x) x$perms)
  }
  
  out = list(home_field = test_results,
             perms = perms)

  return(out)
}

