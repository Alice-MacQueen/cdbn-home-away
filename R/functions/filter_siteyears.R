#' @title filter_siteyears
#' 
#' @description Subset locations by a minimum number of occurrances
#' 
#' @param x data.frame
#' @param site column containing trial location
#' @param year column containing year
#' @param min_times numeric, minimum number of years a site must be present. Default 3. 
#' 
#' @return a data.frame, a subset of x. 

filter_siteyears = function(x, site, year, min_times=3) {
  require(magrittr)
  
  common_sites = c(site, year) %>% 
    x[, .] %>% 
    unique %>%
    .[[site]]
  
  common_sites %<>%
    tapply(., ., length) %>%  # count occurrences
    .[!is.na(.)]
  
  common_sites %<>%  # filter by min_times
    .[.>=min_times]
    

  test = x[, site] %in% names(common_sites)
  out = x[test, ]
  
  return(out)
}
