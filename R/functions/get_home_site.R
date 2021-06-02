#' @title get_home_site
#' 
#' @description extract home location information for each genotype
#' 
#' @param data data.frame from call to `id_home()`
#' @param geno column in `data` with genotype
#' @param site column in `data` with location environment
#' 
#' @return a `data.frame` with columns `geno`, `site`, and `home_years` - the number
#' of years a genotype appeared at a home site.

get_home_site = function(data, geno, site) {
  is_home = data[, 'is_home']
  out = data[is_home, c(geno, site)]
  counts = data[, geno] %>% 
    tapply(., ., length)
  
  out %<>% unique
  out[, 'home_years'] = counts[out[, geno]]
  return(out)
}


#' @title measure_home_distance
#' 
#' @description calculate the distance between home sites for among genotypes
#' 
#' @param data data.frame from a call to `id_home()`
#' @param locations a SpatialPointsDataFrame (recommended) or a data.frame with columns `site`, `lat`, and `long`.
#' @param geno genotype column in `data` and in `locations`
#' @param site spatial environment column in `data`
#' @param lat latitude/y coordinate column in `locations`
#' @param long longitude/x coordinate column in `locations`
#' @param great_circle are coordinates on an ellipse? 
#' 
#' @return a distance matrix. row/column names are genotypes. Values are distances
#' among between home sites of the genotype pairs.
#' 
#' @note If lat/long are in degrees, then set `great_circle` to TRUE. If lat/long
#' are in linear units, set false. Can also pass a SpatialPointsDataFrame, which will
#' override this the `great_circle` setting and coordinates will be derived from the object. 
#' 
#' @seealso package `"sp"`

measure_home_distance = function(data, locations, geno, site, lat=NA, long=NA, great_circle=TRUE) {
  require(sp)
  require(magrittr)
  
  if (grepl("SpatialPoints", class(locations))) {
    distances = spDists(locations)
  } else {
    distances = locations[, c(long, lat)] %>% 
      as.matrix %>% 
      spDists(longlat=great_circle)
  }
  
  colnames(distances) = locations[[site]]
  rownames(distances) = locations[[site]]
  
  data %<>% 
    get_home_site(geno, site)
  
  vars = data[, geno] %>% 
    as.character
  homes = data[, site] %>%
    as.character %>% 
    set_names(vars)
  nvar = length(vars)
  
  out = matrix(nrow = nvar,
               ncol = nvar) %>% 
    set_rownames(vars) %>% 
    set_colnames(vars)
  
  for (i in vars) {
    home_i = homes[i]
    for (j in vars) {
      home_j = homes[j]
      out[i, j] = distances[home_i, home_j]
    }
  }
  
  return(out)
}
