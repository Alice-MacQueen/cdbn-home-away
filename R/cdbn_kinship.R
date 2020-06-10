#' @title Find a kinship matrix using the van Raden method.
#'
#' @description Calculate the kinship matrix using code from GAPIT and bigsnpr
#'     and the methods of VanRaden (2009, J. Dairy Sci. 91:4414???C4423). This
#'     code is a modified version of the same code in GAPIT. Note
#'     that this matrix cannot currently be used with the GWAS methods in
#'     bigsnpr; however, this matrix could be used for other analyses.
#'
#' @param snp A "bigSNP" object; load with bigsnpr::snp_attach().
#' @param ind.row An integer vector of the rows (individuals) to find a
#'     kinship matrix for. Defaults to all rows.
#' @param hasInbred Logical. Does the SNP file contain inbred individuals or
#'     closely related individuals, like siblings? Default is TRUE.
#' @param saveoutput Logical. Should the output be saved to the working
#'     directory?
#'
#' @return A kinship matrix with labeled rows and columns.
#'
#' @import bigsnpr
#' @import bigstatsr
#'
#' @examples
#' \dontrun{
#' K <- cdbn_kinship(snp = snp, saveoutput = TRUE)
#' }
#'
#' @export
cdbn_kinship <- function(snp, ind.row = NA, hasInbred = TRUE,
                          saveoutput = FALSE){
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if(saveoutput == FALSE){
    message(paste0("'saveoutput' is FALSE, so the kinship matrix will not ",
                   "be saved to the working directory."))
  }
  
  G <- snp$genotypes
  if(!is.na(ind.row[1])){
    nInd <- length(ind.row)
  } else {
    nInd <- snp$genotypes$nrow
    ind.row <- 1:nInd
  }
  # Centered (mean of rows subtracted) transposed crossproduct of snp file.
  K <- big_tcrossprodSelf(G, ind.row = ind.row,
                          fun.scaling = big_scale(center = TRUE,
                                                  scale = FALSE))
  
  #Extract diagonals
  i = 1:nInd
  j = (i - 1)*nInd
  index = i + j
  d = K[index]
  DL = min(d)
  DU = max(d)
  floor = min(K[i, i])
  
  K = (K[i, i] - floor)/(DL - floor)
  MD = (DU - floor)/(DL - floor)
  
  if(!is.na(ind.row[1])){
    rownames(K) <- snp$fam$sample.ID[ind.row]
    colnames(K) <- snp$fam$sample.ID[ind.row]
  } else {
    rownames(K) <- snp$fam$sample.ID
    colnames(K) <- snp$fam$sample.ID
  }
  
  if(MD > 2){
    K[index] <- K[index]/(MD-1)+1
  }
  #Handler of inbred
  if(MD < 2 & hasInbred){
    K = 2*K / ((DU - floor) / (DL - floor))
  }
  
  if(saveoutput){
    saveRDS(K, paste0("Kinship_van_Raden_method_", nInd, "_individuals_",
                      ".rds"))
  }
  return(K)
}

#' @title Wrapper for the snp_autoSVD function for the CDBN.
#'
#' @description This is a wrapper to determine population structure for GWAS
#'     for a bigSNP object. Arguments that are recognized by 
#'     bigsnpr::snp_autoSVD can also be specified in this function.
#'
#' @param snp A "bigSNP" object; load with bigsnpr::snp_attach().
#' @param k Integer. The number of principal components to find. Default is 10.
#' @param ncores Integer. Number of cores to use. Default is one.
#' @param saveoutput Logical. Should the output be saved to the working
#'     directory?
#' @param ... Other arguments to \code{\link{snp_autoSVD}}.
#'
#' @return A big_SVD object.
#'
#' @import bigsnpr
#' @import bigstatsr
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when
#' @importFrom tibble enframe
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' svd20 <- cdbn_autoSVD(snp = snp, k = 20, saveoutput = TRUE)
#' }
#'
#' @export
cdbn_autoSVD <- function(snp, k = 10, ncores = 1, saveoutput = FALSE, ...){
  requireNamespace("dots") # devtools::install_github("lcolladotor/dots")
  fun.scaling <- dots::dots(name = 'fun.scaling', value = snp_scaleBinom(),
                            ...)
  thr.r2 <- dots::dots(name = 'thr.r2', value = 0.2, ...)
  size <- dots::dots(name = 'size', value = 100/thr.r2, ...)
  roll.size <- dots::dots(name = 'roll.size', value = 50, ...)
  int.min.size <- dots::dots(name = 'int.min.size', value = 20, ...)
  #alpha.tukey <- dots::dots(name = 'alpha.tukey', value = 0.05, ...)
  #min.mac <- dots::dots(name = 'min.mac', value = 10, ...)
  #max.iter <- dots::dots(name = 'max.iter', value = 5, ...)
  verbose <- dots::dots(name = 'verbose', value = FALSE, ...)
  
  # Argument error checks; Set up numeric chromosome data frame.
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  G <- snp$genotypes
  CHR <- snp$map$chromosome
  POS <- snp$map$physical.pos
  plants <- snp$fam$sample.ID
  
  if(saveoutput == FALSE){
    message(paste0("'saveoutput' is FALSE, so the svd will not be saved to ",
                   "the working directory."))
  }
  
  # Determine population structure
  svd <- snp_autoSVD(G = G, infos.chr = CHR, infos.pos = POS,
                     ncores = ncores, k = k, fun.scaling = fun.scaling,
                     thr.r2 = thr.r2, size = size, roll.size = roll.size,
                     int.min.size = int.min.size, #alpha.tukey = alpha.tukey,
                     #min.mac = min.mac, max.iter = max.iter,
                     verbose = verbose)
  if(saveoutput){
    saveRDS(svd, file = paste0("SVD_", length(plants), "g_", k, "PCs.rds"))
  }
  return(svd)
}
