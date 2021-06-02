#' @title try to load a library, and download if needed
#' 
#' @description mostly helps resolve packrat bugs. Also installs in parallel.
#' 
#' @param lib_string charracter (vector) with package names
#' @param cores integer or logical - install packages in parallel? TRUE uses
#' `parallel::detectCores()`
#' 
#' @returns nothing, but loads packages into the user environment.

get_libraries = function(lib_string, in_parallel=TRUE) {
  
  if (is.logical(in_parallel)) {
    tt = in_parallel & ('parallel' %in% rownames(installed.packages()))
    if (tt) {
      threads = parallel::detectCores() 
    } else {
      threads = 1
    }
    
  } else {
    threads = in_parallel
  }
  
  if (!require(lib_string, character.only=TRUE)) {
    msg = paste('Trying to install', lib_string)
    message(msg)
    
    install.packages(lib_string, Ncpus=threads)
    if (!require(lib_string, character.only=TRUE)) {
      msg = paste("Could not install", lib_string)
      stop(msg)
    }
  }
}
