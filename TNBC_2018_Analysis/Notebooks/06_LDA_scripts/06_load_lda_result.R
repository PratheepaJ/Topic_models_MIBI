#' Load existing LDA estimations and extract posterior samples.
#' 
#' @param file string. File path for the existing LDA estimations
#' 
#' @return list. Large list contains parameters theta and beta estimations
#' @importFrom rstan extract
#' 

load_lda_result <- function(file) {
  # load lda file as stan.fit
  #load(here::here("Output", "Data", "05_LDA_workflow", paste0(file)))
  load(file)
  
  # extract 
  res <- rstan::extract(
    stan.fit,
    permuted = TRUE,
    inc_warmup = FALSE,
    include = TRUE
  )
  

  return(list(stan.fit,
              res))
}