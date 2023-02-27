#' Align ASVs proportions in topics across chain.
#'
#' @param beta array. three dimensional array (iterations * topic * ASVs).
#' @param aligned matrix. A matrix of topic assignment for each chain. output of alignmentMatrix().
#' @inheritParams alignmentMatrix
#'
#' @return array. three dimensional array.
#' @importFrom abind abind
#' @export
#'

betaAligned <- function(
  beta,
  K,
  aligned,
  warm_up_iter = NULL,
  iter = 2000,
  chain = 2){
  
  # determine the iteration used in posterior sampling (subtract warm up iterations)
  if (is.null(warm_up_iter)) {
    iterUse = iter / 2
  } else {
    iterUse = iter - warm_up_iter
  }
  
  # align topics between chains
  # no need to switch array dimnension
  beta_chain <- list()
  beta_chain[[1]] <- beta[1:(iterUse), ,]
  beta_chain[[1]] <- beta_chain[[1]][, aligned[,1], ]
  
  beta_aligned <- beta_chain[[1]]

  for(ch in 2:chain){
    beta_chain[[ch]] <- beta[((ch-1)*(iterUse)+1):(ch*(iterUse)),,]
    beta_chain[[ch]] <- beta_chain[[ch]][,aligned[,ch],]
    
    beta_aligned <- abind::abind(
      beta_aligned,
      beta_chain[[ch]],
      along = 1)
  }

  # beta_aligned <- abind::abind(
  #   beta_chain[[1]],
  #   beta_chain[[2]],
  #   along = 1)
  # for(ch in 3:chain){
  #   beta_aligned <- abind::abind(
  #     beta_aligned,
  #     beta_chain[[ch]],
  #     along = 1)
  # }

  return(beta_aligned)
}
