#' aligned topics across chains
#'
#' @param theta array. three dimensional array (iterations * samples * topic).
#' @param spe SpatialExperiment object.
#' @param K integer. number of topics.
#' @param iterUse integer. number of iterations used in each chain.
#' @param chain integer. number of chains used in MCMC.
#' @param SampleID_name character. variable name that contains sample IDs.
#'
#' @return matrix. A matrix of topic assignment for each chain.
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join filter select
#' @import SpatialExperiment
#' @importFrom reshape2 melt
#' @importFrom stats cor
#' @export
#'

alignmentMatrix <- function(theta,
                            spe,
                            K,
                            iterUse = 2000,
                            chain = 4,
                            SampleID_name = "sample_id"){
  
  Chain <- Topic <- topic.dis <- NULL
  # theta is a 3 dimensional array (iterations * samples * topic)
  
  ## acquire sample identity
  dimnames(theta)[[2]] <- unique(spe$sample_id)
  ## acquire topic number
  dimnames(theta)[[3]] <- c(paste0("Topic_", seq(1,K)))
  
  # array to a dataframe
  ## melt() takes wide-format data and melts into long-format data
  theta_all = reshape2::melt(theta)
  colnames(theta_all) = c("iteration", "Sample", "Topic", "topic.dis")
  theta_all$Chain = paste0("Chain ", rep(seq(1, chain), each = (iterUse/2)))
  theta_all$Topic = factor(theta_all$Topic)
  theta_all$Chain = factor(theta_all$Chain)
  theta_all$Sample = as.character(theta_all$Sample)
  #theta_all$Sample = factor(theta_all$Sample)
  
  # join the SpatialExperimet object column data to theta_all
  ## colData gets the metadata (clinical data)
  # sam = (colData(spe) 
  #        |> data.frame()
  # )
  
  sam = (colData(spe)[1:9]
         |> unique()
         |> data.frame()
  )
  
  theta_all = dplyr::left_join(
    theta_all,
    sam,
    by = c("Sample"= SampleID_name)
  )
  
  
  # align matrix
  aligned <- matrix(nrow = K, ncol = chain)
  aligned[, 1] <- seq(1, K)
  corrTop <- numeric()
  for(j in 1:K) { #Topic of first chain
    chains <- lapply(as.list(1:chain), function(x){
      for(top in 1:K){
        corrTop[top] <- cor(theta_all %>% dplyr::filter(Chain == "Chain 1") %>% dplyr::filter(Topic == paste0("Topic_", j)) %>% dplyr::select(topic.dis),
                            theta_all %>% dplyr::filter(Chain == paste0("Chain ",x)) %>% dplyr::filter(Topic == paste0("Topic_", top)) %>% dplyr::select(topic.dis))
        
      }
      return(which(corrTop == max(corrTop)))
    })
    
    aligned[j, 2:chain] <- unlist(chains)[2:chain]
  }
  
  return(aligned)
}
