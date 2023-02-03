#' aligned topics across chains
#'
#' @param theta array. three dimensional array (iterations * samples * topic).
#' @param ps phyloseq object.
#' @param K integer. number of topics.
#' @param iterUse integer. number of iterations used in each chain after subtracting warm-up samples.
#' @param chain integer. number of chains used in MCMC.
#' @param SampleID_name character. variable name that contains sample IDs.
#'
#' @return matrix. A matrix of topic assignment for each chain.
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join filter select
#' @import phyloseq
#' @importFrom reshape2 melt
#' @importFrom stats cor
#' @export
#'

alignmentMatrix <- function(theta,
                            ps,
                            K,
                            iterUse = 1000,
                            chain = 4,
                            SampleID_name = "SampleID"){
  Chain <- Topic <- topic.dis <- NULL
  # theta is a 3 dimensional array (iterations * samples * topic)
  dimnames(theta)[[2]] <- phyloseq::sample_names(ps) %>%
    as.character()
  dimnames(theta)[[3]] <- c(paste0("Topic_", seq(1,K)))

  # array to a dataframe
  theta_all = reshape2::melt(theta)
  colnames(theta_all) = c("iteration", "Sample", "Topic", "topic.dis")
  theta_all$Chain = paste0("Chain ", rep(seq(1, chain), each = (iterUse)))
  theta_all$Sample = factor(theta_all$Sample)
  theta_all$Topic = factor(theta_all$Topic)

  # join the phyloseq sample data to theta_all
  sam = sample_data(ps) %>%
    data.frame()
  theta_all$Sample = as.character(theta_all$Sample)
  theta_all = dplyr::left_join(
    theta_all,
    sam,
    by =c("Sample"= SampleID_name)
    )
  theta_all$Chain = factor(theta_all$Chain)


  aligned <- matrix(nrow = K, ncol = chain)
  aligned[, 1] <- seq(1, K)
  corrTop <- numeric()
  for(j in 1:K){#Topic of first chain
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
