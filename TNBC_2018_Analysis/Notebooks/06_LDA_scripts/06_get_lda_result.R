#' Perform LDA and save results
#' 
#' @param K integer. Number of topics.
#' @param alpha vector. Hyperparameters for document-topic Dirichlet distribution.
#' @param gamma vector. Hyperparameters for topic-word Dirichlet distribution.
#' @param dtm matrix. Document-Term matrix for the dataset.
#' @param iter integer. Number of iterations used in each chain.
#' @param chains integer. Number of chains used in MCMC.
#' @param file_fold string. File folder to save the LDA result, default set as NULL.
#' @param save_file boolean. Option to save file or not, default set as FALSE.
#' 
#' @return list. Large list contains parameters theta and beta estimations.
#' @importFrom here here
#' @importFrom rstan stan extract
#' 


get_lda_result <- function(K, 
                           alpha, gamma, 
                           dtm, 
                           iter, 
                           chains,
                           file_fold = NULL,
                           save_file = FALSE) {
  
  # setting seed will align topics
  # not applying seed, topics will differ but similar probabilities
  x <- dtm
  # View(x)
  dimnames(x) <- NULL
  
  # theta[d] ~ dirichlet(alpha), alpha pseudo-count for each topic
  # beta[k] ~ dirichlet(gamma), gamma pseudo-count for each ASV in each topic
  
  # n is DTM
  stan.data <- list(
    K = K,
    V = ncol(x),
    D = nrow(x),
    n = x,
    alpha = rep(alpha, K),
    gamma = rep(gamma, ncol(x))
  )
  
  #number of iteration
  iter <- iter
  # file name
  fileN <- paste0("JABES_Endo_all_K_",
                  K,
                  "_ite_",
                  iter,
                  ".RData")
  
  # number of chains
  chains <- chains
  
  # apply dataset
  t1 <- proc.time()
  stan.fit <- rstan::stan(
    file = here::here("Notebooks", "06_LDA_scripts", "06_lda.stan"),
    data = stan.data,
    # must be list object
    iter = iter,
    chains = chains,
    sample_file = NULL,
    diagnostic_file = NULL,
    cores = 4,
    control = list(adapt_delta = 0.9),
    save_dso = TRUE,
    algorithm = "NUTS"
  )
  proc_time <- proc.time() - t1   # processing time
  
  
  # save file with specified name
  if (save_file == TRUE) {
    save(stan.fit, 
         file = here::here("Output", "Data", 
                           paste0(file_fold), paste0(fileN)))
  }
  
  
  res <- rstan::extract(
    stan.fit,
    permuted = TRUE,
    inc_warmup = FALSE,
    include = TRUE
  )
  
  return(list(stan.fit,
              res))
}