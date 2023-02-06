#' Perform LDA topic modelling to SpatialExperiment object and return posterior sampling results and visulization
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 


LDA_analysis <- function(
    spe,
    SampleID_name = "sample_id",
    cellType_name = "mm",
    K,
    alpha,
    gamma,
    iter,
    chain,
    col_names_theta_all = c("iteration", "Sample", "Topic", "topic.dis"),
    col_names_beta_hat = c("iterations", "Topic", "Cell.Type", "beta_h"),
    stan_file_path = here::here("Notebooks", "06_LDA_scripts", "06_lda.stan"),
    load_file = TRUE,
    save_file = FALSE,
    file_default = TRUE,
    file_fold = NULL
) {
  
  
  t0 <- proc.time()
  
  
  
  ################### Doucment-Term matrix ###################
  # Create Doucment-Term matrix
  # return dtm
  
  #dtm = table(spe$SampleID_name, spe$cellType_name)
  dtm = table(
    colData(spe)[SampleID_name],
    colData(spe)[cellType_name]
  )
  
  # Dimension names of the DTM(sample identity and cell type name)
  dim_names <- dimnames(dtm)
  
  
  
  ################### Posterior Sampling ###################
  # Perform LDA and save results or load LDA results
  # return samples, a list
  
  
  # number of iteration and chains
  iter <- iter
  chain <- chain
  
  # file name
  fileN <- paste0("JABES_Endo_all_K_",
                  K,
                  "_ite_",
                  iter,
                  ".RData")
  
  
  # determine file path to save or load LDA posterior rsults
  if (file_default == TRUE) {
    file_path = file = here::here("Output", "Data",
                                  paste0(file_fold), paste0(fileN))
  } else {
    file_path = paste0(file_fold, fileN) # file_fold must end with "/" if file_default = FALSE
  } # if file path is default

  
  
  if (load_file == FALSE) {
    # Compute posterior sampling
    x <- dtm
    print(dim(x))
    
    # set dimensions names of dtm to NULL
    dimnames(dtm) <- NULL
    
    # stan input data
    ## theta[d] ~ dirichlet(alpha), alpha pseudo-count for each topic
    ## beta[k] ~ dirichlet(gamma), gamma pseudo-count for each cell type in each topic
    stan.data <- list(
      K = K,
      V = ncol(x),
      D = nrow(x),
      n = x,
      alpha = rep(alpha, K),
      gamma = rep(gamma, ncol(x))
    )
    
    
    # apply dataset
    t1 <- proc.time()
    stan.fit <- rstan::stan(
      file = stan_file_path,
      data = stan.data,
      # must be list object
      iter = iter,
      chains = chain,
      sample_file = NULL,
      diagnostic_file = NULL,
      cores = 4,
      control = list(adapt_delta = 0.9),
      save_dso = TRUE,
      algorithm = "NUTS"
    )
    proc_time <- proc.time() - t1   # processing time
    

    # save file with specified name to specified path
    if (save_file == TRUE) {
      save(stan.fit,
           file = file_path)
      }
  
  } else {

    load(file = file_path)
      
    } # if compute LDA
  
  # extract posterior results
  samples <- rstan::extract(
    stan.fit,
    permuted = TRUE,
    inc_warmup = FALSE,
    include = TRUE
  )
  
  
  
  
  
  ################### Theta and Beta Estimations ###################
  # extract theta and beta estimations
  # return theta and beta
  
  # theta estimation
  theta <- samples$theta
  # set dimension names
  dimnames(theta)[[2]] <- dim_names[[1]]
  # dimnames(theta[2]) <- dim_names[1]
  dimnames(theta)[[3]] <- c(paste0("Topic_", seq(1,K)))
  
  
  # beta estimation
  beta <- samples$beta
  # set fimension names
  dimnames(beta)[[2]] <- c(paste0("Topic_", seq(1,K)))
  dimnames(beta)[[3]] <- dim_names[[2]]
  
  
  
  ################### Alignment matrix ###################
  # Align topics across chains (must apply when chains > 1)
  # return aligned, alignment matrix
  
  # check whether number of chains are greater than 1
  chain_cond <- isTRUE(chain > 1)
  
  iterUse = iter/2 # number of iterations used in MCMC(subtract warm-up samples) 
  cond = 0 # identifier for future if statements
  
  
  # array to a dataframe
  ## melt() takes wide-format data and melts into long-format data
  theta_all = reshape2::melt(theta)
  
  colnames(theta_all) = col_names_theta_all
  #colnames(theta_all) = c("iteration", "Sample", "Topic", "topic.dis")
  
  # Perfrom alignment
  if (chain_cond == TRUE) {
    
    # theta is a 3 dimensional array (iterations * samples * topic)
    Chain <- Topic <- topic.dis <- NULL
    
    
    # label corresponding chain number
    theta_all$Chain = paste0("Chain ", rep(seq(1, chain), each = (iterUse)))
    
    # factor variables
    # theta_all$Topic = factor(theta_all$Topic)
    # theta_all$Chain = factor(theta_all$Chain)
    
    # convert to character for left_join
    theta_all$Sample = as.character(theta_all$Sample)
    
    # join the SpatialExperimet object column data to theta_all
    ## colData gets the metadata (clinical data)
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
          corrTop[top] <- cor(
            theta_all |> dplyr::filter(Chain == "Chain 1") |> dplyr::filter(Topic == paste0("Topic_", j)) |> dplyr::select(topic.dis),
            theta_all |> dplyr::filter(Chain == paste0("Chain ",x)) |> dplyr::filter(Topic == paste0("Topic_", top)) |> dplyr::select(topic.dis))
          
        }
        return(which(corrTop == max(corrTop)))
      })
      
      aligned[j, 2:chain] <- unlist(chains)[2:chain]
    }
    
    cond = 1
  } # if chains condition is TRUE
  
  
  
  
  
  ################### Align Topic Proportion  ###################
  # align topics between chains
  # return theta_aligned
  
  if (cond == 1) {
    
    # switch samples and Topic dimension in array
    theta <- aperm(theta, c(1,3,2))
    
    theta_chain <- list()
    theta_chain[[1]] <- theta[1:(iterUse), ,]
    theta_chain[[1]] <- theta_chain[[1]][, aligned[,1],]
    
    theta_aligned <- theta_chain[[1]]
    
    for(ch in 2:chain){
      theta_chain[[ch]] <- theta[((ch-1)*(iterUse)+1):(ch*(iterUse)),,]
      theta_chain[[ch]] <- theta_chain[[ch]][,aligned[,ch],]
      
      theta_aligned <- abind::abind(theta_aligned, theta_chain[[ch]], along = 1)
    } # for loop
    
    # switch back samples and Topic dimension in array
    theta_aligned <- aperm(theta_aligned, c(1,3,2))
    
    cond = 2
  } # if cond is 1
  
  
  
  
  ################### Align Cell Type Proportion  ###################
  # align topics between chains
  # return beta_aligned
  
  if (cond == 2) {
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
    } # for loop
    
    cond = 3
  } # if cond is 2
  
  

  ################### Plot Topic Proportion  ###################
  # plot heatmap of topic distribution
  
  Sample <- Topic <- topic.dis <- NULL
  median.topic.dis <- median <- NULL
  
  # factor variables
  theta_all$Topic <- factor(theta_all$Topic)
  theta_all$Sample <- factor(theta_all$Sample)
  #theta_all$pna <- factor(theta_all$pna)
  
  # summary theta distribution
  theta_summary <- ( theta_all |> 
    dplyr::group_by(
      Sample,
      Topic
      #,pna
    ) 
    |> dplyr::summarize(
      median.topic.dis = median(topic.dis)
    ) 
    |> dplyr::ungroup()
    |> dplyr::mutate(
      Topic = factor(
        Topic,
        levels = rev(stringr::str_c("Topic_",1:K)
        )
      )
    )
  )
  
  
  # plot
  p_topic_prop <- ggplot2::ggplot(
    theta_summary,
    ggplot2::aes(
      x= Sample,
      y = Topic)
  )
  
  
  p_topic_prop <- p_topic_prop +
    ggplot2::geom_tile(
      ggplot2::aes(alpha = median.topic.dis)
    ) +
    # ggplot2::facet_grid(
    #   .~Sample,
    #   scale = "free"
    # ) +
    ggplot2::xlab("Sample") +
    # ggplot2::scale_fill_manual(
    #   name = "pna",
    #   values = c("steelblue1","green3")
    # ) +
    ggplot2::scale_alpha(
      name = "median topic \ndistribution"
    ) +
    ggplot2::theme_classic(
      base_size = 20
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      strip.text.x = ggplot2::element_text(angle = 90)
    ) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))
    
  
  
  
  
  ################### Plot Cell Type Proportion  ###################
  # plot heatmap of cell type proportion in each topics
  
  
  Topic <- Cell.Type <- beta_h <- NULL
  beta_median <- median <- NULL
  
  
  beta_hat <- reshape2::melt(beta_aligned) |> as_tibble()
  colnames(beta_hat) <- col_names_beta_hat
  
  
  # factor variables
  beta_hat$Cell.Type <- factor(beta_hat$Cell.Type)
  beta_hat$Topic <- factor(beta_hat$Topic)
  
  
  # compute median of beta
  beta_summary <- (
    beta_hat
    |> dplyr::group_by(Cell.Type, Topic)
    |> dplyr::summarise(
      Cell.Type = Cell.Type[1],
      beta_median = median(beta_h)
    )
  )
  
  beta_summary <- (
    beta_summary
    |> dplyr::mutate(topic = stringr::str_remove(Topic, "Topic_"))
    |> dplyr::mutate(topic = factor(topic, levels = seq(1, K)))
  )
  
  
  p_cellType_prop <- ggplot2::ggplot(beta_summary,
                                     ggplot2::aes(x = topic,
                                                  y = Cell.Type,
                                                  fill = beta_median)) +
    ggplot2::geom_tile() +
    ggplot2::ylab("Cell Type") +
    ggplot2::xlab("Topic") +
    ggplot2::scale_fill_gradientn(name = "Median Cell Type \ndistribution",
                                  colours = c("gray98", "dodgerblue")) +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5))
  
  
  
  
  ################### Return List  ###################
  if (cond == 3) {
    return_list <- list(
      dtm,
      samples,
      theta,
      beta,
      aligned,
      theta_aligned,
      beta_aligned,
      p_topic_prop,
      p_cellType_prop
    )
    
  } else if (cond == 0) {
    return_list <- list(
      dtm,
      samples,
      theta,
      beta,
      p_topic_prop,
      p_cellType_prop
    )
  }
  
  proc_time_all <- proc.time() - t0   # processing time
  print(proc_time_all)
  
  return(return_list)
}