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
    warm_up_iter = NULL,
    iter = 2000,
    chain = 2,
    col_names_theta_all = c("iteration", "Sample", "Topic", "topic.dis"),
    col_names_beta_hat = c("iterations", "Topic", "Cell.Type", "beta_h"),
    stan_file_path = here::here("Notebooks", "06_LDA_scripts", "06_lda.stan"),
    cellTypeIndexToPlot = c(1:16),
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
  
  # determine the iteration used in posterior sampling (subtract warm up iterations)
  if (is.null(warm_up_iter)) {
    iterUse = iter / 2
  } else {
    iterUse = iter - warm_up_iter
  } # number of iterations used in MCMC(subtract warm-up samples) 
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

  
  theta_long <- reshape2::melt(theta_aligned) |> as_tibble()
  colnames(theta_long) = col_names_theta_all
  #theta_long$Chain = paste0("Chain ", rep(seq(1, chain), each = (iterUse)))
  
  
  # factor variables
  #theta_long$Chain <- factor(theta_long$Chain)
  theta_long$Topic <- factor(theta_long$Topic)
  theta_long$Sample <- factor(theta_long$Sample)
  #theta_all$pna <- factor(theta_all$pna)
  
  # summary theta distribution
  theta_summary <- ( theta_long 
    |> dplyr::group_by(
      Sample,
      Topic
      #,pna
    ) 
    |> dplyr::summarize(
      median.topic.dis = median(topic.dis)
    ) 
    |> dplyr::ungroup()
    |> dplyr::mutate(Topic = stringr::str_remove(Topic, "Topic_"))
    |> dplyr::mutate(
      Topic = factor(
        Topic,
        #levels = rev(stringr::str_c("Topic_",1:K)
        levels = seq(1, K))
        )
      )
  
  
  
  # plot
  p_topic_prop <- ggplot2::ggplot(
    theta_summary,
    ggplot2::aes(
      x = Sample,
      y = Topic)
  )
  
  
  p_topic_prop <- p_topic_prop +
    ggplot2::geom_tile(
      ggplot2::aes(fill = median.topic.dis)
    ) +
    ggplot2::ylab("Community") +
    # ggplot2::facet_grid(
    #   .~Sample,
    #   scale = "free"
    # ) +
    ggplot2::xlab("Sample") +
    ggplot2::ylab("Community") +
    ggplot2::scale_fill_gradientn(
      name = "Median Community \ndistribution",
      colors = c("gray98", "#E69F00")
      ) +
    # ggplot2::scale_alpha(
    #   name = "Median Topic \ndistribution"
    # ) +
    ggplot2::theme_classic(
      base_size = 20
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      strip.text.x = ggplot2::element_text(angle = 90)
    ) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))
    
  # source(here::here("Notebooks", 
  #                   "06_LDA_scripts", "06_LDA_analysis.R"))
  # 
  # 
  # p_topic_prop <- plotTopicProportion(
  #   spe = spe,
  #   theta_aligned = theta_aligned,
  #   K = K,
  #   chain = 4,
  #   iter = 2000,
  #   SampleID_name = "sample_id"
  # )
  
  
  
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
    ggplot2::ylab("Cell Phenotype") +
    ggplot2::xlab("Community") +
    ggplot2::scale_fill_gradientn(name = "Median Cell Phenotype \ndistribution",
                                  colours = c("gray98", "dodgerblue")) +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5))
  
  
  
  
  
  ################### Model Diagnostics ###################
  
  Rhat_theta  <- ESS_bulk_theta <- matrix(
    nrow = dim(theta_aligned)[2],
    ncol = dim(theta_aligned)[3]
  )
  
  # ESS_bulk_theta <- matrix(
  #   nrow = dim(theta_aligned)[2],
  #   ncol = dim(theta_aligned)[3]
  # )
  
  
  for(sam in 1:dim(theta_aligned)[2]){
    for(top in 1:dim(theta_aligned)[3]){
      sims_theta <- matrix(
        theta_aligned[ ,sam , top],
        nrow = (iterUse),
        ncol = chain,
        byrow = FALSE
      )
      Rhat_theta[sam, top] <- rstan::Rhat(sims_theta)
      ESS_bulk_theta[sam, top] <- rstan::ess_bulk(sims_theta)
    }
    
  }
  
  Rhat_theta <- as.vector(Rhat_theta)
  ESS_bulk_theta <- as.vector(ESS_bulk_theta)
  
  
  Rhat_beta <- ESS_bulk_beta <- matrix(
    nrow = dim(beta_aligned)[2],
    ncol = dim(beta_aligned)[3]
  )
  # ESS_bulk_beta <- matrix(
  #   nrow = dim(beta_aligned)[2],
  #   ncol = dim(beta_aligned)[3]
  # )
  
  
  for(top in 1:dim(beta_aligned)[2]){
    for(fea in 1:dim(beta_aligned)[3]){
      sims_beta <- matrix(
        beta_aligned[ , top, fea],
        nrow = (iterUse),
        ncol = chain,
        byrow = FALSE)
      Rhat_beta[top, fea] <- rstan::Rhat(sims_beta)
      ESS_bulk_beta[top, fea] <- rstan::ess_bulk(sims_beta)
      
    }
    
  }
  
  Rhat_beta <- as.vector(Rhat_beta)
  ESS_bulk_beta <- as.vector(ESS_bulk_beta)
  
  
  Rhat <- c(Rhat_theta, Rhat_beta)
  
  
  ESS_bulk <- c(ESS_bulk_theta, ESS_bulk_beta)
  
  
  
  # R hat ~ 1.05
  p_rhat <- ggplot2::ggplot(
    data.frame(Rhat = Rhat)
  ) +
    ggplot2::geom_histogram(
      aes(x = Rhat),
      fill = "lavender",
      colour = "black",
      bins = 100
    ) +
    ggplot2::theme(
      plot.title = element_text(hjust = 0.5)
    )  +
    ggplot2::theme_minimal(base_size = 20) +
    ggplot2::xlab("")
  
  
  
  # ESS bulk and ESS tail at least 100 per Markov Chain in order to be reliable and indicate that estimates of respective posterior quantiles are reliable
  
  p_ess_bulk <- ggplot2::ggplot(
    data.frame(ESS_bulk = ESS_bulk)
  ) +
    ggplot2::geom_histogram(
      aes(x = ESS_bulk),
      fill = "lavender",
      colour = "black",
      bins = 100
    ) +
    ggplot2::theme(
      plot.title = element_text(hjust = 0.5)
    )   +
    ggplot2::theme_minimal(
      base_size = 20
    ) +
    ggplot2::xlab("")
  
  
  
  
  ################### Model Assessment ###################
  
  value <- Var2 <- NULL
  
  x <- dtm
  
  # draws from posterior predictive distribution
  x_sim <- samples$x_sim[1:iterUse, , ] # iteration * samples * topics
  
  # choose only the first chain
  max_all <- apply(x_sim[1, , ], 2, max)
  
  # find the maximum of each cell types in each iteration acorss sample 
  for (i in 2:iterUse){
    max_x_sim_i <- apply(
      x_sim[i, ,], 
      2,
      max)
    max_all <- rbind(max_all, max_x_sim_i)
  }
  
  rownames(max_all) <- c(
    paste0(
      "x_max_rep",
      seq(1, iterUse)
    )
  )
  
  colnames(max_all) <- colnames(x)
  #return(max_all)
  
  # subset the interested phenotype
  max_all <- max_all[, cellTypeIndexToPlot]
  max_all_long <- reshape2::melt(max_all)
  
  
  # finding the observed maximum value 
  x_max_cellCount <- data.frame(
    Var1 = rep(
      "x_max_obs",
      dim(x)[2]
    ),
    Var2 = colnames(x),
    count = apply(x, 2, max)
  )
  
  # subset the interested phenotype
  x_max_cellCount <- x_max_cellCount[cellTypeIndexToPlot, ]
  
  
  # plotting
  p_hist <- ggplot2::ggplot(
    data = max_all_long
  ) +
    ggplot2::geom_histogram(
      aes(
        x = value,
        group = Var2
      ),
      color = "#0072B2",
      fill = "#0072B2",
      bins = 50) +
    ggplot2::xlab("maximum")+
    
    ggplot2::facet_wrap(~Var2, nrow = 4, scales = "free_x") +
    
    ggplot2::geom_vline(
      data = x_max_cellCount,
      aes(xintercept = count),
      color = "#CC79A7"
    ) +
    ggplot2::theme_update(
      text = element_text(size = 8)
    )
  
  
  
  
  
  ################### Return List  ###################
  if (cond == 3) {
    return_list <- list(
      DTM = dtm,
      stan.fit = stan.fit,
      pos_samples = samples,
      theta = theta,
      beta = beta,
      alignment_matrix = aligned,
      theta_aligned = theta_aligned,
      beta_aligned = beta_aligned,
      plot_topic_prop = p_topic_prop,
      plot_cellType_prop = p_cellType_prop,
      plot_ESS = p_ess_bulk,
      plot_rhat = p_rhat,
      plot_model_assessment = p_hist
      
    )
    
  } else if (cond == 0) {
    return_list <- list(
      DTM = dtm,
      stan.fit = stan.fit,
      pos_samples = samples,
      theta = theta,
      beta = beta,
      plot_topic_prop = p_topic_prop,
      plot_cellType_prop = p_cellType_prop,
      plot_ESS = p_ess_bulk,
      plot_rhat = p_rhat,
      plt_model_assessment = p_hist
    )
  }
  
  proc_time_all <- proc.time() - t0   # processing time
  print(proc_time_all)
  
  return(return_list)
}