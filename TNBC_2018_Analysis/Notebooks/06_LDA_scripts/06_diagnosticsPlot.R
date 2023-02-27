#' Diagnostics plots - ESS and Rhat.
#'
#' @param theta_aligned array. Three dimensional array returned by thetaAligned().
#' @param beta_aligned array. Three dimensional array returned by betaAligned().
#' @inheritParams alignmentMatrix
#'
#' @return list. ggplots of effective sample size and rhat plots.
#' @importFrom rstan Rhat ess_bulk
#' @importFrom ggplot2 theme_minimal theme geom_histogram ggplot aes element_text xlab
#' @export
#'
#'


diagnosticsPlot <- function(
  theta_aligned,
  beta_aligned,
  warm_up_iter = NULL,
  iter = 2000,
  chain = 4
){
  
  # determine the iteration used in posterior sampling (subtract warm up iterations)
  if (is.null(warm_up_iter)) {
    iterUse = iter / 2
  } else {
    iterUse = iter - warm_up_iter
  }
  
  Rhat_theta <- matrix(
    nrow = dim(theta_aligned)[2],
    ncol = dim(theta_aligned)[3]
  )

  ESS_bulk_theta <- matrix(
    nrow = dim(theta_aligned)[2],
    ncol = dim(theta_aligned)[3]
  )


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


  Rhat_beta <- matrix(
    nrow = dim(beta_aligned)[2],
    ncol = dim(beta_aligned)[3]
  )
  ESS_bulk_beta <- matrix(
    nrow = dim(beta_aligned)[2],
    ncol = dim(beta_aligned)[3]
  )


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

  return(list(p_ess_bulk, p_rhat))
}

