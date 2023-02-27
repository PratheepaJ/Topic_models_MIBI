#' Model assessment for LDA based on the maximum abundance for each taxon.
#'
#' @param stan.fit An instance of stanfit.
#' @param ASVsIndexToPlot An integer vector. Use to select ASVs to show the goodness of fit in histograms.
#' @inheritParams alignmentMatrix
#' @return A ggplot2 object. Histogram of data generated from the posterior estimates and the observed data.
#' @importFrom rstan stan
#' @importFrom ggplot2 facet_wrap geom_vline theme_update aes
#' @import phyloseq
#' @importFrom reshape2 melt
#' @export
#'
modelAssessment <- function(
  #spe,
  dtm,
  #stan.fit = NULL,
  samples,
  warm_up_iter = NULL,
  iter = 2000,
  cellTypeIndexToPlot = c(1:16)
){

  value <- Var2 <- NULL
  
  # determine the iteration used in posterior sampling (subtract warm up iterations)
  if (is.null(warm_up_iter)) {
    iterUse = iter / 2
  } else {
    iterUse = iter - warm_up_iter
  }
  
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
    
    ggplot2::facet_wrap(~Var2, nrow = 4) +
    
    ggplot2::geom_vline(
      data = x_max_cellCount,
      aes(xintercept = count),
      color = "#CC79A7"
    ) +
    ggplot2::theme_update(
      text = element_text(size = 8)
    )
  return(p_hist)

}
