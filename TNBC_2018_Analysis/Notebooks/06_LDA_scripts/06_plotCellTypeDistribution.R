#' Plot the ASVs distribution in each topic.
#'
#' @param beta_aligned List. returned object of betaAligned().
#' @param varnames Character vector. Use to name column of long formatted data frame of ASVs distribution. An example is c("iterations", "topic", "rsv_ix")
#' @param value_name Character. Use to name the summary of posterior samples of ASVs proportions.
#' @param taxonomylevel Character. Use to label ASVs.
#' @param thresholdASVprop Numberic. Use to select ASVs for circular plot.
#' @inheritParams alignmentMatrix
#' @return A ggplot2 object. A circular plot annotated with phylogenetic tree.
#' @importFrom reshape2 melt
#' @import phyloseq
#' @importFrom dplyr left_join mutate arrange summarise group_by
#' @importFrom tibble as_tibble
#' @importFrom stringr str_remove
#' @importFrom ggtree ggtree gheatmap
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom grDevices rainbow
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggnewscale new_scale_fill
#' @importFrom stats median
#' @importFrom tidyr spread
#' @export
#'
#'
#'
#'
#'


plotCellTypeDistribution <- function(
  spe,
  K,
  beta_aligned,
  col_names_beta_hat = c("iterations", "Topic", "Cell.Type", "beta_h")
  #varnames = c("iterations", "Topic", "Cell.Type"),
  #value_name = "beta_h"
) {
  
  Topic <- Cell.Type <- beta_h <- NULL
  beta_median <- median <- NULL
  
  
  dimnames(beta_aligned)[[2]] <- c(paste0("Topic_", seq(1,K)))
  dimnames(beta_aligned)[[3]] <- unique(spe$mm)
  
  
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
  
  
  p <- ggplot2::ggplot(beta_summary,
                       ggplot2::aes(x = topic,
                                    y = Cell.Type,
                                    fill = beta_median)) +
    ggplot2::geom_tile() +
    ggplot2::ylab("Cell Phenotype") +
    ggplot2::xlab("Topic") +
    ggplot2::scale_fill_gradientn(name = "Median Cell \nPhenotype \ndistribution",
                         colors = c("gray98", "dodgerblue")) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}








