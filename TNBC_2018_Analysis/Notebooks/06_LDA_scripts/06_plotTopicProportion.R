#' Plot the topic distribution in all specimens.
#'
#' @param theta_aligned List. returned object of thetaAligned().
#' @param col_names_theta_all Character vector. Column names for long formatted theta. An example is c("iteration", "Sample", "Topic", "topic.dis")
#' @param SampleIDVarInphyloseq Character. SampleID variable name in phyloseq object.
#' @param design Factors of interests.
#' @inheritParams alignmentMatrix
#'
#' @return A ggplot2 object. Histrogram of median of topic proportion in each sample.
#' @import phyloseq
#' @import ggplot2
#' @importFrom dplyr mutate summarize group_by left_join ungroup
#' @importFrom stringr str_c
#' @importFrom stats median
#' @export
#'

plotTopicProportion <- function(
  spe,
  theta_aligned,
  K,
  col_names_theta_all = c("iteration", "Sample", "Topic", "topic.dis"),
  chain = 4,
  iterUse = 1000,
  SampleID_name = "sample_id"
  #design = ~pna
){
  
  Sample <- Topic <- pna <- topic.dis <- NULL
  median.topic.dis <- median <- NULL



  dimnames(theta_aligned)[[2]] <- unique(spe$sample_id)
  dimnames(theta_aligned)[[3]] <- c(paste0("Topic_", seq(1,K)))


  theta_all <- reshape2::melt(theta_aligned)

  colnames(theta_all) <- col_names_theta_all

  theta_all$Chain <- paste0(
    "Chain ",
    rep(seq(1, chain),
        each = iterUse
    )
  )

  sam = (colData(spe)[1:9]
         |> unique()
         |> data.frame()
  )

  theta_all$Sample <- theta_all$Sample |>
    as.character()

  theta_all <- dplyr::left_join(
    theta_all,
    sam,
    by =c("Sample" = SampleID_name)
  )

  theta_all$Chain <- factor(theta_all$Chain)
  theta_all$Topic <- factor(theta_all$Topic)
  theta_all$Sample <- factor(theta_all$Sample)
  #theta_all$pna <- factor(theta_all$pna)


  theta_summary <- theta_all |>
    dplyr::group_by(
      Sample,
      Topic
      #pna
    ) |>
    dplyr::summarize(
      median.topic.dis = median(topic.dis)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      Topic = factor(
        Topic,
        levels = rev(stringr::str_c("Topic_",1:K)
        )
      )
    )

  p <- ggplot2::ggplot(
    theta_summary,
    ggplot2::aes(x= Sample,
        y = Topic)
  )
  
  p <- p+
    ggplot2::geom_tile(
      ggplot2::aes(alpha = median.topic.dis)
    ) +
    # ggplot2::facet_grid(
    #   .~Sample,
    #   scale = "free"
    # ) +
    ggplot2::xlab("Sample") +
    ggplot2::scale_fill_manual(
      name = "Median Topic \ndistribution",
      values = c("gray98", "dodgerblue")
    ) +
    # ggplot2::scale_fill_gradientn(
    #   name = "Median Topic \ndistribution",
    #   colours = c("#999999", "#0072B2")
    # ) +
    ggplot2::scale_alpha(
      name = "Median Topic \ndistribution"
    ) +
    ggplot2::theme_classic(
      base_size = 20
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      strip.text.x = ggplot2::element_text(angle = 90)
    ) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))
  p
}
