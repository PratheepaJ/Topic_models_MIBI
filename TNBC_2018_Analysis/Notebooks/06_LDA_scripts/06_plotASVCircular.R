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

plotASVCircular <- function(
  ps,
  K,
  beta_aligned,
  varnames = c("iterations", "topic", "rsv_ix"),
  value_name = "beta_h",
  taxonomylevel = "Class",
  thresholdASVprop = 0.008
){

  topic <- rsv_ix <- beta_h <- rsv <- NULL
  Phylum <- Class <- Order <- Family <- Genus <- NULL
  beta_median <- median <- NULL


  beta_hat <- beta_aligned %>%
    reshape2::melt(
      varnames = varnames,
      value.name = value_name) %>%
    as_tibble()

  beta_hat$rsv <- rownames(
    tax_table(ps)
  )[beta_hat$rsv_ix]


  taxa <- phyloseq::tax_table(ps) %>%
    as.data.frame() %>%
    as_tibble()

  taxa$rsv <- rownames(
    tax_table(ps)
  )

  beta_hat <- beta_hat %>%
    dplyr::left_join(
      taxa,
      by = "rsv"
    ) %>%
    dplyr::mutate(
      topic = paste("Topic", topic)
    )


  beta_hat$Class <- factor(beta_hat$Class)
  beta_hat$rsv <- factor(beta_hat$rsv)
  beta_hat$rsv_ix <- factor(beta_hat$rsv_ix)
  beta_hat$topic <- factor(beta_hat$topic)


  beta_summary <- beta_hat %>%
    dplyr::group_by(
      rsv_ix,
      topic
    ) %>%
    dplyr::summarise(
      beta_median = median(beta_h),
      rsv = rsv[1],
      Phylum = Phylum[1],
      Class = Class[1],
      Order = Order[1],
      Family = Family[1],
      Genus = Genus[1]
    )


  beta_subset <- beta_summary

  beta_subset$rsv_ix <- rep(
    seq_len(
      nrow(beta_subset) / K
    ),
    each = K
  )


  beta_subset <- beta_subset %>%
    dplyr::arrange(
      rsv_ix,
      topic
    )

  beta_subset <- beta_subset %>%
    dplyr::mutate(
      Class = factor(
        Class,
        levels = unique(beta_subset$Class)
      ),
      Topic = stringr::str_remove(topic, "Topic ")
    )


  beta_subset$Class <- beta_subset$Class %>%
    as.character()

  beta_subset$Class <- ifelse(
    is.na(beta_subset$Class),
    "other",
    beta_subset$Class
  )

  beta_subset$Class <- factor(
    beta_subset$Class
  )


  beta_subset$Topic <- factor(
    beta_subset$Topic,
    levels = seq(1,K) %>%
      as.character()
  )

  ####

  beta_subset_wide <- dplyr::select(
    beta_subset,
    rsv_ix,
    rsv,
    topic,
    Class,
    beta_median
  ) %>%
    as.data.frame()

  beta_subset_wide <- dplyr::select(
    beta_subset_wide,
    -rsv_ix
  )

  beta_subset_wide <- tidyr::spread(
    beta_subset_wide,
    key = topic,
    value = beta_median
  )

  rownames(beta_subset_wide) <- beta_subset_wide$rsv

  beta_subset_wide$Class <- beta_subset_wide$Class %>%
    as.character()

  beta_subset_wide$Class <- ifelse(
    is.na(beta_subset_wide$Class),
    "other",
    beta_subset_wide$Class)

  beta_subset_wide$Class <- factor(
    beta_subset_wide$Class
  )

  # Choose rsv if at least in one topic the feature distribution is greater than the cutoff = 0.008
  choose_rsv_logical <- apply(
    beta_subset_wide[, 3:dim(beta_subset_wide)[2]],
    1,
    function(x){
      sum(x >= thresholdASVprop) >= 1
    }
  )

  choose_rsv <- beta_subset_wide$rsv[which(
    choose_rsv_logical
  )] %>%
    as.character()


  beta_subset_wide <- beta_subset_wide[choose_rsv, ]

  tree <- phyloseq::phy_tree(
    prune_taxa(
      choose_rsv,
      ps
    )
  )

  circ <- ggtree::ggtree(
    tree,
    layout = "circular"
  )

  beta_subset_wide <- beta_subset_wide[tree$tip.label, ]

  beta_subset_wide <- beta_subset_wide[ , c("rsv", "Class", paste0("Topic ", seq(1,K)))]

  color <- randomcoloR::distinctColorPalette(
    length(
      unique(
        beta_subset_wide$Class
      )
    )
  )

  color_topic <- grDevices::rainbow(K)


  df <- data.frame(
    Class = beta_subset_wide$Class
  )

  rownames(df) <- tree$tip.label


  df2 <- beta_subset_wide[, 3:dim(beta_subset_wide)[2]]

  rownames(df2) <- tree$tip.label



  p <- gheatmap(
    p = circ,
    data = df,
    offset = -.1,
    width = .1,
    colnames_angle = 95,
    colnames_offset_y = .5,
    font.size = 5) +
    ggplot2::scale_fill_manual(
      values = color,
      name = "Class"
    )


  for(i in 1:K){
    if(i == 1){
      p <- p +
        ggnewscale::new_scale_fill()
    }
    df2_i <- dplyr::select(
      beta_subset_wide, (i+2)
    )
    rownames(df2_i) <- tree$tip.label
    p <- gheatmap(
      p,
      df2_i,
      offset = i*.08,
      width = .1,
      colnames_angle = 90,
      colnames_offset_y = .25,
      font.size = 6,
      high = "dodgerblue",
      low = "gray98",
      legend_title = expression(beta[k]))
  }

  p <- p +
    theme_minimal(
      base_size = 20
    ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )

  return(p)
}
