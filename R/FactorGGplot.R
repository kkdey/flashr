#' Absolute Factor Loadings (in FLASH) Struture plot with ggplot2 package
#'
#' Make the traditional Structure histogram plot of absolute loadings from FLASH 
#' output
#'
#' @param loadings Absolute loadings for each sample. Usually a
#' sample by factors matrix in the FLASH or any factor analysis model output. 
#' @param annotation A data.frame of two columns: sample_id and label.
#' sample_id is the unique identifying number of each sample (alphanumeric).
#' label is a factor of labels, with levels of the factor
#' arranged in the order of the groups in the Structure (left to right).
#' @param palette A vector of colors assigned to the clusters. First color in
#' the vector is assigned to the cluster labeled as 1, and second color in the
#' vector is assigned to the cluster labeled as 2, etc. The number of colors
#' must be the same or greater than the number of clusters. The clusters not
#' assigned a color are filled with white in the figure. In addition, the
#' recommended choice of color palette here is RColorBrewer, for instance
#' RColorBrewer::brewer.pal(8, "Accent") or RColorBrewwer::brewer.pal(9, "Set1").
#' @param figure_title Title of the plot.
#' @param yaxis_label Axis label for the samples.
#' @param order_sample if TRUE, we order samples in each annotation batch
#' sorted by membership of most representative cluster. If FALSE, we keep
#' the order in the data.
#' @param sample_order_decreasing if order_sample TRUE, then this input
#' determines if the ordering due to main cluster is in ascending or descending
#' order.
#' @param split_line Control parameters for line splitting the batches in the
#' plot.
#' @param axis_tick Control parameters for x-axis and y-axis tick sizes.
#' @param plot_labels A logical parameter, if TRUE the function plots the axis labels.
#'
#' @return Plots the Structure plot visualization of the absolute values of loadings of 
#' FLASH model
#'
#' @import ggplot2
#' @importFrom cowplot ggdraw panel_border switch_axis_position plot_grid
#' @import plyr
#' @import reshape2
#' @export

FactorGGplot <- function(loadings, annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            figure_title = "",
                            yaxis_label = "Factor Types",
                            order_sample = TRUE,
                            sample_order_decreasing = TRUE,
                            split_line = list(split_lwd = 0.2,
                                              split_col = "black"),
                            plot_labels = TRUE,
                            legend_labels = "",
                            scale=TRUE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 3,
                                             axis_label_face = "bold") ) {
  
  if(scale){
    loadings <- apply(loadings,2,function(x) 
    {
      if(sd(x)!=0) {return (x/sd(x))}
      else {return (x)};
    })
  }
  
  # check if the number of colors is same as or more than the number of clusters
  if (dim(loadings)[2] > length(palette)) {
    stop("Color choices is smaller than the number of clusters!")
  }
  
  # check if rownames of loadings are unique
  if(length(unique(rownames(loadings))) != NROW(loadings)) {
    stop("loadings rownames are not unique!")
  }
  
  # check the annotation data.frame
  if (!is.data.frame(annotation))
    stop("annotation must be a data.frame")
  if (!all.equal(colnames(annotation), c("sample_id", "label")) ) {
    stop("annotation data.frame column names must be sample_id and label")
  }
  if ( length(unique(annotation$sample_id)) != NROW(loadings)) {
    stop("sample_id is not unique")
  }
  
  df_ord <- do.call(rbind,
                    lapply(1:nlevels(annotation$label), function(ii) {
                      temp_label <- levels(annotation$label)[ii]
                      temp_df <- loadings[which(annotation$label == temp_label), ]
                      
                      is_single_sample <-
                        ( length(temp_df) == nlevels(annotation$label)|
                            is.null(dim(temp_df)) )
                      # find the dominant cluster in each sample
                      if ( is_single_sample ) {
                        each_sample_order <- which.max(temp_df)
                      } else {
                        each_sample_order <- apply(temp_df, 1, which.max)
                      }
                      
                      # find the dominant cluster across samples
                      sample_order <- as.numeric(attr(table(each_sample_order), "name")[1])
                      
                      if (order_sample == TRUE & !is_single_sample) {
                        # reorder the matrix
                        temp_df_ord <- temp_df[order(temp_df[ , sample_order],
                                                     decreasing = sample_order_decreasing), ]
                      } else {
                        temp_df_ord <- temp_df
                      }
                      temp_df_ord
                    }) )
  
  df_mlt <- reshape2::melt(t(df_ord))
  df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
                                             "Var2" = "document"))
  df_mlt$document <- factor(df_mlt$document)
  df_mlt$topic <- factor(df_mlt$topic)
  
  # set blank background
  ggplot2::theme_set(ggplot2::theme_bw(base_size = 12)) +
    ggplot2::theme_update( panel.grid.minor.x = ggplot2::element_blank(),
                           panel.grid.minor.y = ggplot2::element_blank(),
                           panel.grid.major.x = ggplot2::element_blank(),
                           panel.grid.major.y = ggplot2::element_blank() )
  
  # inflat nubmers to avoid rounding errors
  value_ifl <- 1
  # number of ticks for the weight axis, including 0 and 1
  ticks_number <- 6
  
  # set axis tick positions
  label_count <- table(droplevels(annotation$label))
  label_count_cumsum <- cumsum(table(droplevels(annotation$label)))
  
  label_names <- levels(droplevels(annotation$label))
  label_breaks <- sapply(1:length(label_count), function(i) {
    if (i == 1) {
      if (label_count[i] == 1) bk <- 1
      if (label_count[i] > 1)  bk <- (label_count_cumsum[i] - 0)/2
      return(bk)
    }
    if (i > 1) {
      if (label_count[i] == 1) bk_interval <- 1
      if (label_count[i] > 1 ) {
        bk_interval <- (label_count_cumsum[i] - label_count_cumsum[i-1])/2 }
      bk <- label_count_cumsum[i-1] + bk_interval
      return(bk)
    }
  })
  names(label_breaks) <- label_names
  
  df_mlt_1 <- subset(df_mlt, value > 0);
  df_mlt_2 <- subset(df_mlt, value <=0);
  
  # make ggplot
  a <- ggplot2::ggplot() +
      ggplot2::xlab(yaxis_label) + ggplot2::ylab("") +
      ggplot2::scale_fill_manual(values = palette,
                                 labels = paste(1:dim(loadings)[2], legend_labels)) +
      ggplot2::theme(legend.position = "right",
                   legend.key.size = ggplot2::unit(.2, "cm"),
                   legend.text = ggplot2::element_text(size = 5),
                   ##<-- TBD: center legend title
                   #              legend.title = element_text(hjust = 1),
                   axis.text = ggplot2::element_text(size = axis_tick$axis_label_size,
                                                     face = axis_tick$axis_label_face),
                   axis.ticks.y = ggplot2::element_line(size = axis_tick$axis_ticks_lwd_y),
                   axis.ticks.x = ggplot2::element_line(size = axis_tick$axis_ticks_lwd_x),
                   axis.ticks.length = ggplot2::unit(axis_tick$axis_ticks_length, "cm"),
                   title = ggplot2::element_text(size = 6) ) +
    ggplot2::ggtitle(figure_title) +
    ggplot2::scale_x_discrete(breaks = as.character((levels(df_mlt$document)[round(label_breaks)])),
                              labels = names(label_breaks)) +
    # Add legend title
    ggplot2::labs(fill = "Factors") +
    ggplot2::coord_flip() +
    geom_bar(data = df_mlt_1, ggplot2::aes(x = document, y = value*1, fill = factor(topic)), stat = "identity") +
    geom_bar(data = df_mlt_2, ggplot2::aes(x = document, y = value*1, fill = factor(topic)), stat = "identity")
    cowplot::panel_border(remove = TRUE)
    
  
  # width = 1: increase bar width and in turn remove space
  # between bars
 
  # Add demarcation
  b <- a + ggplot2::geom_vline(
    xintercept = cumsum(table(droplevels(annotation$label)))[
      -length(table(droplevels(annotation$label)))] + .1,
    col = split_line$split_col,
    size = split_line$split_lwd)
  
  if (!plot_labels) {
    b
  } else {
    suppressWarnings(b <- cowplot::ggdraw(cowplot::switch_axis_position((b), axis = "y")))
    b
  }
}
