#' Factor Loadings (in FLASH) Multi-panel Struture plot with ggplot2 package
#'
#' Make the multi-panel bar chart plot of loadings from FLASH 
#' output
#'
#' @param loadings loadings for each sample. Usually a
#' sample by factors matrix in the FLASH or any factor analysis model output. 
#' @param annotation A data.frame of two columns: sample_id and label.
#' sample_id is the unique identifying number of each sample (alphanumeric).
#' label is a factor of labels, with levels of the factor
#' arranged in the order of the labels in the Structure (left to right).
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

FactorGGBar <- function(loadings, annotation,
                         palette = list("mid"="white", 
                                        "low"="red", 
                                        "high"="blue", 
                                        "midpoint"=0),
                         figure_title = "",
                         yaxis_label = "Label type",
                         split_line = list(split_lwd = 0.2,
                                           split_col = "white"),
                         axis_tick = list(axis_ticks_length = .1,
                                          axis_ticks_lwd_y = .1,
                                          axis_ticks_lwd_x = .1,
                                          axis_label_size = 3,
                                          axis_label_face = "bold"),
                        legend_labels = NULL,
                        scale=TRUE,
                        panel=list(panel_rows=2,
                                   panel_title="",
                                   panel_title_fontsize=10,
                                   panel_title_font=3)) {
  
  if(scale){
    loadings <- apply(loadings,2,function(x) 
                                  {
                                      if(sd(x)!=0) {return (x/sd(x))}
                                      else {return (x)}
    })
  }
  
  if(is.null(legend_labels)){
    legend_labels <- rep("", dim(loadings)[2]);
  }
  
  # check if the number of colors is same as or more than the number of clusters
  
  # check if rownames of loadings are unique
  if(length(unique(rownames(loadings))) != NROW(loadings)) {
    stop("loadings rownames are not unique!")
  }
  
  ## check if legend labels size matches with loadings
  
  if(length(legend_labels) != dim(loadings)[2]){
    stop("number of loadings do not match with number of legend labels")
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
  
  graphList <- vector(mode="list");
  ncols <- dim(loadings)[2]
  
  for(n in 1:ncols){
    df_ord <- loadings[,n];
    df_mlt <- reshape2::melt(t(df_ord))
    df_mlt <- plyr::rename(df_mlt, replace = c("Var1" = "topic",
                                               "Var2"="document"))
    
    df_mlt$document <- factor(df_mlt$document)
    df_mlt$topic <- factor(df_mlt$topic)
    
    suppressMessages(graphList[[n]] <- ggplot2::ggplot(df_mlt,
                                      ggplot2::aes(x = document,
                                                   y = value,
                                                   fill =value)) + ggplot2::xlab(yaxis_label) + ggplot2::ylab("") + ggplot2::scale_fill_manual(values = palette) +
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
      ggplot2::ggtitle(paste0("Factor: ", n, legend_labels[n]))  +
      # Add label axis labels
      ggplot2::scale_x_discrete(breaks = as.character((levels(df_mlt$document)[round(label_breaks)])),
                                labels = names(label_breaks))  + geom_bar(stat = "identity",position = "stack",width = 1) + #make the bars
      coord_flip() + #flip the axes so the test names can be horizontal  
      #define the fill color gradient: blue=positive, red=negative
      scale_fill_gradient2(name = "Loading", 
                           high = palette$high, mid = palette$mid, low = palette$low, 
                           midpoint=palette$midpoint, guide=F) + 
      ggplot2::geom_vline(xintercept = cumsum(table(droplevels(annotation$label)))[
        -length(table(droplevels(annotation$label)))] + .1, col = split_line$split_col, size = split_line$split_lwd))
  }
  
  panel$panel_cols <- ceiling(length(graphList)/panel$panel_rows)
  library(gridExtra)
  library(grid)
  suppressWarnings(do.call("grid.arrange", 
          args = list(grobs=graphList,
                      ncol = panel$panel_cols,
                      nrow = panel$panel_rows,
                      top=textGrob(paste0(panel$panel_title),
                                   gp=gpar(fontsize=panel$panel_title_fontsize,
                                           font=panel$panel_title_font)))));
         
}

