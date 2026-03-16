### FISHER ----

fisher_plot_asym_pdf <- function(DF, main="", norm=T, range=c(-2,2), x_axis_title = "", y_axis_title = ""){
  if (is.null(norm)) norm <- F

  if (is.null(range)){
    if (norm) range <- c(-2,2)
    else range <- c(min(DF$lfc[is.finite(DF$lfc)]),max(DF$lfc[is.finite(DF$lfc)]))
  }
  if(norm)DF <- DF %>% mutate(lfc=rangeMinMaxAsym(x = lfc,nmin = -2,nmax=0,pmin = 0,pmax = 2))
  colfunc <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))
  p <- ggplot(DF, aes(x = COL, y = ROW, fill=lfc, group = pv, text = fm)) +      geom_tile(col=NA) +
    theme_minimal()+
    coord_equal(ratio=1) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top",expand=c(-1,0)) + ylab(y_axis_title) + xlab(x_axis_title) +
    theme(title = element_text(size=12,face = "bold"),
          axis.text.x=element_text(angle = 90,size=10, vjust=0, hjust=0.6,face = "bold",margin=margin(20,0,0,0)),
          axis.text.y=element_text(size=10, margin=margin(0,0,0,20),face = "bold"),
          axis.title.x =element_text(size=18, vjust = 20, margin=margin(0,200,0,0),face = "bold"),
          axis.title.y =element_text(size=18,vjust = 0.8, margin=margin(0,0,0,0),face = "bold"),
          plot.margin = unit(c(2,1,2,1), "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +

    geom_text(aes(label = paste0(sign, "\n",intersect)),size=4) +
    ggtitle(main) +
    scale_fill_gradientn(colours = colfunc(10) ,
                         guide = guide_colorbar(barwidth = 1,
                                                title = paste0(ifelse(norm, '_____\nScaled','_____'),'\nLog2 \nOddsRatio\n_____'),
                                                title.theme = element_text(size=8, vjust=0.9, face = "bold",angle = 0),
                                                barheight = 5,
                                                nbin = 40,
                                                ticks = FALSE),
                         breaks=c(-max(abs(range)),0,max(abs(range))),
                         labels=c(-max(abs(range)),0,max(abs(range))),
                         limits=c(-max(abs(range)),max(abs(range))))

  return(p)
}


### GRAPHICAL FUNCTIONS ----

ggPostTreatment <- function(gg_obj,p) {
  if(!is.null(p$xlim))
    gg_obj <- gg_obj + xlim(p$xlim)
  if(!is.null(p$ylim))
    gg_obj <- gg_obj + ylim(p$ylim)

  ## change scale log  ----
  if(!is.null(p$scalelog)){
    if(p$scalelog %in% c('xy', 'yx', 'y'))
      gg_obj <- gg_obj + scale_y_continuous(trans='log10')
    if(p$scalelog %in% c('xy', 'yx', 'x'))
      gg_obj <- gg_obj + scale_x_continuous(trans='log10')
  }

  # add title + labs  ----
  if(!is.null(p$lab_col))
    gg_obj <- gg_obj + labs(col = p$lab_col)
  if(!is.null(p$lab_fill))
    gg_obj <- gg_obj + labs(fill = p$lab_fill)
  if(!is.null(p$lab_x))
    gg_obj <- gg_obj + xlab(p$lab_x)
  if(!is.null(p$lab_y))
    gg_obj <- gg_obj + ylab(p$lab_y)
  if(!is.null(p$title))
    gg_obj <- gg_obj + ggtitle(p$title)

  # Remove all legends from plot
  if(!is.null(p$legend))
    if(p$legend == "F") gg_obj <- gg_obj + theme(legend.position = "none")

  # add colorscale
  if(!is.null(p$scale_col))
    if(!is.null(p$scale_col$breaks))
      gg_obj <- gg_obj + ggplot2::scale_color_manual(values = unlist(p$scale_col$values), breaks = unlist(p$scale_col$breaks))
  else
    gg_obj <- gg_obj + ggplot2::scale_color_manual(values = unlist(p$scale_col$values))


  if(!is.null(p$scale_fill))
    if(!is.null(p$scale_fill$breaks))
      gg_obj <- gg_obj + ggplot2::scale_fill_manual(values = unlist(p$scale_fill$values), breaks = unlist(p$scale_fill$breaks))
  else
    gg_obj <- gg_obj + ggplot2::scale_fill_manual(values = unlist(p$scale_fill$values))

  # gg_obj <- gg_obj + theme(plot.title = element_text(hjust = 0.5))
  # gg_obj <- gg_obj + theme(plot.subtitle = element_text(hjust = 0.5))
  gg_obj + theme_classic()
  gg_obj
}


### OTHERS ----

rangeMinMax <- function(x,min=-1,max=1)
{
  return(min+((x- min(x,na.rm=T))*(max-min)) /(max(x,na.rm=T)-min(x,na.rm=T)))
}

rangeMinMaxAsym <- function(x,nmin=-1,nmax=0,pmin=0,pmax=1)
{
  neg <- x[x<=0]
  pos <- x[x>=0]
  if(length(neg)){
    if(length(unique(neg))-1)
      ret_neg <- nmin+((neg- min(c(-0.0001,neg),na.rm=T))*(nmax-nmin)) /(max(c(-0.0001,neg),na.rm=T)-min(c(-0.0001,neg),na.rm=T))
    else
      ret_neg <- rep(nmin,length(neg))
  }else{
    ret_neg <- NULL
  }
  if(length(pos)){
    if(length(unique(pos))-1)
      ret_pos <- pmin+((pos- min(c(0.0001,pos),na.rm=T))*(pmax-pmin)) /(max(c(0.0001,pos),na.rm=T)-min(c(0.0001,pos),na.rm=T))
    else
      ret_pos <- rep(pmax,length(pos))
  }else{
    ret_pos <- NULL
  }
  idx_neg <- x<0
  idx_pos <- x>=0
  x[idx_neg] <- ret_neg
  x[idx_pos] <- ret_pos
  x
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}


redim_matrix <- function(
    mat,
    target_height = 100,
    target_width = 100,
    summary_func = function(x) mean(x, na.rm = TRUE),
    output_type = 0.0, #vapply style
    n_core = 1 # parallel processing
) {

  if(target_height > nrow(mat) | target_width > ncol(mat)) {
    stop("Input matrix must be bigger than target width and height.")
  }

  seq_height <- round(seq(1, nrow(mat), length.out = target_height + 1))
  seq_width  <- round(seq(1, ncol(mat), length.out = target_width  + 1))

  # complicated way to write a double for loop
  do.call(rbind, parallel::mclapply(seq_len(target_height), function(i) { # i is row
    vapply(seq_len(target_width), function(j) { # j is column
      summary_func(
        mat[
          seq(seq_height[i], seq_height[i + 1]),
          seq(seq_width[j] , seq_width[j + 1] )
        ]
      )
    }, output_type)
  }, mc.cores = n_core))
}
