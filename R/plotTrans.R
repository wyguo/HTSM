plotTrans <- function(gr,
                      exon_color='darkblue',
                      intron_color='black',
                      exon_width=5,
                      intron_width=0.6,
                      num_ticks=10,
                      zoom=NULL){
  
  data2plot <- as.data.frame(gr)
  s <- data2plot$strand[1]
  chr <- data2plot$seqnames[1]
  gene <- data2plot$gene_id[1]
  
  trans <- sort(unique(data2plot$transcript_id))
  y <- seq_along(trans)
  names(y) <- trans
  
  gr_split <- split(gr,gr$transcript_id)
  gr_range <- as.data.frame(unlist(range(gr_split)))
  gr_range$transcript_id <- rownames(gr_range)
  
  ### exons
  data2plot$y <- y[data2plot$transcript_id]
  
  ### introns
  gr_range$y <- y[gr_range$transcript_id]
  
  ### add arrows
  data2arrow <- split(data2plot,data2plot$transcript_id)
  if(s=='+'){
    data2arrow <- lapply(data2arrow, function(x){
      x <- x[-1,]
      data.frame(start=x$start-1,end=x$start,y=x$y)
    })
  } else if (s=='-'){
    data2arrow <- lapply(data2arrow, function(x){
      x <- x[-nrow(x),]
      data.frame(start=x$end+1,end=x$end,y=x$y)
    })
  } else {
    data2arrow <- NULL
  }
  data2arrow <- do.call(rbind,data2arrow)
  
  data2y <- data.frame(label=data2plot$transcript_id,y=data2plot$y) 
  data2y <- data2y[!duplicated(data2y),]
  
  data2x <- seq(from=min(data2plot$start),to=max(data2plot$end),length.out = num_ticks)
  data2x <- as.integer(data2x)
  
  g <- ggplot(data2plot,aes(x=start,y=y,xend=end,yend=y))+
    geom_segment(data = gr_range,aes(x=start,y=y,xend=end,yend=y),size=intron_width)+
    geom_segment(data = data2arrow,aes(x=start,y=y,xend=end,yend=y),size=intron_width,
                 arrow = arrow(length = unit(0.5,"strwidth", "A")),
                 color=intron_color)+
    geom_segment(data = data2plot, size = exon_width,color = exon_color)+
    scale_y_continuous(breaks = data2y$y,labels = data2y$label)+
    scale_x_continuous(breaks = data2x,labels = data2x,position = "top",
                       expand = expansion(mult = c(0.02, 0.05)))+
    theme_bw()+
    labs(title=paste0(chr,': ',gene,' (strand:',s,')'))+ 
    coord_cartesian(ylim = c(min(data2y$y)-0.5,max(data2y$y)+0.5))+
    theme(axis.title = element_blank(),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line.x = element_line(colour = "black"),
          axis.ticks.y = element_blank())
  
  # ggplotly(g)
  if(!is.null(zoom))
    g <- g + facet_zoom(xlim = zoom,zoom.size = 1)
  g
}


plotRanges <- function(x, xlim=x, main=deparse(substitute(x)),
                       col="black", sep=0.5, ...)
{
  height <- 1
  if (is(xlim, "IntegerRanges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
  title(main)
  axis(1)
}

