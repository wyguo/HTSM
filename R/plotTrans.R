plotTrans <- function(gr,
                      exon_color='dodgerblue3',
                      intron_color='black',
                      exon_width=5,
                      intron_width=0.6,
                      x_ticks=10,
                      arrow_ticks=10,
                      zoom=NULL,
                      xlim=NULL,
                      vline=F,
                      vline_color='red',
                      vline_width=0.3){
  if(NROW(gr)==0){
    g <- ggplot() + 
      annotate("text", x = 4, y = 25, label = 'Empty') + 
      theme_void()
    return(g)
  }
  
  type <- ifelse('type' %in% colnames(mcols(gr)),'type','feature')
  gr <- gr[mcols(gr)[,type]=='exon',]
  
  gr$transcript_id <- paste0(gr$gene_id,': ',gr$transcript_id)
  
  gr.list <- split(gr,gr$transcript_id)
  gr.range <- range(gr.list)
  
  data2line <- unlist(gr.range)
  data2line$transcript_id <- names(gr.list)
  data2line <- as.data.frame(data2line)
  
  trans <- sort(unique(data2line$transcript_id))
  trans.idx <- seq_along(trans)
  names(trans.idx) <- trans
  data2line$y <- trans.idx[data2line$transcript_id]
  
  exon2plot <- as.data.frame(gr)
  exon2plot$y <- trans.idx[exon2plot$transcript_id]
  
  loci <- range(gr,ignore.strand=T)
  arrow.idx <- IRanges(start=seq(start(loci),end(loci),length.out = arrow_ticks),width = 2)
  
  idx <- findOverlaps(arrow.idx,gr@ranges)
  # data2arrow1 <- as.data.frame(arrow.idx[idx@from])
  data2arrow1 <- as.data.frame(arrow.idx[idx@from])
  data2arrow1$transcript_id <- gr$transcript_id[idx@to]
  data2arrow1$y <- trans.idx[data2arrow1$transcript_id]
  data2arrow1$strand <- as.character(strand(gr))[idx@to]
  data2arrow1$color <- 'white'
  data2arrow1$start[data2arrow1$strand=='-'] <- data2arrow1$start[data2arrow1$strand=='-']+2
  
  introns <- exon2intron(gr)
  if(!(class(introns) == "GRanges")){
    data2arrow <- data2arrow1
  } else {
    idx <- findOverlaps(arrow.idx,introns@ranges)
    if(NROW(idx)==0){
      data2arrow <- data2arrow1
    } else {
      data2arrow2 <- as.data.frame(arrow.idx[idx@from])
      data2arrow2$transcript_id <- introns$transcript_id[idx@to]
      data2arrow2$y <- trans.idx[data2arrow2$transcript_id]
      data2arrow2$strand <- as.character(strand(introns))[idx@to]
      data2arrow2$color <- 'black'
      data2arrow2$start[data2arrow2$strand=='-'] <- data2arrow2$start[data2arrow2$strand=='-']+2
      data2arrow <- rbind(data2arrow1,data2arrow2)
    }
  }
  
  if(is.null(xlim))
    xlim <- c(min(data2line$start),max(data2line$end))
  
  data2x <- seq(from=xlim[1],to=xlim[2],length.out = x_ticks)
  data2x <- as.integer(data2x)
  
  if(is.null(xlim))
    xlim <- c(min(data2x),max(data2x))
  
  chr <- unique(as.character(seqnames(gr)))
  chr <- paste0(chr,':',paste0(xlim,collapse = '-'))
  
  gene <- unique(gr$gene_id)
  gene <- paste0(gene,collapse = ' | ')
  
  if(!is.null(zoom))
    xlim <- NULL
  
  g <- ggplot(data = data2line,aes(x=start,xend=end,y=y,yend=y))+
    geom_segment(data = data2line,aes(x=start,xend=end,y=y,yend=y),
                 size=intron_width,color = intron_color)+
    geom_segment(data = exon2plot,aes(x=start,xend=end,y=y,yend=y),
                 size=exon_width,color = exon_color)+
    geom_segment(data = data2arrow,aes(color=color),
                 arrow = arrow(length = unit(0.5,"strwidth", "A")))+
    scale_color_manual(values = c('white'='white','black' = intron_color))+
    scale_y_continuous(breaks = data2line$y,labels = data2line$transcript_id)+
    scale_x_continuous(breaks = data2x,labels = data2x,position = "top",
                       expand = expansion(mult = c(0.02, 0.05)))+
    theme_bw()+
    labs(title=chr)+ 
    coord_cartesian(xlim = xlim,ylim = c(min(data2line$y)-0.5,max(data2line$y)+0.5),
                    expand = T)+
    theme(axis.title = element_blank(),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line.x = element_line(colour = "black"),
          axis.ticks.y = element_blank(),
          legend.position = 'none')
    # xlim(xlim)
  
  if(vline){
    vl <- unique(c(start(gr),end(gr)))
    g <- g + geom_vline(xintercept = vl,linetype='dashed',size=vline_width,color=vline_color)
  }
  
  
  # ggplotly(g)
  # if(!is.null(xlim))
  #   g <- g +  coord_cartesian(xlim = xlim)
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

