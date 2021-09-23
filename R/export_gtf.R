#' Export Granges object to gtf file
#' @param gr Grange object 

export_gtf <- function(gr,file2save,
                       source.col='source',
                       feature.col=ifelse('type' %in% colnames(mcols(gr)),'type','feature'),
                       score='.',
                       phase='.',
                       mandatory = c("gene_id", "transcript_id"),
                       return.obj=F,
                       write.obj=T,
                       do_sort=T){
  options(scipen=100)
  if(NROW(gr) == 0){
    write.table(NULL, file=file2save, sep="\t", quote=F, row.names=F, col.names=F)
    message('The gtf is empty')
  } else {
    if(do_sort)
      gr <- sort(gr, by = ~ seqnames + start + end)
    ###check factor
    check <- mcols(gr)
    idx <- sapply(check@listData, is.factor)
    idx <- names(idx)[idx]
    for(i in idx){
      mcols(gr)[,i] <- as.character(mcols(gr)[,i])
    }
    
    idx <- unique(c('gene_id','transcript_id',feature.col,source.col,mandatory))
    if(!(source.col %in% colnames(mcols(gr))))
      mcols(gr)[,source.col] <- 'source'
    mcols(gr) <- mcols(gr)[,idx]
    ## change strands
    strd <- as.character(strand(gr))
    strd[strd=='*'] <- '.'
    score[is.na(score)] <- '.'
    phase[is.na(phase)] <- '.'
    
    gtf <- data.frame(chr=as.character(seqnames(gr)),
                      source=as.character(mcols(gr)[,source.col]),
                      feature=as.character(mcols(gr)[,feature.col]),
                      start=as.numeric(start(gr)),
                      end=as.numeric(end(gr)),
                      score=score,
                      strand=strd,
                      phase=phase,
                      attribute=makeGtfAttributes(df = as.data.frame(mcols(gr)),mandatory = mandatory),
                      stringsAsFactors = F)
    if(write.obj)
      write.table(gtf, file=file2save, sep="\t", quote=F, row.names=F, col.names=F)
    if(return.obj)
      return(gtf)
  }
  
  
}

gr2gtf <- function(gr,file2save,source='source',feature='feature',score='.',phase='.',return.obj=F,
                   mandatory = c("gene_id", "transcript_id")){
  ## change strands
  strd <- as.character(strand(gr))
  strd[strd=='*'] <- '.'
  gtf <- data.frame(chr=as.character(seqnames(gr)),
                    source=source,
                    feature=as.character(mcols(gr)[,feature]),
                    start=as.numeric(start(gr)),
                    end=as.numeric(end(gr)),
                    score=score,
                    strand=strd,
                    phase=phase,
                    attribute=makeGtfAttributes(df = as.data.frame(mcols(gr)),mandatory = mandatory),
                    stringsAsFactors = F)
  gtf
}


#' Convert metadata information of a GRanges object to the attribute of a ftf
#' @param df the metadata information of a GRanges object
#'
makeGtfAttributes <- function(df, cols=NULL,mandatory = c("gene_id", "transcript_id")) {
  if (is.null(cols))
    cols = colnames(df)
  # make sure that gene_id and transcript_id are the first two columns
  # mandatory = c("gene_id", "transcript_id")
  o = match(c(mandatory, setdiff(cols, mandatory)), cols)
  if (any(is.na(o[1:length(mandatory)]))) {
    warning("mandatory gtf attributes gene_id or transcript_id missing")
    o = o[!is.na(o)]
  }
  cols = cols[o]
  
  lab <- apply(df[,cols], 1, function(x){
    x <- as.character(x)
    idx <- which(!is.na(x))
    paste0(cols[idx],' "',x[idx],'"',collapse = '; ')
  })
  
  # lab <- lapply(cols, function(idx){
  #   paste0(idx,' "',df[,idx],'"',sep = '')
  # })
  # 
  # lab <- Reduce(function(...) paste(..., sep = "; "), lab)
  lab
}

gff2gtf <- function(gff){
  strd <- as.character(strand(gff))
  strd[strd=='*'] <- '.'
  gtf <- data.frame(chr=as.character(seqnames(gr)),
                    source=source,
                    feature=as.character(mcols(gr)[,feature]),
                    start=as.numeric(start(gr)),
                    end=as.numeric(end(gr)),
                    score=score,
                    strand=strd,
                    phase=phase,
                    attribute=makeGtfAttributes(df = as.data.frame(mcols(gr))),
                    stringsAsFactors = F)
}

