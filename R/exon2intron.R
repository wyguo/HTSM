#' Get introns from exon coordinates
#' @param exons A GRanges object of exons
#' @param tx2gene A data frame of transcript-gene association
#' @param add_label Logical, whether to add label of seqnames:start-end;strand
#' @return A GRanges object of introns
#'
exon2intron <- function(exons,tx2gene=NULL,add_label=T,sorted=T){
  if(is.null(tx2gene)){
    tx2gene <- DataFrame(transcript_id=exons$transcript_id,gene_id=exons$gene_id)
    tx2gene <- tx2gene[!duplicated(tx2gene),]
    rownames(tx2gene) <- tx2gene$transcript_id
  }
  # message('Split exons')
  exons <- GenomicRanges::split(exons,exons$transcript_id)
  # message('Get introns')
  introns <- psetdiff(unlist(range(exons),use.names=FALSE),exons)
  introns <- unlist(introns)
  if(NROW(introns)==0)
    return('No introns')

  introns$transcript_id <- names(introns)
  introns$gene_id <- tx2gene[introns$transcript_id,'gene_id']
  introns$feature <- 'intron'
  names(introns) <- NULL
  introns$num_introns <- sequence(rle(introns$transcript_id)$length)

  if(sorted)
    introns <- sort(introns,by=~seqnames + start + end)
  if(add_label)
    introns$label <- paste0(seqnames(introns),':',start(introns),'-',end(introns),';',strand(introns))
  return(introns)
}

