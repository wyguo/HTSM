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


intron2SJchain <- function(intron){
  gr <- GenomicRanges::split(intron,intron$transcript_id)
  s <- strand(gr)
  s <- unlist(runValue(s))
  
  chr <- seqnames(gr)
  chr <- unlist(runValue(chr))
  
  sj_chain <- paste0(start(intron),'-',end(intron))
  sj_chain <- GenomicRanges::split(sj_chain,intron$transcript_id)
  sj_chain <- sapply(sj_chain, function(x) paste0(x,collapse = ';'))
  chr <- chr[names(sj_chain)]
  s <- s[names(sj_chain)]
  label <- paste0(chr,':',sj_chain,':',s)
  names(label) <- names(sj_chain)
  label
}

getMultiExonTrans <- function(exon){
  gr <- GenomicRanges::split(exon,exon$transcript_id)
  n <- elementNROWS(gr)
  trans <- names(n)[n>1]
  exon[exon$transcript_id %in% trans,]
}

getMonoExonTrans <- function(exon){
  gr <- GenomicRanges::split(exon,exon$transcript_id)
  n <- elementNROWS(gr)
  trans <- names(n)[n==1]
  exon[exon$transcript_id %in% trans,]
}

getMonoExonGenes <- function(exon){
  gr <- GenomicRanges::split(exon,exon$gene_id)
  n <- elementNROWS(gr)
  genes <- names(n)[n==1]
  exon[exon$gene_id %in% genes,]
}

#getGeneLoci->getGeneRange
getGeneRange <- function(exon){
  gr <- GenomicRanges::split(exon,exon$gene_id)
  gr <- unlist(range(gr))
  gr$gene_id <- names(gr)
  names(gr) <- NULL
  # gr <- sort(gr,by=~seqnames+start+end)
  gr
}

#getTransCallapse->getGeneLoci
getGeneCollapse <- function(exon){
  gr <- GenomicRanges::split(exon,exon$gene_id)
  gr <- reduce(gr)
  gr <- unlist(gr)
  gr$gene_id <- names(gr)
  names(gr) <- NULL
  # gr <- sort(gr,by=~seqnames+start+end)
  gr
}

getOverlapRange <- function(exon){
  gr <- GenomicRanges::split(exon,exon$gene_id)
  gr <- range(gr)
  gr <- unlist(gr)
  gr <- reduce(gr)
  gr$gene_id <- paste0(seqnames(gr),':',start(gr),';',strand(gr))
  gr
}


#getTransLoci->getTransRange
getTransRange <- function(exon){
  gr <- GenomicRanges::split(exon,exon$transcript_id)
  gr <- unlist(range(gr))
  gr$transcript_id <- names(gr)
  names(gr) <- NULL
  # gr <- sort(gr,by=~seqnames+start+end)
  gr
}



cleanMcols <- function(gr){
  meta <- mcols(gr)
  type_id <- ifelse('type' %in% colnames(meta),'type','feature')
  meta <- meta[,c('source',type_id,'score','phase','transcript_id','gene_id')]
  colnames(meta)[2] <- 'type'
  mcols(gr) <- meta
  gr
}

subGR <- function(gr,cood){
  cood <- gsub(',','',cood)
  s <- unlist(strsplit(cood,split = ':|-'))
  gr[seqnames(gr)==s[1] & start(gr) > s[2] & end(gr) < s[3]]
}