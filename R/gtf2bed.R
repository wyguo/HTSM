#' Convert gtf file to bed12 file. 
#' @param gtf_file A gtf file. 
#' @param bed_file A bed file name to save.
#' @return A bed file will be saved to the directory. 
gtf2bed <- function(gtf_file,bed_file=NULL){
  if(is.null(bed_file))
    bed_file <- gsub('.gtf','.bed',gtf_file)
  message('Read gtf file ...')
  gr <- import(gtf_file)
  type <- ifelse('feature' %in% colnames(mcols(gr)),'feature','type')
  if(!(type %in% c('feature','type')))
    stop('The object does not have a "feature" or "type" column of exon labels')
  
  message('Convert to bed file ...')
  exon <- gr[mcols(gr)[,type]=='exon',]
  mapping <- data.frame(gene_id = exon$gene_id,transcript_id = exon$transcript_id)
  mapping <- mapping[!duplicated(mapping),]
  rownames(mapping) <- mapping$transcript_id
  
  x <- split(exon,exon$transcript_id)
  # bed <- rtracklayer::asBED(exon)
  x_range <- range(x)
  if (any(elementNROWS(x_range) != 1L))
    stop("Empty or multi-strand/seqname elements not supported by BED")
  bed <- unlist(x_range, use.names=FALSE)
  mcols(bed) <- mcols(x)
  mcols(bed)$name <- names(x)
  x_ranges <- ranges(unlist(x, use.names=FALSE))
  message('ord_start')
  ord_start <- order(IRanges::togroup(IRanges::PartitioningByEnd(x)), start(x_ranges))
  message('ord_start_finish')
  x_ranges <- shift(x_ranges, 1L - rep(start(bed), elementNROWS(x)))[ord_start]
  # mcols(bed)$blocks <- relist(x_ranges, x)
  blocks <- relist(x_ranges, x)
  
  
  
  # blocks <- bed$blocks
  # bed$blocks <- NULL
  strand_info <- as.character(strand(bed))
  strand_info[strand_info == "*"] <- NA
  
  df <- data.frame(seqnames = seqnames(bed), 
                   start = start(bed) - 1, 
                   end = end(bed), 
                   names = paste0(mapping[bed$name,'gene_id'],';',bed$name), 
                   score = 0, 
                   strand = strand_info, 
                   thickStart = start(bed) - 1, 
                   thickEnd = end(bed), 
                   itemRgb = 1, 
                   blockCount = elementNROWS(blocks), 
                   blockSizes = unlist(lapply(width(blocks), paste, collapse = ","), use.names = FALSE),
                   blockStarts = unlist(lapply(start(blocks) - 1, paste, collapse = ","), use.names = FALSE))
  rownames(df) <- NULL
  
  message('Export bed file ...')
  write.table(df, file = bed_file, sep = "\t", col.names = FALSE,
              row.names = FALSE, quote = FALSE, na = ".")
  
}