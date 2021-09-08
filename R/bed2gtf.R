#' Convert gtf file to bed12 file. 
#' @param bed_file A bed file name to save.
#' @param gtf_file A gtf file. 
#' @return A gtf file will be saved to the directory. 
bed2gtf <- function(bed_file,gtf_file = NULL){
  if(is.null(gtf_file))
    gtf_file <- gsub('.bed','.gtf',bed_file)
  rtd_bed <- import(bed_file)
  rtd_gtf <- unlist(blocks(rtd_bed))
  rtd_gtf$gene_id <- gsub(';.*','',names(rtd_gtf))
  rtd_gtf$transcript_id <- gsub('.*;','',names(rtd_gtf))
  rtd_gtf$feature <- 'exon'
  rtd_gtf <- sort(rtd_gtf,by=~seqnames + start + end)
  export_gtf(gr = rtd_gtf,file2save = gtf_file)
}