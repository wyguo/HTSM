#' Filter collasped transcript datasets
#' @param data_dir Data directory
rtdFilter <- function(data_dir){

  message('|==>   Filter the transcript dataset (ca. 1min): ',Sys.time(),'   <==|')
  message('Step 1: Process SJ and TSS/TES filters')
  start.time <- Sys.time()

  load(file.path(data_dir,'HC_TSS.RData'))
  load(file.path(data_dir,'HC_TES.RData'))
  load(file.path(data_dir,'collasped_bed.RData'))
  load(file.path(data_dir,'sj_trans_filter.RData'))

  if(!file.exists(file.path(data_dir,'HC_TSS.RData')) |
     !file.exists(file.path(data_dir,'HC_TES.RData')) |
     !file.exists(file.path(data_dir,'collasped_bed.RData')) |
     !file.exists(file.path(data_dir,'sj_trans_filter.RData')))
    stop('SJ and TSS/TES analysis must be conducted before RTD filter')

  prefix <- list.files(path = data_dir,
                       pattern = '_collasped.bed$',
                       full.names = F,
                       recursive = T)
  prefix <- gsub('.bed','_filtered',prefix)

  ## Get the labels of HC TSS and TES
  TES_label <- unique(HC_TES$label)
  TSS_label <- unique(HC_TSS$label)

  ## Get transcript from rtd
  gr <- collasped_bed[collasped_bed$transcript_id %in% sj_trans_filter$trans_keep,]

  ## Match labels of TSS
  gr_start <- resize(gr,fix = 'start',width = 1,ignore.strand=F)
  gr_start_label <- paste0(seqnames(gr_start),';',strand(gr_start),';',start(gr_start))
  trans_start <- gr_start$transcript_id[gr_start_label %in% TSS_label]

  ## Match labels of TES
  gr_end <- resize(gr,fix = 'end',width = 1,ignore.strand=F)
  gr_end_label <- paste0(seqnames(gr_end),';',strand(gr_end),';',end(gr_end))
  trans_end <- gr_end$transcript_id[gr_end_label %in% TES_label]

  ## Subset rtd with HC transcripts
  trans <- unique(intersect(trans_start,trans_end))
  rtd_bed <- gr[gr$transcript_id %in% trans,]

  ## get exons of the transcriptome
  rtd_gtf <- unlist(blocks(rtd_bed))
  rtd_gtf$gene_id <- gsub(';.*','',names(rtd_gtf))
  rtd_gtf$transcript_id <- gsub('.*;','',names(rtd_gtf))
  rtd_gtf$feature <- 'exon'
  rtd_gtf <- sort(rtd_gtf,by=~seqnames + start + end)

  message('Step 2: Summarise the final RTD')
  s <- rtdSummary(rtd_gtf)
  write.csv(s,file=file.path(data_dir,paste0(prefix,'_summary.csv')))

  message('Step 3: Export the final RTD')
  # export_gtf(gr = rtd_gtf,file2save = file.path(data_dir,paste0(prefix,'.gtf')))
  export(object = rtd_bed,file.path(data_dir,paste0(prefix,'.bed')),format = 'bed')

  file2merge <- data.frame(file_name=file.path(normalizePath(data_dir,winslash = '/'),
                                               paste0(prefix,'.bed')),
                           cap_flag='capped',
                           merge_priority='1,1,1',
                           source_name=prefix)

  write.table(file2merge, file = file.path(data_dir,'file2merge.txt'),
              sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  ##########################################################################
  rm(collasped_bed)
  rm(gr)
  rm(gr_start)
  rm(gr_end)
  rm(rtd_bed)
  rm(rtd_gtf)
  gc()

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message('Done: ',Sys.time())
  message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units),'\n')
}
