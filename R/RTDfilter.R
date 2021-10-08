#' Filter collasped transcript datasets
#' @param data_dir Data directory
#' @return Results are saved in \code{data_di}
#' 
rtdFilter <- function(data_dir){
  message('|==========> Apply SJ and TSS/TES filters: ',Sys.time(),' <==========|')
  message('Step 1: Load intermediate datasets\n')
  start.time <- Sys.time()
  
  if(!file.exists(file.path(data_dir,'HC_TSS.RData')) |
     !file.exists(file.path(data_dir,'HC_TES.RData')) |
     !file.exists(file.path(data_dir,'trans_bed2filter.RData')) |
     !file.exists(file.path(data_dir,'sj_filter.RData')) |
     !file.exists(file.path(data_dir,'input_files.RData'))
     )
    stop('SJ and TSS/TES analysis must be conducted before RTD filter')

  load(file.path(data_dir,'HC_TSS.RData'))
  load(file.path(data_dir,'HC_TES.RData'))
  load(file.path(data_dir,'trans_bed2filter.RData'))
  load(file.path(data_dir,'sj_filter.RData'))
  load(file.path(data_dir,'input_files.RData'))
  
  message('Step 2: Apply splice junction filter\n')
  ##############################################################
  ###---> filter transcripts with error splice junction
  sj_filter <- unique(unlist(sj_filter))
  
  exons <- unlist(blocks(trans_bed2filter))
  ## add gene and transcript ids
  exons$gene_id <- gsub(';.*','',names(exons))
  exons$transcript_id <- gsub('.*;','',names(exons))
  names(exons) <- NULL
  introns <- exon2intron(exons,sorted = F)
  trans2filter <- unique(introns$transcript_id[introns$label %in% sj_filter])
  
  message('Step 3: Apply TSS and TES filter\n')
  ##############################################################
  ###---> filter transcripts with low quality TSS and TES
  TES_label <- unique(HC_TES$label)
  TSS_label <- unique(HC_TSS$label)
  gr <- trans_bed2filter[!(trans_bed2filter$transcript_id %in% trans2filter),]
  
  ## Match labels of TSS
  gr_start <- resize(gr,fix = 'start',width = 1,ignore.strand=F)
  gr_start_label <- paste0(seqnames(gr_start),';',strand(gr_start),';',start(gr_start))
  trans_start <- gr_start$transcript_id[gr_start_label %in% TSS_label]

  ## Match labels of TES
  gr_end <- resize(gr,fix = 'end',width = 1,ignore.strand=F)
  gr_end_label <- paste0(seqnames(gr_end),';',strand(gr_end),';',end(gr_end))
  trans_end <- gr_end$transcript_id[gr_end_label %in% TES_label]

  ## Subset rtd with HC transcripts
  trans <- intersect(trans_start,trans_end)
  rtd_bed <- gr[gr$transcript_id %in% trans,]
  
  message('Step 4: Save the filtered RTD\n')
  file2save <- input_files$trans_merged[1]
  
  prefix <- paste0(gsub('.bed','',basename(file2save)),'_filtered')
  export(object = rtd_bed,file.path(data_dir,paste0(prefix,'.bed')),format = 'bed')

  file2merge <- data.frame(file_name=file.path(normalizePath(data_dir,winslash = '/'),
                                               paste0(prefix,'.bed')),
                           cap_flag='capped',
                           merge_priority='1,1,1',
                           source_name=prefix)

  write.table(file2merge, file = file.path(data_dir,'tama_merge_final_rtd_file.txt'),
              sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  ##########################################################################

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message('Done: ',Sys.time())
  message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units),'\n')
}
