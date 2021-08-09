#' Transcript start site (TSS) and end site (TES) analysis
#' @param data_dir The data directory of TAMA output.
#' @param cut_off A numeric value of BH adjsuted p-value (fdr), p-value or probablity cut-off.
#' @param cut_at The \code{cut_off} is applied to "fdr" (adjustd p-value), "pval" (p-value) or "prob" (probability).
#' @param TSS_region,TES_region An integer of enriched TSS/TES wobbling region. TSSs/TESs locate within enriched 
#' site+/-site_region are all treated as high confident sites. Default: 50.
#' @param min_read The minimum reads in the window around TSS/TES to support high confident sites. Default: 2.
#' @param bin An integer of window size to aggregate the reads upstream and downstream a TSS/TES.
#' The total reads in [site-bin,site+bin] will be calculated. Default: 5.

TSanalysis <- function(data_dir,cut_off = 0.05,cut_at = 'fdr',
                       TSS_region = 50,TES_region = 50,min_read = 2,bin = 5){
  
  start.time <- Sys.time()
  ###########################################################################
  ### Step 1: Filter transcripts with polyA
  ###########################################################################
  ##---> get polya transcript
  
  message('|==>       TSS and TES analysis (ca. 1.5min): ',Sys.time(),'      <==|')
  message('Step 1: Filter transcripts with polyA')
  file2read <- list.files(path = data_dir,
                          pattern = '_polya.txt$',
                          full.names = T,
                          recursive = T)
  
  if(!file.exists(file2read))
    stop('File "*_polya.txt" is missing in the data directory:\n',data_dir)
  
  polya <- read.delim(file2read)
  cluster_polya <- unique(polya$cluster_id)
  
  ###########################################################################
  ### Step 2: Generate statistics of TSS and TES
  ###########################################################################
  message('Step 2: Generate statistics of TSS and TES')
  file2read <- list.files(path = data_dir,
                          pattern = '_collasped_trans_read.bed$',
                          full.names = T,
                          recursive = T)
  if(!file.exists(file2read))
    stop('File "*_collasped_trans_read.bed" is missing in the data directory:\n',data_dir)
  
  ## input the trans models
  trans_bed <- import(file2read)
  trans_bed$transcript_id <- gsub(';.*','',trans_bed$name)
  trans_bed$gene_id <- gsub('[.].*','',trans_bed$transcript_id)
  trans_bed$cluster_id <- gsub('.*;','',trans_bed$name)
  trans_bed$name <- NULL
  
  ## generate TSS statistics
  TSS_summary <- SummaryTranSite(gr = trans_bed,type = 'TSS',bin = bin)
  
  ## generate TES statistics
  # transcript clusters with polya were removed
  trans_bed_nopolya <- trans_bed[!(trans_bed$cluster_id %in% cluster_polya),]
  TES_summary <- SummaryTranSite(gr = trans_bed_nopolya,type = 'TES',bin = bin)
  
  save(TSS_summary,file = file.path(data_dir,'TSS_summary.RData'))
  save(TES_summary,file = file.path(data_dir,'TES_summary.RData'))
  
  ###########################################################################
  ### Step 3: Extract high confident TSS and TES
  ###########################################################################
  message('Step 3: Extract high confident TSS and TES')
  site_stats <- TSS_summary
  type <- ifelse('TSS' %in% colnames(site_stats),'TSS','TES')
  HC_TSS <- HCsite(site_stats = site_stats,cut_off = cut_off,cut_at = cut_at,
                   site_region = TSS_region,min_read = min_read,type = type)
  
  ## Get HC TES
  site_stats <- TES_summary
  type <- ifelse('TSS' %in% colnames(site_stats),'TSS','TES')
  HC_TES <- HCsite(site_stats = site_stats,cut_off = cut_off,cut_at = cut_at,
                   site_region = TES_region,min_read = min_read,type = type)
  
  save(HC_TSS,file = file.path(data_dir,'HC_TSS.RData'))
  save(HC_TES,file = file.path(data_dir,'HC_TES.RData'))
  
  ##########################################################################
  rm(trans_bed)
  rm(TSS_summary)
  rm(TES_summary)
  rm(HC_TSS)
  rm(HC_TES)
  gc()
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message('Done: ',Sys.time())
  message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units))
  
  
}