#' Transcript start site (TSS) and end site (TES) analysis
#' @param data_dir The data directory of TAMA output.
#' @param cut_off A numeric value of BH adjsuted p-value (fdr), p-value or probablity cut-off.
#' @param cut_at The \code{cut_off} is applied to "fdr" (adjustd p-value), "pval" (p-value) or "prob" (probability).
#' @param TSS_region,TES_region An integer of enriched TSS/TES wobbling region. TSSs/TESs locate within enriched
#' site+/-site_region are all treated as high confident sites. Default: 50.
#' @param min_read The minimum reads in the window around TSS/TES to support high confident sites. Default: 2.
#' @param bin An integer of window size to aggregate the reads upstream and downstream a TSS/TES.
#' The total reads in [site-bin,site+bin] will be calculated. Default: 5.
#' @return Results are saved in \code{data_di}
#' 

TSanalysis <- function(data_dir,cut_off = 0.05,cut_at = 'fdr',
                       TSS_region = 50,TES_region = 50,min_read = 2,bin = 5){

  start.time <- Sys.time()
  message('|=======> Transcript start and end site analysis: ',Sys.time(),' <=======|')
  ###########################################################################
  ### Step 1: Filter transcripts with polyA
  ###########################################################################
  ##---> get polya transcript
  load(file.path(data_dir,'input_files.RData'))
  
  message('Step 1: Get reads with potential poly A truncation')
  ### get read id of polya
  cluster_polya=lapply(input_files$samples, function(s){
    ## input polya
    file2read <- input_files[s,'polya']
    polya <- read.table(file2read,
                        header=T,sep="\t",quote = "\"",dec=".")
    unique(polya$cluster_id)
  })
  cluster_polya <- unique(unlist(cluster_polya))
  
  reads_bed <- lapply(input_files$samples, function(s){
    ## input the reads
    file2read <- input_files[s,'reads']
    bed <- import(file2read)
    bed$transcript_id <- paste0(s,'_',gsub(';.*','',bed$name))
    bed$gene_id <- gsub('[.].*','',bed$transcript_id)
    bed$cluster_id <- gsub('.*;','',bed$name)
    bed$name <- NULL
    bed
  })
  reads_bed <- do.call(c,reads_bed)
  
  ###########################################################################
  ### Step 2: Generate statistics of TSS and TES
  ###########################################################################
  message('\nStep 2: Transcript start site (TSS) analysis')
  ## generate TSS statistics
  TSS_summary <- SummaryTranSite(gr = reads_bed,type = 'TSS',bin = bin)
  
  message('\nStep 3: Transcript end site (TES) analysis')
  ## generate TES statistics
  # transcript clusters with polya were removed
  reads_bed_nopolya <- reads_bed[!(reads_bed$cluster_id %in% cluster_polya),]
  TES_summary <- SummaryTranSite(gr = reads_bed_nopolya,type = 'TES',bin = bin)

  save(TSS_summary,file = file.path(data_dir,'TSS_summary.RData'))
  save(TES_summary,file = file.path(data_dir,'TES_summary.RData'))

  ###########################################################################
  ### Step 3: Extract high confident TSS and TES
  ###########################################################################
  message('\nStep 4: Extract high confident TSS and TES')
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
  rm(reads_bed)
  rm(TSS_summary)
  rm(TES_summary)
  rm(HC_TSS)
  rm(HC_TES)
  gc()

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message('\nDone: ',Sys.time())
  message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units),'\n')


}
