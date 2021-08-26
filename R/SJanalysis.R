#' Splice Junction (SJ) analysis
#' @param data_dir The data direcotory of TAMA output
#' @param genome_fasta The name of genome fasta sequence file. If the file is not in the working directory, please
#' provide the full path.
#' @param sj_overhang An integer of SJ overhang. A SJ is filtered if it has mismatches at its left (-sj_overhang)
#' or right (+sj_overhang) exonic regions.
#' @param ref_rtd An gtf file of reference transcript annotation. If a canonical SJ in the isoseq assembly has 
#' a match in the reference transcript annotation, it will be kept even though it has mismatches in the overhangs.
#'
#' genome_fasta <- 'data/isoseq/Morex_V2_pseudomolecules_and_unplaced_scaffolds_ENA.fasta'
SJanalysis <- function(data_dir,genome_fasta = NULL,sj_overhang = 10,ref_rtd=NULL){
  start.time <- Sys.time()

  if(is.null(genome_fasta))
    stop('Please provide the genome sequence fasta file')

  message('|==>   Splice Junction (SJ) analysis (ca. 6min): ',Sys.time(),'   <==|')
  message('Step 1: Prepare the local density error')
  #######################################################################
  ## Step 1: Prepare the local density error
  #######################################################################

  ###---> read the transcript model for each read
  file2read <- list.files(path = data_dir,
                          pattern = '_local_density_error.txt$',
                          full.names = T,
                          recursive = T)
  if(!file.exists(file2read))
    stop('File "*_local_density_error.txt" is missing in the data directory:\n',data_dir)

  ## only read column 1 and 12
  colClasses <- rep("NULL",13)
  colClasses[c(1,12)] <- "character"
  density_error <- read.delim(file2read,colClasses = colClasses)

  ## Filter read with "na" sj_error_simple
  idx <- which(density_error$sj_error_simple=='na')
  density_error <- density_error[-idx,]

  ## separate sj_error_simple to multiple rows with deliminate ";"
  density_error <- separate_rows(density_error, sj_error_simple,sep = ';')
  density_error$num_introns <- sequence(rle(density_error$cluster_id)$length)

  ###---> read transcript ids
  file2read <- list.files(path = data_dir,
                          pattern = '_collasped_trans_read.bed$',
                          full.names = T,
                          recursive = T)
  if(!file.exists(file2read))
    stop('File "*_collasped_trans_read.bed" is missing in the data directory:\n',data_dir)

  ## only read column 4
  colClasses <- rep("NULL",12)
  colClasses[4] <- "character"

  trans <- read.delim(file2read, header=FALSE,colClasses = colClasses)
  trans <- as.vector(t(trans))

  ## get cluster and transcript ids
  name.idx <- gsub(';.*','',trans)
  names(name.idx) <- gsub('.*;','',trans)

  idx <- which(density_error$cluster_id %in% names(name.idx))
  density_error <- density_error[idx,]

  ## Add transcript id to density error file
  density_error$transcript_id <- name.idx[density_error$cluster_id]
  density_error <- DataFrame(density_error)

  ## transcript id + intron numbering as row names
  rownames(density_error) <- paste0(density_error$transcript_id,'.',density_error$num_introns)
  save(density_error,file=file.path(data_dir,'density_error.RData'))

  message('Step 2: Extract splice junction motifs')
  #######################################################################
  ## Step 2: Extract splice junction motifs
  #######################################################################

  # This is a bed12 format file containing the final collapsed version of your transcriptome.
  file2read <- list.files(path = data_dir,
                          pattern = '_collasped.bed$',
                          full.names = T,
                          recursive = T)
  if(!file.exists(file2read))
    stop('File "*_collasped.bed" is missing in the data directory:\n',data_dir)

  collasped_bed <- import(file2read)
  collasped_bed$gene_id <- gsub(';.*','',collasped_bed$name)
  collasped_bed$transcript_id <- gsub('.*;','',collasped_bed$name)
  save(collasped_bed,file=file.path(data_dir,'collasped_bed.RData'))

  ## get exons of the transcriptome
  exons <- unlist(blocks(collasped_bed))

  ## add gene and transcript ids
  exons$gene_id <- gsub(';.*','',names(exons))
  exons$transcript_id <- gsub('.*;','',names(exons))
  names(exons) <- NULL


  introns <- exon2intron(exons,sorted = F)

  ## get the coordiantes of first two bases of introns
  sj_start <- resize(introns, width=2, fix="start", use.names=T,ignore.strand=T)
  sj_start <- split(sj_start,seqnames(sj_start))

  ## get the coordinates of the last two bases of introns
  sj_end <- resize(introns, width=2, fix="end", use.names=T,ignore.strand=T)
  sj_end <- split(sj_end,seqnames(sj_end))

  ###---> read the genome sequence fasta file
  dna <- readBStringSet(filepath = genome_fasta)
  names(dna) <- gsub(' .*','',names(dna))

  chr <- levels(seqnames(exons))
  names(chr) <- chr
  dna <- dna[chr]

  ## extract sj motifs
  sj_motif <- lapply(chr, function(x){
    # message(x)
    gr_start <- sj_start[[x]]
    if(is.null(gr_start))
      return(NULL)

    ## extract motif of first two bases of introns
    m_start <- extractAt(x = dna[[x]],at = gr_start@ranges)

    # message(x)
    gr_end <- sj_end[[x]]
    if(is.null(gr_end))
      return(NULL)

    ## extract motif of last two bases of introns
    m_end <- extractAt(x = dna[[x]],at = gr_end@ranges)

    ## put them into a data frame, add information of gene id, transcript id, strand
    motif <- data.frame(gene_id=gr_start$gene_id,
                        transcript_id=gr_start$transcript_id,
                        strand=strand(gr_start),
                        num_introns=gr_start$num_introns,
                        sj=gr_start$label,
                        motif=paste0(as.vector(m_start),as.vector(m_end)))
    motif
  })
  # idx <- which(sapply(sj_motif, nrow) > 0)
  # sj_motif <- sj_motif[idx]
  sj_motif <- do.call(rbind,sj_motif)
  sj_motif <- DataFrame(sj_motif)

  ## transcript id + intron numbering as row names
  rownames(sj_motif) <- paste0(sj_motif$transcript_id,'.',sj_motif$num_introns)
  # head(sj_motif)

  idx <- which(strand(sj_motif)=='-')
  motif <- sj_motif[idx,'motif']
  motif_unique <- unique(motif)

  ## reverse complement
  motif_crev <- DNArevComplement(motif_unique)
  names(motif_crev) <- motif_unique
  sj_motif[idx,'motif'] <- motif_crev[motif]

  save(sj_motif,file=file.path(data_dir,'sj_motif.RData'))

  message('Step 3: Filter non-canonical SJs and SJs with mismatches at overhang +/-',sj_overhang)
  #######################################################################
  ## Step 3: Filter non-canonical SJs and SJs with errors at overhang +/-
  #######################################################################
  canonical <- c("GTAG","GCAG","ATAC","CTAC","CTGC","GTAT")
  
  ## ref_rtd
  if(!is.null(ref_rtd)){
    gtf_db <- import(ref_rtd)
    idx <- ifelse('type' %in% colnames(mcols(gtf_db)),'type','feature')
    gtf_db <- gtf_db[mcols(gtf_db)[,idx]=='exon',]
    intron_db <- exon2intron(gtf_db)
    intron_db <- unique(intron_db$label)
    rm(gtf_db)
  } else {
    intron_db <- NULL
  }

  ## match row names of density error and sj motif (transcript_id.intron.number)
  idx <- match(rownames(density_error),rownames(sj_motif))
  motif2add <- sj_motif[idx,c('sj','motif')]
  rownames(motif2add) <- rownames(density_error)
  density_error_motif <- cbind(density_error,motif2add)

  ####################################
  ## filter non-canonical motif
  density_error_motif_filtered <- density_error_motif[density_error_motif$motif %in% canonical,]

  ## Filter SJs with mismatch errors;
  p <- paste0(c(rep('_',sj_overhang),'>',rep('_',sj_overhang)),collapse = '')
  idx1 <- grep(p,density_error_motif_filtered$sj_error_simple)
  idx2 <- which(density_error_motif_filtered$sj %in% intron_db)
  idx <- union(idx1,idx2)
  
  density_error_motif_filtered <- density_error_motif_filtered[idx,]

  save(density_error_motif_filtered,file=file.path(data_dir,'density_error_motif_filtered.RData'))

  ## rtd include mono-exon transcripts, must use "not in" instead of "in"
  idx <- which(sj_motif$sj %in% density_error_motif_filtered$sj)

  ## trans filter
  trans_remove <- sj_motif$transcript_id[-idx]
  trans_keep <- setdiff(exons$transcript_id,trans_remove)
  sj_trans_filter <- list(trans_remove=trans_remove,trans_keep=trans_keep)
  save(sj_trans_filter,file=file.path(data_dir,'sj_trans_filter.RData'))

  # rtd <- exons[which(exons$transcript_id %in% trans_keep),]
  # rtd$feature <- 'exon'
  # save(rtd,file=file.path(data_dr,'rtd.RData'))

  ##########################################################################
  rm(density_error)
  rm(collasped_bed)
  rm(exons)
  rm(introns)
  rm(sj_start)
  rm(sj_end)
  rm(dna)
  rm(sj_motif)
  rm(density_error_motif_filtered)
  gc()


  message('Done: ',Sys.time())
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units),'\n')

}
