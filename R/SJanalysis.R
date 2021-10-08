#' Splice Junction (SJ) analysis
#' @param data_dir The data direcotory of TAMA output
#' @param genome_fasta The name of genome fasta sequence file. If the file is not in the working directory, please
#' provide the full path.
#' @param sj_overhang An integer of SJ overhang. A SJ is filtered if it has mismatches at its left (-sj_overhang)
#' or right (+sj_overhang) exonic regions.
#' @param sjdatabase A tab separated file of splice junction information. If a canonical SJ in the isoseq assembly has 
#' a match in the database, it will be kept even though it has mismatches in the overhangs.
#' @return Results are saved in \code{data_di}
#' 
SJanalysis <- function(data_dir,genome_fasta = NULL,sj_overhang = 10,sjdatabase=NULL){
  start.time <- Sys.time()
  message('|==========> Splice junction analysis: ',Sys.time(),' <==========|')
  message('Step1: Check alignment error +/-',sj_overhang,' bases around splice junctions')
  
  ##############################################################
  ###---> sj error analysis
  file_error <- list.files(path = data_dir,
                           pattern = '_local_density_error.txt$',
                           full.names = T,
                           recursive = T)
  samples <- names(file_error) <- gsub('_.*','',basename(file_error))
  
  
  file_polya <- list.files(path = data_dir,
                           pattern = '_polya.txt$',
                           full.names = T,
                           recursive = T)
  names(file_polya) <- gsub('_.*','',basename(file_polya))
  
  file_reads <- list.files(path = data_dir,
                           pattern = '_trans_read.bed$',
                           full.names = T,
                           recursive = T)
  names(file_reads) <- gsub('_.*','',basename(file_reads))
  
  file_trans <- list.files(path = data_dir,
                           pattern = '_collasped.bed$',
                           full.names = T,
                           recursive = T)
  names(file_trans) <- gsub('_.*','',basename(file_trans))
  
  input_files <- data.frame(samples=samples,
                            sj_error=file_error[samples],
                            polya=file_polya[samples],
                            reads=file_reads[samples],
                            trans=file_trans[samples],
                            row.names = samples)
  
  if(nrow(input_files)==1){
    input_files$trans_merged <- input_files$trans
  } else {
    file_trans_merged <- list.files(path = data_dir,
                                    pattern = 'sample_merged.bed$',
                                    full.names = T,
                                    recursive = T)
    
    if(length(file_trans_merged)==0 & nrow(input_files)>0)
      stop('It seems that the dataset has multiple samples. Please run tama merge to merge
       transcripts from multipele samples.')
    
    input_files$trans_merged <- file_trans_merged
  }
  
  write.csv(input_files,file = file.path(data_dir,'input_files.csv'),row.names = F)
  save(input_files,file = file.path(data_dir,'input_files.RData'))
  
  ##################################################################
  ###---> get sj label
  # s <- 's1'
  sj_trans <- lapply(samples, function(s){
    file2read <- input_files[s,'trans']
    reads_bed <- import(file2read)
    # reads_bed$gene_id <- gsub(';.*','',reads_bed$name)
    # reads_bed$transcript_id <- gsub('.*;','',reads_bed$name)
    
    ## get exons of the transcriptome
    exons <- unlist(blocks(reads_bed))
    
    ## add gene and transcript ids
    exons$gene_id <- gsub(';.*','',names(exons))
    exons$transcript_id <- gsub('.*;','',names(exons))
    names(exons) <- NULL
    introns <- exon2intron(exons,sorted = F)
    label <- data.frame(gene_id=paste0(s,'_',introns$gene_id),
                        transcript_id=paste0(s,'_',introns$transcript_id),
                        intron_id=paste0(s,'_',introns$transcript_id,'.',introns$num_introns),
                        sj=introns$label)
    label[!duplicated(label),]
  })
  names(sj_trans) <- samples
  
  ##################################################################
  ###---> sj database
  ## sjdatabase
  if(!is.null(sjdatabase)){
    if(!file.exists(sjdatabase))
      stop("The sj database file does not exist in the data direcotry.")
    message('  ->Reading SJ database: ',basename(sjdatabase),'\n')
    sj <- read.table(sjdatabase, header=FALSE,sep = '\t',quote = "\"",dec='.',fill = T,comment.char = "")
    sj <- sj[,1:4]
    sj <- sj[!duplicated(sj),]
    colnames(sj) <- c('seqnames','start','end','strand')
    sj$strand[sj$strand==0] <- '*'
    sj$strand[sj$strand==1] <- '+'
    sj$strand[sj$strand==2] <- '-'
    intron_db <- paste0(sj$seqnames,':',sj$start,'-',sj$end,';',sj$strand)
    intron_db <- unique(intron_db)
  } else {
    intron_db <- NULL
  }
  
  ##################################################################
  ###---> sj error
  p <- paste0(c(rep('_',sj_overhang),'>',rep('_',sj_overhang)),collapse = '')
  
  sj2keep <- list()
  sj2remove <- list()
  for(s in samples){
    file2read <- input_files[s,'sj_error']
    density_error <- read.table(file2read,header = T,sep="\t",quote = "\"",dec=".")
    density_error <- density_error[,c("cluster_id","scaff_name","start_pos","end_pos","strand","sj_error_simple" )]
    
    ## Filter read with "na" sj_error_simple
    idx <- which(density_error$sj_error_simple=='na')
    density_error <- density_error[-idx,]
    
    ## separate sj_error_simple to multiple rows with deliminate ";"
    density_error <- separate_rows(density_error, sj_error_simple,sep = ';')
    density_error$num_introns <- sequence(rle(density_error$cluster_id)$length)
    # density_error$sample <- s
    
    
    ###########
    ## only read column 4
    file2read <- input_files[s,'reads']
    colClasses <- rep("NULL",12)
    colClasses[4] <- "character"
    
    trans <- read.table(file2read, header=FALSE,sep="\t",quote = "\"",dec=".",
                        colClasses = colClasses)
    trans <- as.vector(t(trans))
    
    ## get cluster and transcript ids
    name.idx <- gsub(';.*','',trans)
    names(name.idx) <- gsub('.*;','',trans)
    density_error$transcript_id <- paste0(s,'_',name.idx[density_error$cluster_id])
    density_error$intron_id <- paste0(density_error$transcript_id,'.',density_error$num_introns)
    
    ## add sj labels
    sj_db <- sj_trans[[s]]
    idx <- match(density_error$intron_id,sj_db$intron_id)
    density_error$sj <- sj_db$sj[idx]
    
    #########
    idx1 <- grep(p,density_error$sj_error_simple)
    idx2 <- which(density_error$sj %in% intron_db)
    idx2keep <- union(idx1,idx2)
    kep <- unique(density_error$sj[idx2keep])
    rem <- setdiff(density_error$sj,kep)
    sj2keep <- c(sj2keep,setNames(object = list(kep),nm = s))
    sj2remove <- c(sj2remove,setNames(object = list(rem),nm = s))
  }
  sj_error <- setdiff(Reduce(union,sj2remove),Reduce(union,sj2keep))
  
  message('Step2: Check nocanonical splice junctions')
  ###############################################################
  ###---> check canonical sj
  # if(!file.exists(file2read))
  #   stop('File "*_collasped.bed" is missing in the data directory:\n',data_dir)
  file2read <- input_files$trans_merged[1]
  trans_bed2filter <- import(file2read)
  trans_bed2filter$gene_id <- gsub(';.*','',trans_bed2filter$name)
  trans_bed2filter$transcript_id <- gsub('.*;','',trans_bed2filter$name)
  save(trans_bed2filter,file=file.path(data_dir,'trans_bed2filter.RData'))
  
  ## get exons of the transcriptome
  exons <- unlist(blocks(trans_bed2filter))
  
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
  
  message('  ->Reading reference genome: ',basename(genome_fasta),'\n')
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
  
  sj_motif <- do.call(rbind,sj_motif)
  
  ## transcript id + intron numbering as row names
  # rownames(sj_motif) <- paste0(sj_motif$transcript_id,'.',sj_motif$num_introns)
  
  idx <- which(sj_motif$strand=='-')
  motif <- sj_motif[idx,'motif']
  motif_unique <- unique(motif)
  
  ## reverse complement
  motif_crev <- DNArevComplement(motif_unique)
  names(motif_crev) <- motif_unique
  sj_motif[idx,'motif'] <- motif_crev[motif]
  
  canonical <- c("GTAG","GCAG","ATAC","CTAC","CTGC","GTAT")
  sj_noncanonical <- unique(sj_motif$sj[!(sj_motif$motif %in% canonical)])
  
  sj_filter <- list(sj_error=sj_error,sj_noncanonical=sj_noncanonical)
  save(sj_filter,file=file.path(data_dir,'sj_filter.RData'))
  
  message('Done: ',Sys.time())
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units),'\n')
  # return(sj_filter)
}
