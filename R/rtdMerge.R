#' Merge inference RTD to reference RTD
#' @inf_file The gtf file of a inference RTD,e.g. a RTD assembled from short reads
#' @ref_file The gtf file of a high confident reference RTD, e.g. a RTD assembled from long reads. The transcripts of the 
#' reference RTD will be all kept in the final results.
#' @prefix_ref A prefix attached to the gene and transcript ids to distinguish their origin of "ref" RTD. 
#' Default: "ref"
#' @prefix_inf A prefix attached to the gene and transcript ids to distinguish their origin of "inf" RTD. 
#' Default: "inf"
#' @chimeric_tolerance If the loci overlap of two genes < chimeric_tolerance, they are treated as two seperate gene models.
#' @data_dir The directory to save the results
#' @return A merged RTD saved in the \code{data_di}
#' @example 
#' data_dir <- 'data/Morex_V3'
#' inf_file <- 'data/Morex_V3/Morex_V3_short_reads_assembly.gtf'
#' ref_file <- 'data/Morex_V3/Morex_V3_Iso_tama_merged50.gtf'
#' prefix_ref <- 'isoseq'
#' prefix_inf <- 'rnaseq'
#' chimeric_tolerance <- n
#' library(rtracklayer)
#' library(GenomicRanges)
#' rtdMerge(inf_file = inf_file,ref_file = ref_file,
#'          prefix_ref = 'isoseq',prefix_inf = 'rnaseq',
#'          chimeric_tolerance = 10,data_dir = data_dir)

rtdMerge <- function(inf_file,
                     ref_file,
                     prefix_ref='ref',
                     prefix_inf='inf',
                     chimeric_tolerance = 0.05,
                     data_dir){
  start.time <- Sys.time()
  message('Step 1: Load gtf files')
  ## load gtf files
  ref <- import(ref_file)
  inf <- import(inf_file)
  
  if(grepl(pattern = '[.]bed',inf_file)){
    inf <- unlist(blocks(inf))
    inf$gene_id <- gsub(';.*','',names(inf))
    inf$transcript_id <- gsub('.*;','',names(inf))
    inf$type <- 'exon'
    names(inf) <- NULL
    inf <- sort(inf,by=~seqnames + start + end)
  }
  
  if(grepl(pattern = '[.]bed',ref_file)){
    ref <- unlist(blocks(ref))
    ref$gene_id <- gsub(';.*','',names(ref))
    ref$transcript_id <- gsub('.*;','',names(ref))
    ref$type <- 'exon'
    names(ref) <- NULL
    ref <- sort(ref,by=~seqnames + start + end)
  }
  
  ref <- cleanMcols(ref)
  inf <- cleanMcols(inf)
  
  ref <- ref[ref$type=='exon',]
  inf <- inf[inf$type=='exon',]
  
  ## add prefix to gene ids
  ref$gene_id <- paste0(prefix_ref,'_',ref$gene_id)
  ref$transcript_id <- paste0(prefix_ref,'_',ref$transcript_id)
  ref$source <- prefix_ref
  
  inf$gene_id <- paste0(prefix_inf,'_',inf$gene_id)
  inf$transcript_id <- paste0(prefix_inf,'_',inf$transcript_id)
  inf$source <- prefix_inf
  
  message('Step 2: Identify novel transcript in the inference RTD')
  inf_trans_list <- list()
  
  message('  -> Process multi-exon transcripts')
  ####################################################
  ###---> Multi-exon transcripts
  
  ## generate unique labels of SJ chain
  intron_ref <- exon2intron(ref)
  intron_ref_sj <- intron_ref$label 
  
  intron_inf <- exon2intron(inf)
  intron_inf_sj <- intron_inf$label
  
  ## get the transcript names with novel SJ labels
  idx <- which(!(intron_inf_sj %in% intron_ref_sj))
  trans_multiexon_novel <- unique(intron_inf[idx]$transcript_id)

  
  ## multi-exon transcripts in inf which are not novel
  trans_multiexon_filtered <- setdiff(names(intron_inf_sj),trans_multiexon_novel)
  inf_trans_list <- c(inf_trans_list,trans_multiexon_filtered=list(trans_multiexon_filtered))
  inf_trans_list <- c(inf_trans_list,trans_multiexon_novel=list(trans_multiexon_novel))
  
  
  message('  -> Process mono-exon transcripts')
  ## get mono exon transcripts
  ref_mono <- getMonoExonTrans(ref)
  inf_mono <- getMonoExonTrans(inf)
  
  ## check whether inf mono exon trans has overlap of exons of ref, introic is excluded. 
  trans_mono <- unique(inf_mono$transcript_id)
  hits <- suppressWarnings(findOverlaps(query = inf_mono,
                                        subject = ref,ignore.strand=FALSE))
  
  ## get novel mono exon trans with no overlap
  trans_overlap <- inf_mono[queryHits(hits),]$transcript_id
  trans_monoexon_novel <- setdiff(trans_mono,trans_overlap)
  trans_monoexon_filtered <- setdiff(trans_mono,trans_monoexon_novel)
  
  inf_trans_list <- c(inf_trans_list,trans_monoexon_novel=list(trans_monoexon_novel))
  inf_trans_list <- c(inf_trans_list,trans_monoexon_filtered=list(trans_monoexon_filtered))
  
  ## novel transcripts in inf rtd, multi-exon + mono-exon
  trans_novel<- union(trans_multiexon_novel,trans_monoexon_novel)
  inf2keep <- inf[inf$transcript_id %in% trans_novel]
  
  message('Step 3: Merge inf to ref RTD and rename transcript and gene ids')
  rtd <- suppressWarnings(c(ref,inf2keep))
  rtd <- sort(rtd,by=~seqnames+start+end)
  rtd$gene_id_old <- rtd$gene_id
  rtd$transcript_id_old <- rtd$transcript_id
  # export_gtf(gr = rtd,file2save = file.path(data_dir,'tmp_gtf.gtf'))
  
  message('  -> Process intronic genes')
  ## overlap between gene ranges and introns of collapsed exons
  gene_range <- getGeneRange(rtd)
  gene_collapse <- getGeneCollapse(rtd)
  gene_collapse$transcript_id <- gene_collapse$gene_id
  gene_collapse$type <- 'exon'
  intron <- exon2intron(gene_collapse)
  
  hits <- findOverlaps(query = gene_range,subject = intron)
  idx1 <- pintersect(gene_range[queryHits(hits)],intron[subjectHits(hits)])
  idx2 <- gene_range[queryHits(hits)]
  intronic <- idx2[idx1==idx2]
  genes_intronic <- unique(intronic$gene_id)
  
  message('  -> Process chimeric check')
  ###---> find overlap gene ranges and gene overlap range
  hits <- findOverlaps(query = gene_range,subject = gene_range)
  idx <- which(queryHits(hits)!=subjectHits(hits))
  hits <- hits[idx,]
  gr1 <- gene_range[queryHits(hits)]
  gr2 <- gene_range[subjectHits(hits)]
  overlaps <- pintersect(gr1,gr2)
  p1 <- width(overlaps)/width(gr1)
  p2 <- width(overlaps)/width(gr2)
  
  idx <- which(p1 < chimeric_tolerance & p2 < chimeric_tolerance)
  genes_nochimeric <- union(gr1[idx]$gene_id,gr2[idx]$gene_id)
  
  genes2recover <- union(genes_intronic,genes_nochimeric)
  
  message('  -> Process overlap genes models')
  ###---> find overlap gene ranges and gene overlap range
  gene_overlap <- getOverlapRange(rtd)
  hits <- findOverlaps(query = gene_range,subject = gene_overlap)
  genes_rename <- gene_overlap[subjectHits(hits)]$gene_id
  names(genes_rename) <- gene_range[queryHits(hits)]$gene_id
  
  recover <- paste0(seqnames(gene_range),':',start(gene_range),';',strand(gene_range))
  names(recover) <- gene_range$gene_id
  genes_rename[genes2recover] <- recover[genes2recover]

  message('  -> Rename gene id')
  rtd$gene_id <- genes_rename[rtd$gene_id_old]
  
  message('  -> Rename transcript id')
  mapping <- data.frame(seqnames=seqnames(rtd),
                        transcript_id_old=rtd$transcript_id,
                        gene_id=rtd$gene_id,
                        source=rtd$source)
  mapping <- mapping[!duplicated(mapping),]
  rownames(mapping) <- mapping$transcript_id_old
  mapping <- mapping[order(mapping$gene_id,decreasing = F),]
  mapping$gene_id <- rep(1:length(unique(mapping$gene_id)),rle(mapping$gene_id)$lengths)
  
  mapping$transcript_id <- sequence(rle(mapping$gene_id)$lengths)
  mapping$gene_id <- paste0(mapping$seqnames,'G',mapping$gene_id)
  mapping$transcript_id <- paste0(mapping$source,'_',mapping$gene_id,'.',mapping$transcript_id)
  
  rtd$transcript_id <- mapping[rtd$transcript_id_old,'transcript_id']
  rtd$gene_id <- mapping[rtd$transcript_id_old,'gene_id']
  
  message('Step 4: Summary the merged RTD')
  
  ###---> transcript per gene
  n <- as.character(c(1:20,'>20'))
  diversity <- lapply(list(ref=ref,inf=inf,merged=rtd), function(gr){
    trans2gene <- data.frame(TXNAME=gr$transcript_id,GENEID=gr$gene_id)
    trans2gene <- trans2gene[!duplicated(trans2gene),]
    idx <- table(table(trans2gene$GENEID))
    if(max(as.numeric(names(idx))) > 20){
      idx1 <- idx[as.numeric(names(idx))<=20]
      idx2 <- sum(idx[as.numeric(names(idx))>20])
      names(idx2) <- '>20'
      idx <- c(idx1,idx2)
    } 
    num <- rep(0,length(n))
    names(num) <- n
    num[names(idx)] <- idx
    num
  })
  
  transpergene <- do.call(cbind,diversity)
  transpergene <- data.frame(TransPerGene=rownames(transpergene),transpergene,row.names = NULL)
  colnames(transpergene) <- c('TransPerGene',prefix_ref,prefix_inf,'merged')
  write.csv(transpergene,
            file=file.path(data_dir,paste0(prefix_ref,' and ', prefix_inf,' merged rtd transcript per gene.csv')),
            row.names = F)
  
  data2plot <- rbind(
    data.frame(TransPerGene=transpergene[,1],Number=transpergene[,2],RTDs=colnames(transpergene)[2]),
    data.frame(TransPerGene=transpergene[,1],Number=transpergene[,3],RTDs=colnames(transpergene)[3]),
    data.frame(TransPerGene=transpergene[,1],Number=transpergene[,4],RTDs=colnames(transpergene)[4])
  )
  data2plot$TransPerGene <- factor(data2plot$TransPerGene,levels = transpergene$TransPerGene)
  data2plot$RTDs <- factor(data2plot$RTDs,levels = c(prefix_ref,prefix_inf,'merged'))
  
  g <- ggplot(data2plot,aes(x=TransPerGene,y=Number,fill=RTDs))+
    geom_bar(stat='identity',position = position_dodge())+
    geom_text(aes(label = Number,color=RTDs),position = position_dodge(width = 1),
              hjust = -0.2,angle=90,size=2.5)+
    theme_bw()+
    labs(title = 'Transcript number per gene',x='Transcripts in a gene',y='Genen number')+
    coord_cartesian(ylim = c(0,max(data2plot$Number)*1.1))
  
  png(file.path(data_dir,'Transcript number per gene.png'),width = 11,height = 5,res = 300,units = 'in')
  print(g)
  dev.off()
  
  pdf(file.path(data_dir,'Transcript number per gene.pdf'),width = 11,height = 5)
  print(g)
  dev.off()
  
  ###---> basic statistics
  stat_merged <- rtdSummary(gr = rtd)
  stat_ref <- rtdSummary(gr = ref)
  stat_inf <- rtdSummary(gr = inf)
  stat <- cbind(stat_ref,stat_inf,stat_merged)
  colnames(stat) <- c(prefix_ref,prefix_inf,'merged')
  
  write.csv(stat,file=file.path(data_dir,paste0(prefix_ref,' and ', prefix_inf,' merged rtd summary.csv')))
  message('Step 5: Save the results')
  
  export_gtf(gr = rtd,
             file2save = file.path(data_dir,paste0(prefix_ref,' and ', prefix_inf,' merged rtd.gtf')))
  save(rtd,file = file.path(data_dir,paste0(prefix_ref,' and ', prefix_inf,' merged rtd.RData')))
  save(inf_trans_list,file=file.path(data_dir,'inf_trans_list.RData'))
  
  #######################
  ###---> save to bed file
  exon <- rtd
  mapping <- data.frame(gene_id = exon$gene_id,transcript_id = exon$transcript_id)
  mapping <- mapping[!duplicated(mapping),]
  rownames(mapping) <- mapping$transcript_id
  
  exon <- GenomicRanges::split(exon,exon$transcript_id)
  bed <- rtracklayer::asBED(exon)

  blocks <- bed$blocks
  bed$blocks <- NULL
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
  write.table(df, file = file.path(data_dir,paste0(prefix_ref,' and ', prefix_inf,' merged rtd.bed')), 
              sep = "\t", col.names = FALSE,
              row.names = FALSE, quote = FALSE, na = ".")
  
  message('Done: ',Sys.time())
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units),'\n')
}