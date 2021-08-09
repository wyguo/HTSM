rtdSummary <- function(gr,big.mark = ","){

  type <- ifelse('feature' %in% colnames(mcols(gr)),'feature','type')
  if(!(type %in% c('feature','type')))
    stop('The object does not have a "feature" or "type" column of exon labels')
  exon <- gr[mcols(gr)[,type]=='exon',]
  ###multi-isoform genes
  mapping <- data.frame(TXNAME=exon$transcript_id,GENEID=exon$gene_id)
  mapping <- mapping[!duplicated(mapping),]
  genes_multi <- sum(table(mapping$GENEID) > 1)

  intron <- exon2intron(exon)

  gr.reduce <- reduce(exon)
  genome_cov <- sum(width(gr.reduce))

  exon_num <- NROW(exon)
  exon_len <- width(exon)
  exon_len_min <- min(exon_len)
  exon_len_max <- max(exon_len)
  exon_len_total <- sum(exon_len)
  exon_len_ave <- mean(exon_len)
  exon_len_median <- median(exon_len)

  intron_num <- NROW(intron)
  intron_len <- width(intron)
  intron_len_min <- min(intron_len)
  intron_len_max <- max(intron_len)
  intron_len_total <- sum(intron_len)
  intron_len_ave <- mean(intron_len)
  intron_len_median <- median(intron_len)

  exon.list <- GenomicRanges::split(exon,exon$transcript_id)

  genes_num <- length(unique(exon$gene_id))
  trans_num <- length(unique(exon$transcript_id))

  genes_mono_num <- table(elementNROWS(GenomicRanges::split(exon,exon$gene_id))==1)['TRUE']
  trans_mono_num <- table(elementNROWS(GenomicRanges::split(exon,exon$transcript_id))==1)['TRUE']


  trans_len <- sum(width(exon.list))
  trans_len_min <- min(trans_len)
  trans_len_max <- max(trans_len)
  trans_len_total <- sum(trans_len)
  trans_len_ave <- mean(trans_len)
  trans_len_median <- median(trans_len)
  trans_len_N50 <- N50(trans_len)
  trans_len_N90 <- N90(trans_len)

  trans_per_gene <- round(trans_num/genes_num,3)
  exon_per_trans <- round(exon_num/trans_num,3)
  intron_per_trans <- round(intron_num/trans_num,3)
  Number <- data.frame(
    `Genome covered bases` = as.numeric(genome_cov),
    `Gene number` = as.numeric(genes_num),
    `Multi-isoform gene number`= as.numeric(genes_multi),
    `Mono-exon gene number` = as.numeric(genes_mono_num),
    `Multi-exon gene number` = as.numeric(genes_num-genes_mono_num),
    `Transcript number`= as.numeric(trans_num),
    `Mono-exon transcript number`= as.numeric(trans_mono_num),
    `Multi-exon transcript number`= as.numeric(trans_num-trans_mono_num),
    `Transcript number per gene`= as.numeric(trans_per_gene),
    # `Transcript median length`=trans_len_median,
    # `Transcript min length`=trans_len_min,
    # `Transcript max length`=trans_len_max,
    `Exon number`= as.numeric(exon_num),
    `Exon number per transcript`= as.numeric(exon_per_trans),
    # `Exon total length`=exon_len_total,
    `Exon average length`= as.numeric(exon_len_ave),
    # `Exon median length`=exon_len_median,
    # `Exon min length`=exon_len_min,
    # `Exon max length`=exon_len_max,
    # `Exonic coverged bases on the genome`=genome_cov,
    `Intron number`= as.numeric(intron_num),
    `Intron number per transcript`= as.numeric(intron_per_trans),
    # `Intron total length`=intron_len_total,
    `Intron average length`= as.numeric(intron_len_ave),
    `Transcript N50`= as.numeric(trans_len_N50),
    `Transcript N90`= as.numeric(trans_len_N90),
    `Transcript average length (exonic)`= as.numeric(trans_len_ave),
    `Transcript total length (exonic)`= as.numeric(trans_len_total),
    # `Intron median length`=intron_len_median,
    # `Intron min length`=intron_len_min,
    # `Intron max length`=intron_len_max,
    check.names = F
  )


  Number <- prettyNum(Number,big.mark=big.mark,scientific=FALSE)
  Number <- data.frame(Number,check.names = F)
  Number
}
