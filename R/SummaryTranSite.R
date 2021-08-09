#' Summary transcript start site (TSS) and transcript end site (TES)
#' @details  Count the number of TSS and TES on each chromosome.
#' @param gr A GRanges object.
#' @param type "TSS" for transcript start site and "TES" for transcript end site.
#' @param bin An integer of window size to aggregate the reads upstream and downstream of a TSS/TES.
#' The total reads in [site-bin,site+bin] will be calculated.
#' @return A data frame of TSS or TES count.
#' 
SummaryTranSite <- function(gr,type=c('TSS','TES'),bin=5){
  
  # start.time <- Sys.time()
  ##########################################################################
  
  message('Count transcript ',type, ' sites')
  fix <- ifelse(type=='TSS','start','end')
  
  ## get transcript site
  trans <- resize(gr,fix = fix,width = 1,ignore.strand=F)
  
  ## get gene models
  genes <- data.frame(seqnames=seqnames(trans),
                      strand=strand(trans),
                      site=start(trans),
                      genes=trans$gene_id,stringsAsFactors = F)
  genes <- genes[!duplicated(genes),]
  
  ## dictionary of gene names for look up
  genes_label <- genes$genes
  names(genes_label) <- paste0(genes$seqnames,';',genes$strand,';',genes$site)
  
  ## manipulate the transcript site
  
  ###---> get the read number 
  num <- split(end(trans),paste0(seqnames(trans),';',strand(trans)))
  num <- lapply(names(num), function(i){
    # message(i)
    x <- num[[i]]
    if(length(x)==0)
      return(NULL)
    n <- table(x)
    x <- data.frame(seqnames=i,x=as.numeric(names(n)),freq=as.vector(n))
    x
  })
  num <- do.call(rbind,num)
  colnames(num) <- c('seqnames',type,'read')
  
  ###---> add gene names
  num$label <- paste0(num$seqnames,';',num[,type])
  num$gene_id <- genes_label[num$label]
  
  
  ################################
  ## sliding the site, count the total reads in bin
  num_split <- split(num,num$seqnames)
  sliding <- lapply(names(num_split), function(i){
    # message(i)
    x <- num_split[[i]]
    site <- x[,type]
    n <- rep((2*bin+1),nrow(x))
    sliding <- data.frame(
      site=rep(site,each=(2*bin+1)),
      interval=sequence(n,from = site-bin)
    )
    sliding <- sliding[sliding$interval %in% site,]
    count <- x$read
    names(count) <- x[,type]
    sliding$count <- count[as.character(sliding$interval)]
    result <- DataFrame(site=site,
                        count=rowsum(x = sliding$count,group = sliding$site),
                        row.names = NULL)
    colnames(result) <- c('site','count')
    result$label <- x$label
    result
  })
  sliding <- do.call(rbind,sliding)
  rownames(sliding) <- sliding$label
  num[,'read_in_bin'] <- sliding[num$label,'count']
  
  ## remove lable
  # num$label <- NULL
  num <- separate(data = num,col = 'seqnames',
                  into = c('seqnames','strand'),sep = ';')
  
  ## add binomial testing statistics
  read_sum <- ave(num$read,num$gene_id,FUN = sum)
  read_mean <- ave(num$read,num$gene_id,FUN = mean)
  site_num <- ave(num$read,num$gene_id,FUN = length)
  read <- num$read
  
  num$gene_site_num <- site_num
  num$gene_read_sum <- read_sum
  num$gene_read_mean <- read_mean
  num$prob <- dbinom(x = read,size = read_sum,prob = 1/site_num)
  num$pval <- pbinom(q = read-1,size = read_sum,prob = 1/site_num,lower.tail = F)
  num$fdr <- p.adjust(p = num$pval,method = 'BH')
  
  num <- DataFrame(num)
  
  ##########################################################################
  # end.time <- Sys.time()
  # time.taken <- end.time - start.time
  # message('Done!!!')
  # message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units))

  
  num
}
