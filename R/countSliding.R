#' Count the TSS or TES number in sliding window
#' @param x A data frame with a column of "TSS" or "TES" and a column of "count" shows the number
#' of TSSs or ETSs.
#' TES     count
#' 80077         2
#' 80078         2
#' 80079         3
#' 80081         4
#' 80082        18
#' @param w A odd number of window size for sliding.
#' @param type One of "TSS" and "TES".
#' 
countSliding <- function(x,wsize=11,type=c('TTS','TES')){
  type <- match.arg(type,c('TTS','TES'))
  site <- x[,type]
  n <- rep(wsize,nrow(x))
  sliding <- data.frame(
    site=rep(site,each=wsize),
    interval=sequence(n,from = site-(wsize-1)/2)
  )
  sliding <- sliding[sliding$interval %in% site,]
  count <- x$read
  names(count) <- x[,type]
  sliding$count <- count[as.character(sliding$interval)]
  result <- data.frame(site=site,
                       count=rowsum(x = sliding$count,group = sliding$site),
                       row.names = NULL)
  colnames(result) <- c(type,'count_in_window')
  result
}