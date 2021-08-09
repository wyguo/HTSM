
N50 <- function(lengths){
  lengths <- as.numeric(sort(lengths, decreasing=TRUE))
  index <- which(cumsum(lengths) / sum(lengths) >= 50/100)[1]
  return(lengths[index])
}

N90 <- function(lengths){
  lengths <- as.numeric(sort(lengths, decreasing=TRUE))
  index <- which(cumsum(lengths) / sum(lengths) >= 90/100)[1]
  return(lengths[index])
}