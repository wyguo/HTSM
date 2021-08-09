#' Reverse complement of DNA 
#' @param exons A string or a vector of strings
#' @return A reversed complement string or vector of strings
DNArevComplement <- function(s){
  s <- IRanges::reverse(s)
  chartr("ATCG", "TAGC", s)
}
