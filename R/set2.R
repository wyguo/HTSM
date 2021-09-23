set2 <- function (x, y){
  x.only <- setdiff(x, y)
  xy <- intersect(x, y)
  y.only <- setdiff(y, x)
  results <- list(x.only = x.only, xy = xy, y.only = y.only)
  attributes(results) <- list(x.only = "x-x&y", xy = "x&y", 
                              y.only = "y-x&y")
  names(results) <- c("x.only", "xy", "y.only")
  return(results)
}