
message('Loading R packages ... ')
suppressMessages(library("optparse"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("RTDBox"))

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = '*.R')) {
    #if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    #if(trace) cat("/n")
  }
}
sourceDir(path = 'R',encoding = 'UTF-8')
sourceDir(path = '/home/wguo/scratch/pantrans_isoseq/code/R',encoding = 'UTF-8')

option_list <- list(
  make_option(opt_str = c('-i','--input'),help = 'Input bed or gtf file'),
  make_option(opt_str = c('-f','--format'),help = 'Format of the input file. Options: "gtf" and "bed"'),
  make_option(opt_str = c('-o','--ouput'),help = 'Output csv file name to save the summary.')
)

opt_parser <- OptionParser(option_list=option_list)
# print_help(opt_parser)
opt <- parse_args(opt_parser)

file2read <- opt$input
f <- opt$format
file2save <- opt$ouput

if(is.null(f))
  f <- gsub('.*[.]','',file2read)

if(is.null(file2save)){
  file2save <- paste0(gsub('[.].*','',file2read),'_summary.csv')
}
  

gr <- import(file2read)
if(f=='bed'){
  gtf <- unlist(blocks(gr))
  gtf$gene_id <- gsub(';.*','',names(gtf))
  gtf$transcript_id <- gsub('.*;','',names(gtf))
  gtf$feature <- 'exon'
} else {
  if('type' %in% colnames(mcols(gr))){
    gtf <- gr[gr$type=='exon',]
  } else if('feature' %in% colnames(mcols(gr))){
    gtf <- gr[gr$feature=='exon',]
  } else {
    stop('Please provide proper format of input file')
  }
}

message('Summarising the RTD ... ')
s <- rtdSummary(gtf)
write.csv(s,file = file2save)



  

