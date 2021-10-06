
message('Loading R packages ... ')
suppressMessages(library("optparse"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("HTSM"))

# sourceDir <- function(path, trace = TRUE, ...) {
#   for (nm in list.files(path, pattern = '*.R')) {
#     #if(trace) cat(nm,":")
#     source(file.path(path, nm), ...)
#     #if(trace) cat("/n")
#   }
# }
# sourceDir(path = 'R',encoding = 'UTF-8')


option_list <- list(
  make_option(opt_str = c('-g','--gtf'),help = 'The gtf file name to save.'),
  make_option(opt_str = c('-b','--bed'),help = 'A bed file to convert gtf')
)

opt_parser <- OptionParser(option_list=option_list)
# print_help(opt_parser)
opt <- parse_args(opt_parser)

bed_file <- opt$bed
gtf_file <- opt$gtf

message('Converting the format')
gtf2bed(gtf2bed,bed_file = bed_file)


