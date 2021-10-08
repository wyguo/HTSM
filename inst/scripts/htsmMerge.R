message('Loading R packages ... ')
suppressMessages(library("optparse"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("HTSM"))

option_list <- list(
  make_option(opt_str = c('-i','--inf'),
              help = 'Input inference RTD gtf or bed file. Format will be determined from the file extension.'),
  make_option(opt_str = c('-r','--ref'),
              help = 'Input reference RTD gtf or bed file. Format will be determined from the file extension.'),
  make_option(opt_str = c('-n','--preinf'),default = 'inf',
              help = 'Prefix to append to transcript ids from inference RTD. Default: inf'),
  make_option(opt_str = c('-m','--preref'),default = 'ref',
              help = 'Prefix to append to transcript ids from reference RTD. Default: ref'),
  make_option(opt_str = c('-c','--chimeric'),type = 'double',default = 0.05,
              help = 'A percentage value of chimeric gene model tolerance. For two overlapped gene models, if
                the overlap length < c% of both gene lengths, they are treated as two seprate gene models,
                otherwise, they are given the same gene id.'),
  make_option(opt_str = c('-o','--output'),
              help = 'The output directory.')
)

opt_parser <- OptionParser(option_list=option_list)
# print_help(opt_parser)
opt <- parse_args(opt_parser)

file2read <- opt$input
f <- opt$format
file2save <- opt$ouput

if(is.null(opt$output)){
  data_dir <- getwd()
} else {
  data_dir <- opt$output
}

if(is.null(opt$inf)){
  print_help(opt_parser)
  stop('Please provide the inference RTD gtf file')
}

if(is.null(opt$ref)){
  print_help(opt_parser)
  stop('Please provide the reference RTD gtf file')
}

inf_file <- opt$inf
ref_file <- opt$ref

prefix_ref <- opt$preref
if(is.null(prefix_ref))
  prefix_ref <- 'htms'

prefix_inf <- opt$preinf
if(is.null(prefix_inf))
  prefix_inf <- 'rtd'

chimeric_tolerance <- opt$chimeric


rtdMerge(inf_file = inf_file,ref_file = ref_file,
         prefix_ref = prefix_ref,prefix_inf = prefix_inf,
         chimeric_tolerance = chimeric_tolerance,
         data_dir = data_dir)





  

