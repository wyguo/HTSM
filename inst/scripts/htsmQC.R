message('Loading R packages ... ')
suppressMessages(library("optparse"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("Biostrings"))
suppressMessages(library("tidyr"))
suppressMessages(library("HTSM"))

options(stringsAsFactors = F,scipen = 99)

# data_dir <- 'data/isoseq'
# genome_fasta <- 'data/isoseq/Morex_V2_pseudomolecules_and_unplaced_scaffolds_ENA.fasta'
# sj_overhang <- 10
# cut_off <- 0.05
# cut_at <- 'fdr'
# min_read <- 2
# TSS_region <- 50
# TES_region <- 50
# bin <- 5

option_list <- list(
  make_option(opt_str = c('-d','--directory'),type = NULL,default = NULL,
              help = 'The data directory of TAMA output'),
  make_option(opt_str = c('-g','--genome'),
              help = 'Genome sequence fasta file'),
  make_option(opt_str = c('-s','--sjOverhang'),type = 'integer',default = 10,
              help = 'Length of overhang upstream and downstrean of SJ to check match errors.
                Default: 10'),
  make_option(opt_str = c('-c','--cutoff'),type = 'double',default = 0.05,
              help = 'Cut-off of "prob", "pval" or "fdr" for TSS and TES enrichment. Default: 0.05'),
  make_option(opt_str = c('-a','--cutat'),type = 'character',default = 'fdr',
              help = 'Statistics to determine enriched TSS/TES. Options are "fdr" for BH adjusted 
                p-value, "pval" for p-value and "prob" for probablity. Default: fdr.'),
  make_option(opt_str = c('-t','--TSSregion'),type = 'integer',default = 50,
              help = 'The upstream and downstream region that a significantly enriched TSS 
                can wobble. Default: 50.'),
  make_option(opt_str = c('-e','--TESregion'),type = 'integer',default = 50,
              help = 'The upstream and downstream region that a significantly enriched TES 
                can wobble. Default: 50.'),
  make_option(opt_str = c('-b','--bin'),type = 'integer',default = 5,
              help = 'An integer of unstream and downstream window size to aggregate the reads 
                around a TSS/TES. The total reads in [site-bin,site+bin] will be calculated. 
                Default: 5.'),
  make_option(opt_str = c('-r','--minreads'),type = 'integer',default = 2,
              help = 'The minimum reads in the window around TSS/TES to support high confident 
                sites. Default: 2.'),
  make_option(opt_str = c('-j','--sjdatabase'),default = NULL,
              help = 'Optional. A tab separated file of splice junction (SJ) of which the first column is
                chromosome names, second column is start coordinate of SJ, third column is the end coordinate
                of SJ and fourth column is strand. For example, the SJ.out.tab file from STAR mapping can be
                used here. The first four columns will be taken as input for the analysis. Default: NULL')
)
opt_parser <- OptionParser(option_list=option_list)
# print_help(opt_parser)
opt <- parse_args(opt_parser)



if(is.null(opt$genome)){
  print_help(opt_parser)
  stop('Please provide the genome sequence fasta file')
}

if(is.null(opt$directory)){
  data_dir <- getwd()
} else {
  data_dir <- opt$directory
}

genome_fasta <- opt$genome
sj_overhang <- opt$sjOverhang
cut_off <- opt$cutoff
cut_at <- opt$cutat
TSS_region <- opt$TSSregion
TES_region <- opt$TESregion
bin <- opt$bin
min_read <- opt$minreads
sjdatabase <- opt$sjdatabase

cat('data_dir =',data_dir,'\n')
cat('genome_fasta =',genome_fasta,'\n')
cat('sj_overhang =',sj_overhang,'\n')
cat('cut_off =',cut_off,'\n')
cat('cut_at =',cut_at,'\n')
cat('TSS_region =',TSS_region,'\n')
cat('TES_region =',TES_region,'\n')
cat('bin =',bin,'\n')
cat('min_read =',min_read,'\n')
cat('sjdatabase =',sjdatabase,'\n')

############################################
###---> SJ analysis
SJanalysis(data_dir = data_dir,
           genome_fasta = genome_fasta,
           sj_overhang = sj_overhang,
           sjdatabase = sjdatabase)

###---> TSS/TES analysis
TSanalysis(data_dir = data_dir,
           cut_off = cut_off,
           cut_at = cut_at,
           min_read = min_read,
           TSS_region = TSS_region,
           TES_region = TES_region,
           bin = bin)

###---> filter the RTD
rtdFilter(data_dir = data_dir)

