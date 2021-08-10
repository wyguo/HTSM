## RTDBox-isoseq

Install packages
----------------

```
###---> Install RTDBox
install.packages('devtools') ## if not installed
devtools::install_github("wyguo/RTDBox")

###---> packages of dependencies
bioconductor.package.list <- c('rtracklayer','GenomicRanges','GenomicRanges','Biostrings')
for(i in bioconductor.package.list){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if(!(i %in% rownames(installed.packages()))){
    message('Installing package: ',i)
    BiocManager::install(i,dependencies = T)
  } else next
}
```

Input data
----------

The following files are mandatory for the analysis

  - A genome fasta file
  - The directory of TAMA output. The directory must includes:
      + prefix_local_density_error.txt
      + prefix_collasped_trans_read.bed
      + prefix_collasped.bed
      + prefix__polya.txt

Output data
-----------

Output to run TAMA merge

  - file2merge.txt: file list to merge in TAMA merge, see "-f" in https://github.com/GenomeRIK/tama/wiki/Tama-Merge
  - prefix_collasped_filtered.bed: The filtered collapsed version of your RTD.
  - prefix_collasped_filtered_summary.csv: The summary of the filtered collapsed RTD.

Intermediate output in .RData object, which you can manipulate with R scripts. 
  - HC_TES.RData
  - HC_TSS.RData
  - TES_summary.RData
  - TSS_summary.RData
  - sj_trans_filter.RData
  - density_error_motif_filtered.RData
  - sj_motif.RData
  - collasped_bed.RData
  - density_error.RData


Run pipeline with R code
------------------------

```
library("RTDBox")
library("rtracklayer")
library("GenomicRanges")
library("Biostrings")
library("tidyr")

genome_fasta <- 'genome.fasta'
data_dir <- 'tama_output'
sj_overhang <- 10
cut_off <- 0.05
cut_at <- 'fdr'
min_read <- 2
TSS_region <- 50
TES_region <- 50
bin <- 5

###---> SJ analysis
# ?SJanalysis for help information

SJanalysis(data_dir = data_dir,genome_fasta = genome_fasta,sj_overhang = sj_overhang)

###---> TSS/TES analysis
# ?TSanalysis for help information

TSanalysis(data_dir = data_dir,cut_off = cut_off,cut_at = cut_at,
           min_read = min_read,TSS_region = TSS_region,
           TES_region = TES_region,bin = bin)

###---> Filter RTD
# ?RTDfilter for help information

RTDfilter(data_dir = data_dir)

```

Run pipeline with bash script
-----------------------------

```
genome_fasta='genome.fasta'
data_dir='tama_output'
sj_overhang=10
cut_off=0.05
cut_at='fdr'
min_read=2
TSS_region=50
TES_region=50
bin=5

## help
# Rscript /mnt/shared/scratch/wguo/apps/conda/envs/isoseq/lib/R/library/RTDBox/scripts/isofilter.R --help
# Note: isofilter.R is in the scripts folder where you installed the RTDBox package

Rscript /mnt/shared/scratch/wguo/apps/conda/envs/isoseq/lib/R/library/RTDBox/scripts/isofilter.R \
-d $data_dir \
-g $genome_fasta \
-s $sj_overhang \
-c $cut_off \
-a $cut_at \
-r $min_read \
-t $TSS_region \
-e $TES_region \
-b $bin


```
