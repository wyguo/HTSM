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

Run pipeline with R code
------------------------

```
library("RTDBox")
library("rtracklayer")
library("GenomicRanges")
library("Biostrings")
library("tidyr")

genome_fasta <- '/home/wguo/projects/jhi/barley/202012_panTranscriptome/genomes/ftp.ipk-gatersleben.de/SHAPE_I/Morex_V2_pseudomolecules_and_unplaced_scaffolds_ENA.fasta'
data_dir <- '/home/wguo/scratch/pantrans_isoseq/result/test/morex'
sj_overhang <- 10
cut_off <- 0.05
cut_at <- 'fdr'
min_read <- 2
TSS_region <- 50
TES_region <- 50
bin <- 5

###---> SJ analysis
SJanalysis(data_dir = data_dir,genome_fasta = genome_fasta,sj_overhang = sj_overhang)

###---> TSS/TES analysis
TSanalysis(data_dir = data_dir,cut_off = cut_off,cut_at = cut_at,
           min_read = min_read,TSS_region = TSS_region,
           TES_region = TES_region,bin = bin)

###---> Filter RTD
RTDfilter(data_dir = data_dir)

```

Run pipeline with bash script
-----------------------------

```
genome_fasta='/home/wguo/projects/jhi/barley/202012_panTranscriptome/genomes/ftp.ipk-gatersleben.de/SHAPE_I/Morex_V2_pseudomolecules_and_unplaced_scaffolds_ENA.fasta'
data_dir='/home/wguo/scratch/pantrans_isoseq/result/test/morex'
sj_overhang=10
cut_off=0.05
cut_at='fdr'
min_read=2
TSS_region=50
TES_region=50
bin=5

## help
# Rscript /mnt/shared/scratch/wguo/apps/conda/envs/isoseq/lib/R/library/RTDBox/scripts/isofilter.R --help

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
