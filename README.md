## High resolution Transcriptome using Single Molecule PacBio Iso-seq data (HTSM)

Install packages
----------------

```
###---> Install HTSM
install.packages('devtools') ## if not installed
devtools::install_github("wyguo/HTSM")

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
      + prefix_polya.txt

Output data
-----------

Output to run TAMA merge

  - file2merge.txt: file list to merge in TAMA merge, see "-f" in https://github.com/GenomeRIK/tama/wiki/Tama-Merge
  - prefix_collasped_filtered.bed: The filtered collapsed version of your RTD.
  - prefix_collasped_filtered_summary.csv: The summary of the filtered collapsed RTD.

Intermediate output in .RData object, which you can manipulate with R scripts. See https://github.com/GenomeRIK/tama/wiki/Tama-Collapse for TAMA collapse file information. 
  - collasped_bed.RData: the final collapsed RTD from TAMA collapse is saved in .RData object. 
  - density_error.RData: the prefix_local_density_error file in TAMA collapse. 
  - density_error_motif_filtered.RData: the prefix_local_density_error after finter non-canonical and error splice junctions. 
  - sj_motif.RData: Splice junction motifs.
  - sj_trans_filter.RData: Splice junction motifs after filter non-canonical and error splice junctions. 
  - TSS_summary.RData: The Binomial testing statistics of TSS. 
  - TES_summary.RData: The Binomial testing statistics of TES. 
  - HC_TSS.RData: TSS of high confidence.
  - HC_TES.RData: TES of high confidence.

Run pipeline with R code
------------------------

```
library("HTSM")
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
# Rscript /mnt/shared/scratch/wguo/apps/conda/envs/isoseq/lib/R/library/HTSM/scripts/isofilter.R --help
# Note: isofilter.R is in the scripts folder where you installed the HTSM package

Rscript /mnt/shared/scratch/wguo/apps/conda/envs/isoseq/lib/R/library/HTSM/scripts/isofilter.R \
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
