# SpliceWizResources

This repository contains the following resources:

* Simulation data for SpliceWiz
  * Input `pro` and `par` input files for flux-simulator
    * These can be used to generate raw sequencing data for a simulation of differential alternative splicing
  * SpliceWiz processBAM() output files run on alignment BAM files of the abovementioned simulation
* Mappability Exclusion resources for SpliceWiz

### Usage

To use the following resources for SpliceWiz:

```
git clone https://github.com/alexchwong/SpliceWizResources.git
cd SpliceWizResources
```

NB: this git repository is approximately 640 Mb

### To install SpliceWiz

Please refer to the installation instructions viewable at https://github.com/alexchwong/SpliceWiz

Briefly, in R:

```{r}
devtools::install_github("alexchwong/SpliceWiz")
```

## Simulation data for SpliceWiz

flux-simulator `pro` and `par` files for the simulation dataset for SpliceWiz. Also, it contains SpliceWiz processBAM output files of the simulation.

The included `pro` and `par` files were used by flux-simulator to simulate raw sequencing FASTA files, which were then converted to paired-end sequencing FASTQ files. The resulting FASTQ files were aligned to hg38 (Ensembl release 94) genome. Alignment BAM files were then processed by SpliceWiz's `processBAM()` command to generate the output files.

Please refer to the documentation for flux-simulator for use of these files. The Ensembl GRCh38 (release 94) genome / gene annotations were used as the reference for flux-simulator.

Also, it contains the output files of SpliceWiz's processBAM() 
function, performed on the BAM files that were produced by alignment of 
flux-simulator generated sequencing data.

### Usage (in R/RStudio)

```{r}
library(SpliceWiz)
setwd("./simulation")

# Generate reference file (hg38 v94)
# Requires Bioconductor 3.14 or higher
buildRef(
    reference_path = "./Reference",
    fasta = "AH65745",
    gtf = "AH64631",
    genome_type = "hg38"
)

# Collate the SpliceWiz output files into an experiment
expr <- findSpliceWizOutput("./pb_output")
collateData(expr, "Reference", "NxtSE")

# Import experiment as NxtSE
se <- makeSE("NxtSE", realize = TRUE)

# Ensures NxtSE knows the proper path names of COV files
covfile(se) <- expr$cov_file 

# Annotate the samples
colData(se)$Biology <- rep(c("A", "B"), each = 3)

# Use SpliceWiz's optimized filters
se.filtered <- se[applyFilters(se),]

# Limma-based differential analysis
require(limma)
res_limma <- ASE_limma(se.filtered, "Biology", "A", "B")

# Example heatmap
library(pheatmap)
mat <- makeMatrix(se.filtered, res_limma$EventName[1:20])
pheatmap(mat, annotation_col = as.data.frame(colData(se.filtered)))

# Example coverage plot
p <- plotCoverage(se.filtered, res_limma$EventName[1], 
    condition = "Biology", tracks = c("A", "B"), stack_tracks = TRUE)

as_ggplot_cov(p) # displays ggplot (static plot)
p$final_plot     # displays plotly object (interactive plot)
```

## Mappability Exclusion Resources for SpliceWiz

This repository also contains links to Mappability Exclusion resources for SpliceWiz. These files are intended for those who are running SpliceWiz on Bioconductor 3.13 or earlier.

### Usage (in R)

```{r}
library(SpliceWiz)

# Set path for SpliceWiz reference
reference_path = "./Reference"
FTP <- "ftp://ftp.ensembl.org/pub/release-94/"

buildRef(
    reference_path = reference_path,
    fasta = paste0(FTP, "fasta/homo_sapiens/dna/",
        "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
    gtf = paste0(FTP, "gtf/homo_sapiens/",
        "Homo_sapiens.GRCh38.94.chr.gtf.gz"),
    genome_type = "hg38",
    MappabilityRef = "Mappability/hg38.MappabilityExclusion.bed.Rds"
)
```
