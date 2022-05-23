# SpliceWizResources

This repository contains flux-simulator `pro` and `par` files for the simulation dataset for SpliceWiz. Also, it contains SpliceWiz processBAM output files
of the simulation.

The included `pro` and `par` files were used by flux-simulator to simulate raw sequencing FASTA files, which were then converted to paired-end sequencing FASTQ files. The resulting FASTQ files were aligned to hg38 (Ensembl release 94) genome. Alignment BAM files were then processed by SpliceWiz's `processBAM()` command to generate the output files.

This repository also contains links to Mappability Exclusion resources for SpliceWiz. These files are intended for those who are running SpliceWiz on Bioconductor 3.13 or earlier.

Usage:
```{r}
devtools::install_github("alexchwong/SpliceWiz")
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
    MappabilityRef = "hg38.MappabilityExclusion.bed.Rds"
)
```
