# SpliceWizResources

This repository contains links to Mappability Exclusion resources for SpliceWiz. It is intended for those who are running SpliceWiz on Bioconductor 3.13 or earlier.

Usage:
```{r}
library(SpliceWiz)
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
