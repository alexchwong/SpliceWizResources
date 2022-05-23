# Simulation dataset for SpliceWiz

This folder contains flux-simulator input files for generating the alternative 
splicing simulation. Also, it contains the output files of SpliceWiz's processBAM() 
function, performed on the BAM files that were produced by alignment of 
flux-simulator generated sequencing data.

To view the SpliceWiz results:

```{r}
library(SpliceWiz)

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

# Annotate the samples
colData(se)$Biology <- rep(c("A", "B"), each = 3)

# Use SpliceWiz's optimized filters
se.filtered <- se[applyFilters(se),]

# Limma-based differential analysis
require(limma)
limma_res <- ASE_limma(se.filtered, "Biology", "A", "B")
```
