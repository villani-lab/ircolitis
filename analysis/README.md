# Overview

We used R scripts for the analyses in this project.

Our analysis output files are provided in the `output/` folder.

The following scripts are most relevant for each figure in the published
article:

Figure 1

    R/14-figures.R

Figures 2, 3, 4, 5

    R/14-figures.R
    R/plot-analysis.R
    R/15-tcr.R
    R/17-blood.R

Figures 6, 7

    R/05-test-cluster-communication.R
    R/covarying-composition.R


# Raw and processed data files

scRNA-seq gene expression data is available at NCBI GEO accession [GSE206301].

[GSE206301]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206301

The NCBI GEO accession is made up of three different datasets:

1. Blood immune cells (`blood_cells`)
2. Tissue immune cells (`tissue_cells`)
3. Tissue epithelial and mesenchymal nuclei (`tissue_nuclei`)

The quickest way to get started is to download the `.h5ad` files that
correspond to each of the major cell lineages in our study.

For example, we can access the tissue immune cell data by downloading the
corresponding `.h5ad` files:

```bash
# Colon Tissue:
# B cells, CD4 T cells, CD8 T cells, Myeloid cells
url=ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206299/suppl
wget ${url}/GSE206299_ircolitis-tissue-b.h5ad.gz
wget ${url}/GSE206299_ircolitis-tissue-cd4.h5ad.gz
wget ${url}/GSE206299_ircolitis-tissue-cd8.h5ad.gz
wget ${url}/GSE206299_ircolitis-tissue-myeloid.h5ad.gz
```

```bash
# Epithelial and Mesenchymal nuclei
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE206300&format=file&file=GSE206300%5Fircolitis%2Dtissue%2Depithelial%2Eh5ad%2Egz
url=ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206300/suppl
wget ${url}/GSE206299_ircolitis-tissue-epithelial.h5ad.gz
```

```bash
# Blood PBMCs:
# B cells, CD4 T cells, CD8 T cells, Myeloid cells
url=ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206298/suppl
wget ${url}/GSE206299_ircolitis-blood-b.h5ad.gz
wget ${url}/GSE206299_ircolitis-blood-cd4.h5ad.gz
wget ${url}/GSE206299_ircolitis-blood-cd8.h5ad.gz
wget ${url}/GSE206299_ircolitis-blood-myeloid.h5ad.gz
```

Decompress the files and then read them with Python or R:

```bash
# Decompress .h5ad.gz to .h5ad
gzip --decompress GSE206299_ircolitis-tissue-b.h5ad.gz
```

Use the [anndata](https://pypi.org/project/anndata/) Python package to read the file:

```Python
import anndata
ad = anndata.read_h5ad("GSE206299_ircolitis-tissue-b.h5ad")
```

Use the [anndata](https://CRAN.R-project.org/package=anndata ) R package to read the file:

```R
library(anndata)
ad <- read_h5ad("GSE206299_ircolitis-tissue-b.h5ad")
```
 
