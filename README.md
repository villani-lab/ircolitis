Welcome! Here, we present our study of immunotherapy-related colitis (ircolitis).

## Read the paper :mortar_board: 

Please read and cite our original research article:

- Thomas MF, Slowikowski K, Manakongtreecheep K, Sen P, Tantivit J, Nasrallah M,
et al. <b><a rel="noopener" target="#"
href="https://doi.org/10.1101/2021.09.17.460868">Altered interactions between
circulating and tissue-resident CD8 T cells with the colonic mucosa define
colitis associated with immune checkpoint inhibitors.</a></b> bioRxiv. 2021. doi:10.1101/2021.09.17.460868


## Explore the data :microscope: 
<table>
<tr>
<td width="33%">
<a href="https://villani.mgh.harvard.edu/ircolitis/app/?ds=a12_4_4_t4_cd8_1_2&gene=PDCD1">
<img src="https://user-images.githubusercontent.com/209714/182650768-3d646624-6655-489d-b5ea-f4309d0d1931.png"></img>
</a>
</td>
<td>
<b>Cell Clusters</b>

Metadata variables and gene expression in two-dimensional embeddings.

Tissue immune cells:
- CD8 T cells, CD4 T cells, Myeloid cells, B cells

Tissue epithelial and mesenchymal nuclei:
- 23 types of epithelial and mesenchymal nuclei

Blood immune cells:
- CD8 T cells, CD4 T cells, Myeloid cells, B cells

<a href="https://villani.mgh.harvard.edu/ircolitis/app/?ds=a12_4_4_t4_cd8_1_2&gene=PDCD1">View Cell Clusters</a>
</td>
</tr>
<tr>
<td>
<a href="https://villani.mgh.harvard.edu/ircolitis/de?gene=IL26">
<img src="https://user-images.githubusercontent.com/209714/182650864-02e76b1e-cad5-4a36-ba10-668c8f956cdd.png"></img>
</a>
</td>
<td>
<b>Gene Contrasts</b>

Differential expression statistics for all genes:
- 3 contrasts
    - irColitis Case vs Control
    - Case PD1/CTLA4 vs Case PD1
    - Control PD1 vs Control None
- 9 major cell lineages
- 105 cell clusters

<a href="https://villani.mgh.harvard.edu/ircolitis/de?gene=IL26">View Gene Contrasts</a>
</td>
</tr>
</table>


## Read the source code &#x1F4BB;

This repository includes two folders:

[analysis/][analysis]
- The R and Python source code for all of the analysis.

[website/][website]
- The HTML, CSS, and Javascript source code for an interactive website to view the data and results.

[analysis]: https://github.com/villani-lab/ircolitis/tree/master/analysis
[website]: https://github.com/villani-lab/ircolitis/tree/master/website


## Download the data &#x1F4BE;

We provide:

- Analysis outputs
- Raw and processed data files


### Analysis outputs

Please find analysis outputs in the folder [analysis/output](analysis/output)

See [this page](https://villani.mgh.harvard.edu/ircolitis/tables.html) to learn
more about what is contained in each file.

The contents of this GitHub repo are permanently archived at Zenodo.


### Raw and processed data files

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
# Tissue B cells, CD4 T cells, CD8 T cells, Myeloid cells
url=ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206299/suppl
wget ${url}/GSE206299_ircolitis-tissue-b.h5ad.gz
wget ${url}/GSE206299_ircolitis-tissue-cd4.h5ad.gz
wget ${url}/GSE206299_ircolitis-tissue-cd8.h5ad.gz
wget ${url}/GSE206299_ircolitis-tissue-myeloid.h5ad.gz
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

