<a href="https://zenodo.org/badge/latestdoi/419842722"><img src="https://zenodo.org/badge/419842722.svg" alt="DOI"></a>

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
<img src="https://github.com/villani-lab/ircolitis/assets/209714/2299353e-cc51-4929-b0a6-b672842435ee"></img>
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

<table><tr><td><a href="https://villani.mgh.harvard.edu/ircolitis/app/?ds=a12_4_4_t4_cd8_1_2&gene=PDCD1">View Cell Clusters</a> :microscope:</a></td></tr></table>
</td>
</tr>
<tr>
<td>
<a href="https://villani.mgh.harvard.edu/ircolitis/gene-contrasts">
<img src="https://github.com/villani-lab/ircolitis/assets/209714/7ffd151c-3540-4002-9141-339aef2988f8"></img>
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

<table><tr><td><a href="https://villani.mgh.harvard.edu/ircolitis/gene-contrasts">View Gene Contrasts :microscope:</a></td></tr></table>
</td>
</tr>
</table>


## Read the source code &#x1F4BB;

This repository includes three main folders:

[analysis/R][R]
- The R source code for our analyses.

[analysis/output][output]
- The output results files from our analyses.

[website/][website]
- The HTML, CSS, and Javascript source code for an interactive website to view the data and results.

[R]: https://github.com/villani-lab/ircolitis/tree/main/analysis/R
[output]: https://github.com/villani-lab/ircolitis/tree/main/analysis/output
[website]: https://github.com/villani-lab/ircolitis/tree/main/website


## Download the data &#x1F4BE;

The raw and processed scRNA-seq gene expression files are available at NCBI GEO
[GSE206301].

Sequencing reads are available at dbGAP accession [phs003418.v1.p1].

Please see [analysis/output/README.md][output] for instructions on how
to access the files.

[GSE206301]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206301
[phs003418.v1.p1]: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003418.v1.p1

