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
[GSE206301]. Please see [analysis/README.md](analysis) for instructions on how
to access the files.

[GSE206301]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206301

