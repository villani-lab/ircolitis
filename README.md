<a href="https://zenodo.org/badge/latestdoi/419842722"><img src="https://zenodo.org/badge/419842722.svg" alt="DOI"></a>

Welcome! Here, we present our study of immunotherapy-related colitis (ircolitis).

## Read the paper :mortar_board: 

Please read and cite our original research article:

- Thomas MF, Slowikowski K, Manakongtreecheep K, Sen P, Samanta N, Tantivit J, et al. <a rel="noopener" target="#" href="https://doi.org/10.1038/s41591-024-02895-x">Single-cell transcriptomic analyses reveal distinct immune cell contributions to epithelial barrier dysfunction in checkpoint inhibitor colitis.</a> <i>Nature Medicine</i>. 2024; 1–14. doi:10.1038/s41591-024-02895-x

<details>
    <summary>BibTeX</summary>
    <pre>
@ARTICLE{Thomas2024-lk,
  title     = "{Single-cell transcriptomic analyses reveal distinct immune cell
               contributions to epithelial barrier dysfunction in checkpoint
               inhibitor colitis}",
  author    = "Thomas, Molly Fisher and Slowikowski, Kamil and
               Manakongtreecheep, Kasidet and Sen, Pritha and Samanta, Nandini
               and Tantivit, Jessica and Nasrallah, Mazen and Zubiri, Leyre and
               Smith, Neal P and Tirard, Alice and Ramesh, Swetha and Arnold,
               Benjamin Y and Nieman, Linda T and Chen, Jonathan H and
               Eisenhaure, Thomas and Pelka, Karin and Song, Yuhui and Xu,
               Katherine H and Jorgji, Vjola and Pinto, Christopher J and
               Sharova, Tatyana and Glasser, Rachel and Chan, Puiyee and
               Sullivan, Ryan J and Khalili, Hamed and Juric, Dejan and Boland,
               Genevieve M and Dougan, Michael and Hacohen, Nir and Li, Bo and
               Reynolds, Kerry L and Villani, Alexandra-Chloé",
  journal   = "Nature Medicine",
  publisher = "Nature Publishing Group",
  pages     = "1--14",
  abstract  = "Immune checkpoint inhibitor (ICI) therapy has revolutionized
               oncology, but treatments are limited by immune-related adverse
               events, including checkpoint inhibitor colitis (irColitis).
               Little is understood about the pathogenic mechanisms driving
               irColitis, which does not readily occur in model organisms, such
               as mice. To define molecular drivers of irColitis, we used
               single-cell multi-omics to profile approximately 300,000 cells
               from the colon mucosa and blood of 13 patients with cancer who
               developed irColitis (nine on anti-PD-1 or anti-CTLA-4 monotherapy
               and four on dual ICI therapy; most patients had skin or lung
               cancer), eight controls on ICI therapy and eight healthy
               controls. Patients with irColitis showed expanded mucosal Tregs,
               ITGAEHi CD8 tissue-resident memory T cells expressing CXCL13 and
               Th17 gene programs and recirculating ITGB2Hi CD8 T cells.
               Cytotoxic GNLYHi CD4 T cells, recirculating ITGB2Hi CD8 T cells
               and endothelial cells expressing hypoxia gene programs were
               further expanded in colitis associated with anti-PD-1/CTLA-4
               therapy compared to anti-PD-1 therapy. Luminal epithelial cells
               in patients with irColitis expressed PCSK9, PD-L1 and
               interferon-induced signatures associated with apoptosis,
               increased cell turnover and malabsorption. Together, these data
               suggest roles for circulating T cells and epithelial–immune
               crosstalk critical to PD-1/CTLA-4-dependent tolerance and barrier
               function and identify potential therapeutic targets for
               irColitis. Single-cell multi-omic analysis of 300,000 cells from
               29 patients representing peripheral immune cells and colon
               mucosal immune, epithelial and mesenchymal cells reveals
               crosstalk between circulating and tissue-resident immune cells
               with epithelial cells in checkpoint inhibitor colitis and
               identifies potential therapeutic targets.",
  month     =  may,
  year      =  2024,
  doi       = "10.1038/s41591-024-02895-x",
  issn      = "1546-170X,1546-170X",
  language  = "en"
}
    </pre>
</details>

## Explore the data :microscope: 
<table>
<tr>
<td width="33%">
<a href="https://villani.mgh.harvard.edu/ircolitis/app/?ds=a12_4_4_t4_cd8_1_2&gene=PDCD1">
<img src="https://github.com/villani-lab/ircolitis/assets/209714/34b32b5e-9000-4c3f-b0a8-2b487fc19e89"></img>
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
<img src="https://github.com/villani-lab/ircolitis/assets/209714/5bfd2ffe-018b-4680-ab69-efb90cff0599"></img>
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

