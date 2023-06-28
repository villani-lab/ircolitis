-   [Overview](#overview)
-   [Single-cell expression data](#single-cell-expression-data)
-   [Pseudobulk expression data](#pseudobulk-expression-data)
    -   [Pseudobulk at the level of cell
        lineage](#pseudobulk-at-the-level-of-cell-lineage)
    -   [Pseudobulk at the level of cell
        cluster](#pseudobulk-at-the-level-of-cell-cluster)
-   [Analysis results](#analysis-results)
    -   [Cell cluster differential abundance statistics for Case vs
        Control](#cell-cluster-differential-abundance-statistics-for-case-vs-control)
    -   [Cell cluster differential abundance statistics for
        anti-PD-1/anti-CTLA-4 vs anti-PD-1 within Cases, not including
        Controls](#cell-cluster-differential-abundance-statistics-for-anti-pd-1anti-ctla-4-vs-anti-pd-1-within-cases-not-including-controls)
    -   [Cell cluster markers](#cell-cluster-markers)
    -   [Differential gene expression statistics for three
        contrasts](#differential-gene-expression-statistics-for-three-contrasts)
    -   [Cell-cell communication analysis with means of gene
        pairs](#cell-cell-communication-analysis-with-means-of-gene-pairs)
    -   [Cell-cell communication analysis with Spearman correlation of
        gene
        percentages](#cell-cell-communication-analysis-with-spearman-correlation-of-gene-percentages)
    -   [Correlation analysis of cluster
        abundances](#correlation-analysis-of-cluster-abundances)

Overview
========

Here, we provide the analysis output files that we used in our
publication:

-   Thomas MF, Slowikowski K, Manakongtreecheep K, Sen P, Tantivit J,
    Nasrallah M, et al. **Altered interactions between circulating and
    tissue-resident CD8 T cells with the colonic mucosa define colitis
    associated with immune checkpoint inhibitors.** bioRxiv. 2021.
    p. 2021.09.17.460868.
    <a href="doi:10.1101/2021.09.17.460868" class="uri">doi:10.1101/2021.09.17.460868</a>

Please cite our publication if you use this data:

    @UNPUBLISHED{Thomas2021,
      title    = "{Altered interactions between circulating and tissue-resident CD8
                  T cells with the colonic mucosa define colitis associated with
                  immune checkpoint inhibitors}",
      author   = "Thomas, Molly Fisher and Slowikowski, Kamil and
                  Manakongtreecheep, Kasidet and Sen, Pritha and Tantivit, Jessica
                  and Nasrallah, Mazen and Smith, Neal P and Ramesh, Swetha and
                  Zubiri, Leyre and Tirard, Alice and Arnold, Benjamin Y and
                  Nieman, Linda T and Chen, Jonathan H and Eisenhaure, Thomas and
                  Pelka, Karin and Xu, Katherine H and Jorgji, Vjola and Pinto,
                  Christopher J and Sharova, Tatyana and Glasser, Rachel and Chan,
                  Elaina Puiyee and Sullivan, Ryan J and Khalili, Hamed and Juric,
                  Dejan and Boland, Genevieve M and Dougan, Michael and Hacohen,
                  Nir and Reynolds, Kerry L and Li, Bo and Villani,
                  Alexandra-Chlo{\'e}",
      journal  = "bioRxiv",
      pages    = "2021.09.17.460868",
      month    =  sep,
      year     =  2021,
      language = "en",
      doi      = "10.1101/2021.09.17.460868"
    }

Please contact [Kamil Slowikowski](mailto:kslowikowski@mgh.harvard.edu)
with comments or questions about these files.

Single-cell expression data
===========================

scRNA-seq gene expression data is available at NCBI GEO accession
[GSE206301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206301),
which is a SuperSeries that is composed of the following SubSeries:

-   [GSE206299](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206299)
    Tissue immune cells (`tissue_cells`)
-   [GSE206300](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206300)
    Tissue epithelial and mesenchymal nuclei (`tissue_nuclei`)
-   [GSE206298](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206298)
    Blood immune cells (`blood_cells`)

The quickest way to get started is to download the `.h5ad` files that
correspond to each of the major cell lineages in our study.

For example, we can access the tissue immune cell data by downloading
the corresponding `.h5ad` files:

    # Colon Tissue:
    # B cells, CD4 T cells, CD8 T cells, Myeloid cells
    url=ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206299/suppl
    wget ${url}/GSE206299_ircolitis-tissue-b.h5ad.gz
    wget ${url}/GSE206299_ircolitis-tissue-cd4.h5ad.gz
    wget ${url}/GSE206299_ircolitis-tissue-cd8.h5ad.gz
    wget ${url}/GSE206299_ircolitis-tissue-myeloid.h5ad.gz

    # Epithelial and Mesenchymal nuclei
    url=ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206300/suppl
    wget ${url}/GSE206300_ircolitis-tissue-epithelial.h5ad.gz

    # Blood PBMCs:
    # B cells, CD4 T cells, CD8 T cells, Myeloid cells
    url=ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206298/suppl
    wget ${url}/GSE206298_ircolitis-blood-b.h5ad.gz
    wget ${url}/GSE206298_ircolitis-blood-cd4.h5ad.gz
    wget ${url}/GSE206298_ircolitis-blood-cd8.h5ad.gz
    wget ${url}/GSE206298_ircolitis-blood-myeloid.h5ad.gz

Decompress the files and then read them with Python or R:

    # Decompress .h5ad.gz to .h5ad
    gzip --decompress GSE206299_ircolitis-tissue-b.h5ad.gz

Use the [anndata](https://pypi.org/project/anndata/) Python package to
read the file:

    import anndata
    ad = anndata.read_h5ad("GSE206299_ircolitis-tissue-b.h5ad")

Use the [anndata](https://CRAN.R-project.org/package=anndata) R package
to read the file:

    library(anndata)
    ad <- read_h5ad("GSE206299_ircolitis-tissue-b.h5ad")

We provide raw counts and log2CPM values in each of the `.h5ad` files in
the [anndata](https://anndata.readthedocs.io/en/latest/) format.

Within each file, the `obs` table contains the relevant metadata
(e.g. cell cluster names, PC scores, UMAP coordinates, etc.).

The TCR and BCR sequencing data is also contained in the `obs` table.
(Look for columns like `TRA_cdr3` and `IGH_cdr3`.)

Read the data in R:

    a1 <- anndata::read_h5ad("tissue-cd8.h5ad")
    a1

    ## AnnData object with n_obs × n_vars = 25341 × 28165
    ##     obs: 'cell', 'channel', 'case', 'class', 'class2', 'class_short', 'donor', 'facs_sorting', 'facs_gate', 'n_counts', 'n_features', 'mito_counts', 'mito_pct', 'drug', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19', 'PC20', 'PC21', 'PC22', 'PC23', 'PC24', 'UMAP1', 'UMAP2', 'leiden0.5', 'leiden0.644', 'leiden0.789', 'leiden0.933', 'leiden1.08', 'leiden1.22', 'leiden1.37', 'leiden1.51', 'leiden1.66', 'leiden1.8', 'leiden', 'barcode', 'TRAV', 'TRAD', 'TRAJ', 'TRAC', 'TRA_cdr3', 'TRA_cdr3_nt', 'TRA_length', 'TRA_reads', 'TRA_umis', 'TRBV', 'TRBD', 'TRBJ', 'TRBC', 'TRB_cdr3', 'TRB_cdr3_nt', 'TRB_length', 'TRB_reads', 'TRB_umis', 'TRB_cdr3_trim', 'TRA_cdr3_trim', 'TRAV2', 'TRBV2', 'TRAV_cdr1', 'TRAV_cdr2', 'TRAV_cdr25', 'TRBV_cdr1', 'TRBV_cdr2', 'TRBV_cdr25', 'has_tcr', 'IGHV', 'IGHD', 'IGHJ', 'IGHC', 'IGH_cdr3', 'IGH_cdr3_nt', 'IGH_length', 'IGH_reads', 'IGH_umis', 'IGKV', 'IGKD', 'IGKJ', 'IGKC', 'IGK_cdr3', 'IGK_cdr3_nt', 'IGK_length', 'IGK_reads', 'IGK_umis', 'IGLV', 'IGLD', 'IGLJ', 'IGLC', 'IGL_cdr3', 'IGL_cdr3_nt', 'IGL_length', 'IGL_reads', 'IGL_umis', 'has_bcr', 'cluster'
    ##     var: 'mean', 'sd', 'percent', 'gene', 'exclude', 'include', 'fitted', 'residuals', 'rank'
    ##     uns: 'de', 'knn', 'mcv', 'pca', 'pca_h'

Read the data in Python:

    import anndata
    a1 = anndata.read_h5ad("tissue-cd8.h5ad")
    a1

    ## AnnData object with n_obs × n_vars = 25341 × 28165
    ##     obs: 'cell', 'channel', 'case', 'class', 'class2', 'class_short', 'donor', 'facs_sorting', 'facs_gate', 'n_counts', 'n_features', 'mito_counts', 'mito_pct', 'drug', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19', 'PC20', 'PC21', 'PC22', 'PC23', 'PC24', 'UMAP1', 'UMAP2', 'leiden0.5', 'leiden0.644', 'leiden0.789', 'leiden0.933', 'leiden1.08', 'leiden1.22', 'leiden1.37', 'leiden1.51', 'leiden1.66', 'leiden1.8', 'leiden', 'barcode', 'TRAV', 'TRAD', 'TRAJ', 'TRAC', 'TRA_cdr3', 'TRA_cdr3_nt', 'TRA_length', 'TRA_reads', 'TRA_umis', 'TRBV', 'TRBD', 'TRBJ', 'TRBC', 'TRB_cdr3', 'TRB_cdr3_nt', 'TRB_length', 'TRB_reads', 'TRB_umis', 'TRB_cdr3_trim', 'TRA_cdr3_trim', 'TRAV2', 'TRBV2', 'TRAV_cdr1', 'TRAV_cdr2', 'TRAV_cdr25', 'TRBV_cdr1', 'TRBV_cdr2', 'TRBV_cdr25', 'has_tcr', 'IGHV', 'IGHD', 'IGHJ', 'IGHC', 'IGH_cdr3', 'IGH_cdr3_nt', 'IGH_length', 'IGH_reads', 'IGH_umis', 'IGKV', 'IGKD', 'IGKJ', 'IGKC', 'IGK_cdr3', 'IGK_cdr3_nt', 'IGK_length', 'IGK_reads', 'IGK_umis', 'IGLV', 'IGLD', 'IGLJ', 'IGLC', 'IGL_cdr3', 'IGL_cdr3_nt', 'IGL_length', 'IGL_reads', 'IGL_umis', 'has_bcr', 'cluster'
    ##     var: 'mean', 'sd', 'percent', 'gene', 'exclude', 'include', 'fitted', 'residuals', 'rank'
    ##     uns: 'de', 'knn', 'mcv', 'pca', 'pca_h'

Use the file \[ensembl\_id-symbol.tsv\] to convert between gene
identifiers in this project:

    fread("ensembl_id-symbol.tsv")

    ##             ensembl_id       symbol
    ##     1: ENSG00000121410         A1BG
    ##     2: ENSG00000268895     A1BG-AS1
    ##     3: ENSG00000148584         A1CF
    ##     4: ENSG00000175899          A2M
    ##     5: ENSG00000245105      A2M-AS1
    ##    ---                             
    ## 33534: ENSG00000203995       ZYG11A
    ## 33535: ENSG00000162378       ZYG11B
    ## 33536: ENSG00000159840          ZYX
    ## 33537: ENSG00000074755        ZZEF1
    ## 33538: ENSG00000272920 hsa-mir-1253

Pseudobulk expression data
==========================

We provide `.h5` files with log2CPM values in `matrix` with one row per
gene and one column per sample. And a corresponding table in `obs` that
describes each sample.

Pseudobulk at the level of cell lineage
---------------------------------------

The expression data in log2CPM:

    h5 <- rhdf5::h5read("pseudobulk_donor.h5", "matrix")
    log2cpm <- Matrix::sparseMatrix(
      dims   = h5$shape,
      i      = as.numeric(h5$indices),
      p      = as.numeric(h5$indptr),
      x      = as.numeric(h5$data),
      index1 = FALSE
    )
    rownames(log2cpm) <- h5$features
    colnames(log2cpm) <- h5$barcodes
    dim(log2cpm)

    ## [1] 33538   249

    log2cpm[1:5,1:5]

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##                 Tissue CD8 T cells|MC_1 Tissue CD8 T cells|MC_2 Tissue CD8 T cells|MC_9
    ## ENSG00000243485                    .                      .                        .   
    ## ENSG00000238009                    .                      0.546                    .   
    ## ENSG00000241599                    .                      .                        .   
    ## ENSG00000235146                    .                      .                        .   
    ## ENSG00000237491                    5.36                   5.870                    5.13
    ##                 Tissue CD8 T cells|SIC_100 Tissue CD8 T cells|SIC_109
    ## ENSG00000243485                      .                           .   
    ## ENSG00000238009                      0.549                       1.42
    ## ENSG00000241599                      .                           .   
    ## ENSG00000235146                      .                           .   
    ## ENSG00000237491                      5.401                       5.24

The table with one row per sample (notice `id` is in the format
`{analysis}|{donor}`):

    obs <- h5read("pseudobulk_donor.h5", "obs")
    str(obs)

    ## 'data.frame':    249 obs. of  7 variables:
    ##  $ donor  : chr [1:249(1d)] "MC_1" "MC_2" "MC_9" "SIC_100" ...
    ##  $ case   : chr [1:249(1d)] "Control" "Control" "Control" "Case" ...
    ##  $ drug   : chr [1:249(1d)] "None" "None" "None" "PD-1" ...
    ##  $ sex    : chr [1:249(1d)] "F" "M" "F" "F" ...
    ##  $ dataset: chr [1:249(1d)] "Tissue CD8 T cells" "Tissue CD8 T cells" "Tissue CD8 T cells" "Tissue CD8 T cells" ...
    ##  $ tissue : chr [1:249(1d)] "Tissue" "Tissue" "Tissue" "Tissue" ...
    ##  $ id     : chr [1:249(1d)] "Tissue CD8 T cells|MC_1" "Tissue CD8 T cells|MC_2" "Tissue CD8 T cells|MC_9" "Tissue CD8 T cells|SIC_100" ...

Complete contents:

    h5ls("pseudobulk_donor.h5")

    ##     group     name       otype   dclass     dim
    ## 0       /   matrix   H5I_GROUP                 
    ## 1 /matrix barcodes H5I_DATASET   STRING     249
    ## 2 /matrix     data H5I_DATASET    FLOAT 3567136
    ## 3 /matrix features H5I_DATASET   STRING   33538
    ## 4 /matrix  indices H5I_DATASET  INTEGER 3567136
    ## 5 /matrix   indptr H5I_DATASET  INTEGER     250
    ## 6 /matrix    shape H5I_DATASET  INTEGER       2
    ## 7       /      obs H5I_DATASET COMPOUND     249

Pseudobulk at the level of cell cluster
---------------------------------------

The expression data in log2CPM:

    h5 <- rhdf5::h5read("pseudobulk_donor_cluster.h5", "matrix")
    log2cpm <- Matrix::sparseMatrix(
      dims   = h5$shape,
      i      = as.numeric(h5$indices),
      p      = as.numeric(h5$indptr),
      x      = as.numeric(h5$data),
      index1 = FALSE
    )
    rownames(log2cpm) <- h5$features
    colnames(log2cpm) <- h5$barcodes
    dim(log2cpm)

    ## [1] 33538  2958

    log2cpm[1:5,1:5]

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##                 Tissue CD8 T cells|1|MC_1 Tissue CD8 T cells|2|MC_1 Tissue CD8 T cells|3|MC_1
    ## ENSG00000243485                      .                         .                            .
    ## ENSG00000238009                      .                         .                            .
    ## ENSG00000241599                      .                         .                            .
    ## ENSG00000235146                      .                         .                            .
    ## ENSG00000237491                      1.25                      1.64                         .
    ##                 Tissue CD8 T cells|4|MC_1 Tissue CD8 T cells|5|MC_1
    ## ENSG00000243485                      .                         .   
    ## ENSG00000238009                      .                         .   
    ## ENSG00000241599                      .                         .   
    ## ENSG00000235146                      .                         .   
    ## ENSG00000237491                      2.35                      1.57

The table with one row per sample (notice `id` is in the format
`{analysis}|{cluster}|{donor}`):

    obs <- h5read("pseudobulk_donor_cluster.h5", "obs")
    str(obs)

    ## 'data.frame':    2958 obs. of  8 variables:
    ##  $ cluster: chr [1:2958(1d)] "1" "2" "3" "4" ...
    ##  $ donor  : chr [1:2958(1d)] "MC_1" "MC_1" "MC_1" "MC_1" ...
    ##  $ case   : chr [1:2958(1d)] "Control" "Control" "Control" "Control" ...
    ##  $ drug   : chr [1:2958(1d)] "None" "None" "None" "None" ...
    ##  $ sex    : chr [1:2958(1d)] "F" "F" "F" "F" ...
    ##  $ dataset: chr [1:2958(1d)] "Tissue CD8 T cells" "Tissue CD8 T cells" "Tissue CD8 T cells" "Tissue CD8 T cells" ...
    ##  $ tissue : chr [1:2958(1d)] "Tissue" "Tissue" "Tissue" "Tissue" ...
    ##  $ id     : chr [1:2958(1d)] "Tissue CD8 T cells|1|MC_1" "Tissue CD8 T cells|2|MC_1" "Tissue CD8 T cells|3|MC_1" "Tissue CD8 T cells|4|MC_1" ...

Complete contents:

    h5ls("pseudobulk_donor_cluster.h5")

    ##     group     name       otype   dclass      dim
    ## 0       /   matrix   H5I_GROUP                  
    ## 1 /matrix barcodes H5I_DATASET   STRING     2958
    ## 2 /matrix     data H5I_DATASET    FLOAT 24693797
    ## 3 /matrix features H5I_DATASET   STRING    33538
    ## 4 /matrix  indices H5I_DATASET  INTEGER 24693797
    ## 5 /matrix   indptr H5I_DATASET  INTEGER     2959
    ## 6 /matrix    shape H5I_DATASET  INTEGER        2
    ## 7       /      obs H5I_DATASET COMPOUND     2958

Analysis results
================

Cell cluster differential abundance statistics for Case vs Control
------------------------------------------------------------------

Which cell clusters have greater relative abundance in irColitis Cases
compared to Controls?

Columns:

-   `analysis`: The dataset being tested.
-   `cluster`: The cluster being tested.
-   `value`: The name of the coefficient in the model.
-   `lrt_p`: The likelihoood ratio test p-value comparing two models.
-   `lrt_fdr`: The false discovery rate (adjusted p-value) for the
    likelihood ratio test.
-   `est`: Estimate of the coefficient.
-   `sd`: Coefficient standard deviation.
-   `z`: Coefficient Z-score.
-   `p`: Coefficient P-value.
-   `est_low`: Lower bound of the 95% confidence interval of the
    coefficient.
-   `est_high`: Upper bound of the 95% confidence interval of the
    coefficient.
-   `OR`: The logisitic regression odds ratio.
-   `OR_2.5`: The left side of the 95% confidence interval.
-   `OR_97.5`: The right side of the 95% confidence interval.

The likelihood ratio test compares two models:

    Model 0: in_cluster ~ 1 + sex + (1|donor)
    Model 1: in_cluster ~ 1 + sex + case + (1|donor)

When fitting this model, the variable `case` has two levels `Case` and
`Control`.

    d <- fread("cluster-abundance_case-vs-control.tsv")
    d

    ##                      analysis cluster       value    lrt_p  lrt_fdr     est     sd      z         p
    ##   1:       Tissue CD8 T cells       1 (Intercept) 1.71e-03 5.12e-03 -0.9113 0.1464  -6.23  4.79e-10
    ##   2:       Tissue CD8 T cells       1        sexM 1.71e-03 5.12e-03 -0.0698 0.1623  -0.43  6.67e-01
    ##   3:       Tissue CD8 T cells       1        Case 1.71e-03 5.12e-03 -0.5611 0.1615  -3.47  5.14e-04
    ##   4:       Tissue CD8 T cells       2 (Intercept) 2.66e-01 2.91e-01 -1.8476 0.1606 -11.50  1.30e-30
    ##   5:       Tissue CD8 T cells       2        sexM 2.66e-01 2.91e-01  0.2199 0.1763   1.25  2.12e-01
    ##  ---                                                                                               
    ## 359: Luoma Tissue CD4 T cells       9        Case 1.22e-01 1.49e-01  0.2348 0.1464   1.60  1.09e-01
    ## 360: Luoma Tissue CD4 T cells      10 (Intercept) 1.31e-02 1.80e-02 -3.4681 0.0985 -35.19 2.68e-271
    ## 361: Luoma Tissue CD4 T cells      10        Case 1.31e-02 1.80e-02  0.4143 0.1484   2.79  5.25e-03
    ## 362: Luoma Tissue CD4 T cells      11 (Intercept) 3.77e-06 4.15e-05 -5.5139 0.1965 -28.06 3.11e-173
    ## 363: Luoma Tissue CD4 T cells      11        Case 3.77e-06 4.15e-05  1.2603 0.2345   5.38  7.64e-08
    ##      est_low est_high      OR  OR_2.5 OR_97.5
    ##   1: -1.1982   -0.624 0.40200 0.30173 0.53557
    ##   2: -0.3879    0.248 0.93256 0.67849 1.28176
    ##   3: -0.8777   -0.244 0.57059 0.41575 0.78311
    ##   4: -2.1624   -1.533 0.15762 0.11505 0.21595
    ##   5: -0.1257    0.565 1.24591 0.88186 1.76024
    ##  ---                                         
    ## 359: -0.0521    0.522 1.26471 0.94923 1.68504
    ## 360: -3.6612   -3.275 0.03118 0.02570 0.03782
    ## 361:  0.1234    0.705 1.51338 1.13136 2.02439
    ## 362: -5.8990   -5.129 0.00403 0.00274 0.00592
    ## 363:  0.8008    1.720 3.52663 2.22731 5.58392

    # Which Tissue CD8 T cell clusters are differentially abundant?
    d %>%
      filter(analysis == "Tissue CD8 T cells", value == "Case") %>%
      arrange(lrt_p)

    ##               analysis cluster value    lrt_p  lrt_fdr    est    sd      z        p est_low
    ##  1: Tissue CD8 T cells       4  Case 4.24e-05 0.000509  0.832 0.174  4.776 1.79e-06  0.4905
    ##  2: Tissue CD8 T cells       7  Case 4.67e-04 0.002801  0.557 0.133  4.186 2.85e-05  0.2963
    ##  3: Tissue CD8 T cells       3  Case 7.52e-04 0.003007  0.556 0.152  3.656 2.56e-04  0.2581
    ##  4: Tissue CD8 T cells       1  Case 1.71e-03 0.005115 -0.561 0.162 -3.473 5.14e-04 -0.8777
    ##  5: Tissue CD8 T cells      12  Case 1.06e-02 0.025469  0.843 0.308  2.737 6.20e-03  0.2393
    ##  6: Tissue CD8 T cells      11  Case 3.44e-02 0.063479  0.503 0.229  2.199 2.79e-02  0.0546
    ##  7: Tissue CD8 T cells      10  Case 3.70e-02 0.063479  0.519 0.231  2.248 2.46e-02  0.0666
    ##  8: Tissue CD8 T cells       6  Case 5.06e-02 0.075933 -0.278 0.135 -2.055 3.99e-02 -0.5426
    ##  9: Tissue CD8 T cells       9  Case 7.26e-02 0.096814 -0.253 0.139 -1.820 6.88e-02 -0.5263
    ## 10: Tissue CD8 T cells       8  Case 1.41e-01 0.169070 -0.500 0.331 -1.508 1.32e-01 -1.1489
    ## 11: Tissue CD8 T cells       2  Case 2.66e-01 0.290676  0.197 0.175  1.124 2.61e-01 -0.1466
    ## 12: Tissue CD8 T cells       5  Case 4.22e-01 0.421721 -0.146 0.179 -0.814 4.16e-01 -0.4974
    ##     est_high    OR OR_2.5 OR_97.5
    ##  1:   1.1733 2.298  1.633   3.233
    ##  2:   0.8181 1.746  1.345   2.266
    ##  3:   0.8546 1.744  1.294   2.350
    ##  4:  -0.2445 0.571  0.416   0.783
    ##  5:   1.4465 2.323  1.270   4.248
    ##  6:   0.9509 1.653  1.056   2.588
    ##  7:   0.9718 1.681  1.069   2.643
    ##  8:  -0.0129 0.758  0.581   0.987
    ##  9:   0.0195 0.776  0.591   1.020
    ## 10:   0.1498 0.607  0.317   1.162
    ## 11:   0.5412 1.218  0.864   1.718
    ## 12:   0.2055 0.864  0.608   1.228

Cell cluster differential abundance statistics for anti-PD-1/anti-CTLA-4 vs anti-PD-1 within Cases, not including Controls
--------------------------------------------------------------------------------------------------------------------------

Which cell clusters are associated with dual vs single therapy within
cases? (Controls are excluded from this analysis.)

Columns:

-   `analysis`: The dataset being tested.
-   `cluster`: The cluster being tested.
-   `value`: The name of the coefficient in the model.
-   `lrt_p`: The likelihoood ratio test p-value comparing two models.
-   `lrt_fdr`: The false discovery rate (adjusted p-value) for the
    likelihood ratio test.
-   `est`: Estimate of the coefficient.
-   `sd`: Coefficient standard deviation.
-   `z`: Coefficient Z-score.
-   `p`: Coefficient P-value.
-   `est_low`: Lower bound of the 95% confidence interval of the
    coefficient.
-   `est_high`: Upper bound of the 95% confidence interval of the
    coefficient.
-   `OR`: The logisitic regression odds ratio.
-   `OR_2.5`: The left side of the 95% confidence interval.
-   `OR_97.5`: The right side of the 95% confidence interval.

The likelihood ratio test compares two models:

    Model 0: in_cluster ~ 1 + sex + (1|donor)
    Model 1: in_cluster ~ 1 + sex + drug + (1|donor)

The variable `drug` has two levels `anti-PD-1/anti-CTLA-4` and
`anti-PD-1`.

    fread("cluster-abundance_drug-within-cases.tsv")

    ##                analysis cluster           value  lrt_p lrt_fdr    est    sd       z         p
    ##   1: Tissue CD8 T cells       1     (Intercept) 0.1804   0.309 -1.441 0.149  -9.660  4.45e-22
    ##   2: Tissue CD8 T cells       1 drugPD-1/CTLA-4 0.1804   0.309 -0.331 0.235  -1.409  1.59e-01
    ##   3: Tissue CD8 T cells       1            sexM 0.1804   0.309  0.133 0.257   0.519  6.04e-01
    ##   4: Tissue CD8 T cells      10     (Intercept) 0.0521   0.208 -4.147 0.155 -26.799 3.36e-158
    ##   5: Tissue CD8 T cells      10 drugPD-1/CTLA-4 0.0521   0.208  0.526 0.229   2.298  2.16e-02
    ##  ---                                                                                         
    ## 311:  Blood CD8 T cells       8 drugPD-1/CTLA-4 0.3179   0.445  0.481 0.464   1.035  3.01e-01
    ## 312:  Blood CD8 T cells       8            sexM 0.3179   0.445  0.602 0.447   1.347  1.78e-01
    ## 313:  Blood CD8 T cells       9     (Intercept) 0.1762   0.335 -3.091 0.262 -11.793  4.23e-32
    ## 314:  Blood CD8 T cells       9 drugPD-1/CTLA-4 0.1762   0.335 -0.680 0.512  -1.328  1.84e-01
    ## 315:  Blood CD8 T cells       9            sexM 0.1762   0.335  0.849 0.464   1.829  6.74e-02
    ##      est_low est_high
    ##   1: -1.7334   -1.149
    ##   2: -0.7919    0.130
    ##   3: -0.3697    0.636
    ##   4: -4.4508   -3.844
    ##   5:  0.0774    0.975
    ##  ---                 
    ## 311: -0.4296    1.391
    ## 312: -0.2739    1.478
    ## 313: -3.6053   -2.578
    ## 314: -1.6840    0.324
    ## 315: -0.0607    1.759

Cell cluster markers
--------------------

Which genes are good markers for each cell cluster?

These results were computed using the limma R package.

Columns:

-   `analysis`: The dataset being tested.
-   `contrast`: The cluster being tested (one versus all).
-   `gene`: The gene symbol.
-   `auc`: The area under the receiver operator curve for this gene
    predicting cell membership to this cluster.
-   `ensembl_id`: The Ensembl gene ID.
-   `log_fc`: The log2 fold-change.
-   `ci_l`: The left side of the 95% confidence interval.
-   `ci_r`: The right side of the 95% confidence interval.
-   `ave_expr`: The mean expression of the gene in log2CPM.
-   `t`: The t-statistic.
-   `p_value`: The limma p-value.
-   `adj_p_value`: The false discovery rate (adjusted p-value).
-   `b`: The log odds ratio.

<!-- -->

    fread("de-ova.tsv.gz")

    ##                    analysis contrast       gene   auc      ensembl_id    log_fc    ci_l   ci_r
    ##       1: Tissue CD8 T cells 1 vs all      ODAPH 0.502 ENSG00000174792  1.63e-01  0.1240 0.2020
    ##       2: Tissue CD8 T cells 1 vs all AC068987.3 0.500 ENSG00000260473  5.31e-02  0.0396 0.0666
    ##       3: Tissue CD8 T cells 1 vs all AC120193.1 0.510 ENSG00000253535  4.79e-01  0.3530 0.6060
    ##       4: Tissue CD8 T cells 1 vs all      HPSE2 0.501 ENSG00000172987  7.29e-02  0.0516 0.0942
    ##       5: Tissue CD8 T cells 1 vs all     OR52I1 0.500 ENSG00000232268  9.95e-02  0.0683 0.1310
    ##      ---                                                                                      
    ## 2459151:  Blood CD8 T cells 9 vs all     MBTPS2 0.497 ENSG00000012174 -6.66e-05 -0.2940 0.2940
    ## 2459152:  Blood CD8 T cells 9 vs all  LINC00685 0.499 ENSG00000226179  4.13e-05 -0.2160 0.2160
    ## 2459153:  Blood CD8 T cells 9 vs all     NLGN4Y 0.500 ENSG00000165246  2.33e-05 -0.1440 0.1440
    ## 2459154:  Blood CD8 T cells 9 vs all AL513320.1 0.500 ENSG00000238260 -1.54e-05 -0.1200 0.1200
    ## 2459155:  Blood CD8 T cells 9 vs all    CCDC189 0.499 ENSG00000196118  1.63e-05 -0.1540 0.1540
    ##          ave_expr         t  p_value adj_p_val     b
    ##       1:  0.03080  8.170000 7.00e-15  1.61e-10 23.20
    ##       2:  0.00442  7.740000 1.33e-13  1.53e-09 20.40
    ##       3:  0.14700  7.440000 9.23e-13  7.09e-09 18.50
    ##       4:  0.01090  6.740000 7.28e-11  4.19e-07 14.30
    ##       5:  0.01630  6.280000 1.08e-09  4.99e-06 11.70
    ##      ---                                            
    ## 2459151:  0.70800 -0.000446 1.00e+00  1.00e+00 -6.41
    ## 2459152:  0.35600  0.000376 1.00e+00  1.00e+00 -6.41
    ## 2459153:  0.13200  0.000318 1.00e+00  1.00e+00 -6.41
    ## 2459154:  0.10600 -0.000253 1.00e+00  1.00e+00 -6.41
    ## 2459155:  0.21700  0.000209 1.00e+00  1.00e+00 -6.41

Differential gene expression statistics for three contrasts
-----------------------------------------------------------

Which genes are differentially expressed for each of these three
contrasts?

-   Case vs Control
-   CasePD1CTLA4 vs CasePD1
-   ControlPD1 vs ControlNone

These results were computed using the limma R package.

Columns:

-   `analysis`: The dataset being tested.
-   `cluster`: The name of the cell cluster. All cells if it ends with
    `-all`.
-   `contrast`: One of the three contrasts (Case vs Control,
    CasePD1CTLA4 vs CasePD1, ControlPD1 vs ControlNone).
-   `gene`: The gene symbol.
-   `ensembl_id`: The Ensembl gene ID.
-   `log_fc`: The log2 fold-change.
-   `ci_l`: The left side of the 95% confidence interval.
-   `ci_r`: The right side of the 95% confidence interval.
-   `ave_expr`: The mean expression of the gene in log2CPM.
-   `t`: The t-statistic.
-   `p_value`: The limma p-value.
-   `adj_p_value`: The false discovery rate (adjusted p-value).
-   `b`: The log odds ratio.

<!-- -->

    fread("de-contrasts.tsv.gz")

    ##                    analysis   cluster       gene percent    log_fc   ci_l  ci_r ave_expr         t
    ##       1: Tissue CD8 T cells T-CD8-all       IL26   3.910  4.69e+00  3.490 5.880    2.690  8.06e+00
    ##       2: Tissue CD8 T cells T-CD8-all     DEPDC1   2.030  2.89e+00  2.080 3.700    3.400  7.29e+00
    ##       3: Tissue CD8 T cells T-CD8-all      IL17A   1.030  4.85e+00  3.470 6.230    2.670  7.22e+00
    ##       4: Tissue CD8 T cells T-CD8-all AC010424.3   0.564  2.86e+00  2.040 3.680    1.260  7.14e+00
    ##       5: Tissue CD8 T cells T-CD8-all     DLGAP5   3.220  3.54e+00  2.520 4.570    4.310  7.09e+00
    ##      ---                                                                                          
    ## 3805056:  Blood CD8 T cells   B-CD8-9       RPS3      NA  6.57e-04 -0.905 0.907    7.440  1.52e-03
    ## 3805057:  Blood CD8 T cells   B-CD8-9     ANP32B      NA  1.01e-03 -1.500 1.500    3.650  1.41e-03
    ## 3805058:  Blood CD8 T cells   B-CD8-9     ZNF365      NA -6.52e-04 -1.190 1.180    0.761 -1.15e-03
    ## 3805059:  Blood CD8 T cells   B-CD8-9     RPL18A      NA  2.69e-04 -0.901 0.901    7.070  6.24e-04
    ## 3805060:  Blood CD8 T cells   B-CD8-9     HIGD1A      NA  1.18e-05 -1.480 1.480    1.130  1.67e-05
    ##           p_value adj_p_val     b      ensembl_id                contrast
    ##       1: 1.25e-08  0.000184  9.70 ENSG00000111536         Case vs Control
    ##       2: 8.11e-08  0.000316  7.98 ENSG00000024526         Case vs Control
    ##       3: 9.76e-08  0.000316  7.81 ENSG00000112115         Case vs Control
    ##       4: 1.19e-07  0.000316  7.63 ENSG00000282925         Case vs Control
    ##       5: 1.35e-07  0.000316  7.51 ENSG00000126787         Case vs Control
    ##      ---                                                                 
    ## 3805056: 9.99e-01  0.999000 -5.89 ENSG00000149273 CasePD1CTLA4 vs CasePD1
    ## 3805057: 9.99e-01  0.999000 -5.89 ENSG00000136938 CasePD1CTLA4 vs CasePD1
    ## 3805058: 9.99e-01  0.999000 -5.89 ENSG00000138311 CasePD1CTLA4 vs CasePD1
    ## 3805059: 1.00e+00  1.000000 -5.89 ENSG00000105640 CasePD1CTLA4 vs CasePD1
    ## 3805060: 1.00e+00  1.000000 -5.89 ENSG00000181061 CasePD1CTLA4 vs CasePD1

Cell-cell communication analysis with means of gene pairs
---------------------------------------------------------

This table has Case vs Control statistics computed by limma for the
means of 1610 gene pairs, for each pair of cell types.

Columns for cell type 1:

-   `c1`: The cell type being tested for gene 1.
-   `e1`: The Ensembl gene ID of gene 1.
-   `donors1`: The number of donors with expression of gene 1.
-   `g1`: The name of gene 1.

Columns for cell type 2:

-   `c2`: The cell type being tested for gene 2.
-   `e2`: The Ensembl gene ID of gene 2.
-   `donors2`: The number of donors with expression of gene 2.
-   `g2`: The name of gene 2.

Columns from limma linear model results comparing cases versus controls:

-   `logFC`: The log2 fold-change.
-   `CI.L`: The left side of the 95% confidence interval.
-   `CI.R`: The right side of the 95% confidence interval.
-   `AveExpr`: The mean expression of the gene in log2CPM.
-   `t`: The t-statistic.
-   `P.Value`: The limma p-value.
-   `adj.P.Val`: The false discovery rate (adjusted p-value).
-   `B`: The log odds ratio.

Other columns:

-   `id`: A unique identifier in the format `{c2} {g2} {c1} {g1}`.
-   `percent1`: The percent of cells with expression of gene 1.
-   `percent2`: The percent of cells with expression of gene 2.
-   `fdr`: The false discovery rate (adjusted p-value) for gene pairs
    within each cluster pair.
-   `epi`: Does the interaction involve epithelial nuclei?
-   `inter`: Does the interaction involve different cell types?

<!-- -->

    fread("ccc-limma-cell-lineages.tsv.gz")

    ##         c1              e1 donors1        g1  c2              e2 donors2       g2     logFC   CI.L
    ##     1:   M ENSG00000243646      20    IL10RB CD8 ENSG00000111536      17     IL26  3.614530  2.825
    ##     2: CD8 ENSG00000184451      21     CCR10 CD8 ENSG00000156234      14   CXCL13  3.011380  2.310
    ##     3:   E ENSG00000243646      26    IL10RB CD8 ENSG00000111536      17     IL26  2.403400  1.794
    ##     4:   E ENSG00000116329      25     OPRD1 CD8 ENSG00000156234      14   CXCL13  2.930270  2.216
    ##     5:   M ENSG00000163702      21    IL17RC CD8 ENSG00000112115      15    IL17A  3.565190  2.685
    ##    ---                                                                                            
    ## 13235:   B ENSG00000138448      27     ITGAV   B ENSG00000108821      12   COL1A1  0.000380 -0.437
    ## 13236:   B ENSG00000240505      27 TNFRSF13B   B ENSG00000102524      25 TNFSF13B -0.000556 -1.030
    ## 13237:   B ENSG00000124145      11      SDC4   E ENSG00000137801      26    THBS1 -0.000494 -0.931
    ## 13238:   B ENSG00000155760      17      FZD7   E ENSG00000085741      26    WNT11 -0.000162 -0.409
    ## 13239:   B ENSG00000106799      27    TGFBR1 CD8 ENSG00000164404      15     GDF9 -0.000127 -0.574
    ##         CI.R AveExpr         t  P.Value adj.P.Val     B                     id percent1 percent2
    ##     1: 4.404    2.71  9.274530 3.64e-11  9.09e-10 15.37      CD8 IL26 M IL10RB  0.04683 0.039067
    ##     2: 3.713    1.90  8.747690 5.50e-10  1.38e-08 12.77   CD8 CXCL13 CD8 CCR10  0.00312 0.016456
    ##     3: 3.013    2.44  7.988070 1.49e-09  1.87e-08 11.68      CD8 IL26 E IL10RB  0.00974 0.039067
    ##     4: 3.645    1.84  8.360580 1.88e-09  4.69e-08 11.56     CD8 CXCL13 E OPRD1  0.00286 0.016456
    ##     5: 4.446    2.84  8.260210 2.65e-09  4.64e-08 11.23     CD8 IL17A M IL17RC  0.04951 0.010299
    ##    ---                                                                                          
    ## 13235: 0.438    1.31  0.001770 9.99e-01  9.99e-01 -6.04       B COL1A1 B ITGAV  0.01279 0.000768
    ## 13236: 1.029    4.27 -0.001093 9.99e-01  9.99e-01 -6.17 B TNFSF13B B TNFRSF13B  0.23912 0.011871
    ## 13237: 0.930    3.48 -0.001082 9.99e-01  9.99e-01 -6.65         E THBS1 B SDC4  0.00052 0.173278
    ## 13238: 0.409    1.33 -0.000810 9.99e-01  1.00e+00 -5.56         E WNT11 B FZD7  0.00245 0.006829
    ## 13239: 0.574    1.70 -0.000444 1.00e+00  1.00e+00 -4.64      CD8 GDF9 B TGFBR1  0.02163 0.001026
    ##             fdr   epi inter
    ##     1: 1.27e-08 FALSE  TRUE
    ##     2: 2.40e-07 FALSE FALSE
    ##     3: 5.64e-07  TRUE  TRUE
    ##     4: 5.64e-07  TRUE  TRUE
    ##     5: 4.62e-07 FALSE  TRUE
    ##    ---                     
    ## 13235: 9.99e-01 FALSE FALSE
    ## 13236: 9.99e-01 FALSE FALSE
    ## 13237: 9.99e-01  TRUE  TRUE
    ## 13238: 9.99e-01  TRUE  TRUE
    ## 13239: 1.00e+00 FALSE  TRUE

Cell-cell communication analysis with Spearman correlation of gene percentages
------------------------------------------------------------------------------

This table has Spearman correlation coefficients for 1441 gene pairs,
for each pair of cell types.

These results were computed using the limma R package.

Columns for cell type 1:

-   `c1`: The cell type being tested for gene 1.
-   `g1`: The name of gene 1.
-   `donors1`: The number of donors with expression of gene 1.
-   `percent1`: The percent of cells with expression of gene 1.
-   `mean1`: The mean log2CPM of gene 1.

Columns for cell type 2:

-   `c2`: The cell type being tested for gene 2.
-   `g2`: The name of gene 2.
-   `donors2`: The number of donors with expression of gene 2.
-   `percent2`: The percent of cells with expression of gene 2.
-   `mean2`: The mean log2CPM of gene 2.

Other columns:

-   `estimate`: Spearman correlation.
-   `p.value`: Spearman correlation p-value.
-   `fdr`: The false discovery rate (adjusted p-value) for gene pairs
    within each cluster pair.
-   `epi`: Does the interaction involve epithelial nuclei?
-   `inter`: Does the interaction involve different cell types?

<!-- -->

    fread("ccc-spearman-cell-lineages.tsv.gz")

    ##         c1      g1 donors1  c2    g2 donors2 estimate  p.value percent1    mean1 percent2    mean2
    ##     1: CD4   PDIA3      26 CD8  CALR      27    0.866 1.10e-08 0.336965 0.364139 0.499941 0.655255
    ##     2: CD4   KCNA3      26 CD8  IL16      27    0.840 7.74e-08 0.074123 0.069246 0.332465 0.348667
    ##     3: CD8   ITGAL      27 CD4 ICAM3      26    0.832 1.33e-07 0.283375 0.290258 0.403503 0.444684
    ##     4: CD8    CD8A      27   B HLA-A      27    0.822 1.42e-07 0.647686 1.263167 0.834382 1.796378
    ##     5: CD4   ITGB2      26 CD8 ICAM3      27    0.809 5.65e-07 0.352410 0.398149 0.395249 0.424643
    ##    ---                                                                                            
    ## 10097: CD8   SCN8A      11   B FGF13      11   -1.000 1.00e+00 0.000829 0.000555 0.000421 0.000443
    ## 10098:   E   LRRC4      15 CD8 NTNG2      16    0.000 1.00e+00 0.000783 0.000628 0.001263 0.001155
    ## 10099:   B    CCR1      25   E CCL13      12    0.000 1.00e+00 0.011325 0.013007 0.000379 0.000415
    ## 10100:   E IL22RA2      13 CD4  IL22      12    1.000 1.00e+00 0.000281 0.000259 0.005171 0.006594
    ## 10101: CD8   IL2RG      27 CD4  IL21      23    0.000 1.00e+00 0.265893 0.281055 0.033235 0.037720
    ##             fdr   epi inter                  id
    ##     1: 4.37e-06 FALSE  TRUE  CD8 CALR CD4 PDIA3
    ##     2: 1.54e-05 FALSE  TRUE  CD8 IL16 CD4 KCNA3
    ##     3: 4.94e-05 FALSE  TRUE CD4 ICAM3 CD8 ITGAL
    ##     4: 7.64e-05 FALSE  TRUE    B HLA-A CD8 CD8A
    ##     5: 7.47e-05 FALSE  TRUE CD8 ICAM3 CD4 ITGB2
    ##    ---                                         
    ## 10097: 1.00e+00 FALSE  TRUE   B FGF13 CD8 SCN8A
    ## 10098: 1.00e+00  TRUE  TRUE   CD8 NTNG2 E LRRC4
    ## 10099: 1.00e+00  TRUE  TRUE      E CCL13 B CCR1
    ## 10100: 1.00e+00  TRUE  TRUE  CD4 IL22 E IL22RA2
    ## 10101: 1.00e+00 FALSE  TRUE  CD4 IL21 CD8 IL2RG

Correlation analysis of cluster abundances
------------------------------------------

This table has Spearman correlation coefficients for 2145 cell cluster
pairs.

Columns for cell type 1:

-   `c1`: The cell type being tested for gene 1.
-   `g1`: The name of gene 1.
-   `donors1`: The number of donors with expression of gene 1.
-   `percent1`: The percent of cells with expression of gene 1.
-   `mean1`: The mean log2CPM of gene 1.

Columns for cell type 2:

-   `c2`: The cell type being tested for gene 2.
-   `g2`: The name of gene 2.
-   `donors2`: The number of donors with expression of gene 2.
-   `percent2`: The percent of cells with expression of gene 2.
-   `mean2`: The mean log2CPM of gene 2.

Other columns:

-   `estimate`: Spearman correlation.
-   `p.value`: Spearman correlation p-value.
-   `fdr`: The false discovery rate (adjusted p-value) for gene pairs
    within each cluster pair.
-   `epi`: Does the interaction involve epithelial nuclei?
-   `inter`: Does the interaction involve different cell types?

<!-- -->

    fread("spearman-cluster-abundance.tsv.gz")

    ##         group cluster1 cluster2 n_donors estimate  p.value      fdr
    ##    1: Control      E-6     MP-6        4    1.000 0.00e+00 0.00e+00
    ##    2: Control     MP-6     MP-7        3    1.000 0.00e+00 0.00e+00
    ##    3:    Both      B-4      B-6       27    0.946 9.18e-14 1.97e-10
    ##    4:    Both      B-3      B-5       26    0.914 7.07e-11 7.58e-08
    ##    5:    Both      B-1      B-3       26   -0.901 3.56e-10 2.55e-07
    ##   ---                                                              
    ## 6431: Control      E-4    CD4-3       12    0.000 1.00e+00 1.00e+00
    ## 6432: Control     E-19     MP-6        4    0.000 1.00e+00 1.00e+00
    ## 6433: Control     E-20     MP-6        4    0.000 1.00e+00 1.00e+00
    ## 6434: Control     B-12     E-15       11    0.000 1.00e+00 1.00e+00
    ## 6435: Control     B-10   CD8-11       13    0.000 1.00e+00 1.00e+00
