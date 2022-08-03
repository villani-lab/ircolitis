We provide access to all of the data files that support the analyses in our publication:

-   Thomas MF, Slowikowski K, Manakongtreecheep K, Sen P, Tantivit J, Nasrallah M, et al. Altered interactions between circulating and tissue-resident CD8 T cells with the colonic mucosa define colitis associated with immune checkpoint inhibitors. bioRxiv. 2021. p. 2021.09.17.460868. <doi:10.1101/2021.09.17.460868>

Please contact [Kamil Slowikowski](mailto:kslowikowski@mgh.harvard.edu) with any questions.

## Analysis results

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
file
</th>
<th style="text-align:left;">
description
</th>
<th style="text-align:left;">
rows
</th>
<th style="text-align:right;">
cols
</th>
<th style="text-align:left;">
size
</th>
<th style="text-align:left;">
md5
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
cluster-abundance\_case-vs-control.tsv
</td>
<td style="text-align:left;">
Cell cluster differential abundance statistics for the model 'in\_cluster ~ 1 + sex + case + (1|donor)'
</td>
<td style="text-align:left;">
363
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
81.7 K
</td>
<td style="text-align:left;">
8a6bc704dc20bec9fb9730b59aa3ca7c
</td>
</tr>
<tr>
<td style="text-align:left;">
cluster-abundance\_drug-within-cases.tsv
</td>
<td style="text-align:left;">
Cell cluster differential abundance statistics for the model 'in\_cluster ~ 1 + sex + drug + (1|donor)' fit only to data from cases, not including controls
</td>
<td style="text-align:left;">
315
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
55.1 K
</td>
<td style="text-align:left;">
a8202fc2cab816f3bbae8975335febb9
</td>
</tr>
<tr>
<td style="text-align:left;">
de-ova.tsv.gz
</td>
<td style="text-align:left;">
Cell cluster differential gene expression statistics
</td>
<td style="text-align:left;">
2,459,155
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
53.1 M
</td>
<td style="text-align:left;">
0db04846f4e226543e4803b78751bb88
</td>
</tr>
<tr>
<td style="text-align:left;">
de-contrasts.tsv.gz
</td>
<td style="text-align:left;">
Differential gene expression statistics for three contrasts
</td>
<td style="text-align:left;">
3,805,060
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
82.0 M
</td>
<td style="text-align:left;">
5ed19f7c4ad9c3b8dbcbfa171a1599db
</td>
</tr>
<tr>
<td style="text-align:left;">
ccc-limma-cell-lineages.tsv.gz
</td>
<td style="text-align:left;">
Cell-cell communication analysis with means of gene pairs
</td>
<td style="text-align:left;">
13,239
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:left;">
792.3 K
</td>
<td style="text-align:left;">
5cfe86191df4c75fe734cb4157da0b8e
</td>
</tr>
<tr>
<td style="text-align:left;">
ccc-spearman-cell-lineages.tsv.gz
</td>
<td style="text-align:left;">
Cell-cell communication analysis with Spearman correlation
</td>
<td style="text-align:left;">
10,101
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:left;">
544.7 K
</td>
<td style="text-align:left;">
3c97d32ddca5b87ec29a10a352fd18fa
</td>
</tr>
<tr>
<td style="text-align:left;">
spearman-cluster-abundance.tsv.gz
</td>
<td style="text-align:left;">
Correlation analysis of cell cluster abundances
</td>
<td style="text-align:left;">
6,435
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
89.1 K
</td>
<td style="text-align:left;">
59974a12e5af282a8a9e3c4e843fbcfc
</td>
</tr>
</tbody>
</table>

## Pseudobulk expression data

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
file
</th>
<th style="text-align:left;">
description
</th>
<th style="text-align:left;">
samples
</th>
<th style="text-align:left;">
genes
</th>
<th style="text-align:left;">
size
</th>
<th style="text-align:left;">
md5
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
pseudobulk\_donor.h5
</td>
<td style="text-align:left;">
Pseudobulk expression and metadata at the level of cell lineages
</td>
<td style="text-align:left;">
249
</td>
<td style="text-align:left;">
33,538
</td>
<td style="text-align:left;">
21.3 M
</td>
<td style="text-align:left;">
871ebe6fd6f6e648ad7cdd3f44aa91a4
</td>
</tr>
<tr>
<td style="text-align:left;">
pseudobulk\_donor\_cluster.h5
</td>
<td style="text-align:left;">
Pseudobulk expression and metadata at the level of cell clusters
</td>
<td style="text-align:left;">
2,958
</td>
<td style="text-align:left;">
33,538
</td>
<td style="text-align:left;">
128.0 M
</td>
<td style="text-align:left;">
0b84b1b27c8879a3da960834aad27f3b
</td>
</tr>
</tbody>
</table>

## Single-cell expression data

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
file
</th>
<th style="text-align:left;">
description
</th>
<th style="text-align:left;">
cells
</th>
<th style="text-align:left;">
genes
</th>
<th style="text-align:left;">
size
</th>
<th style="text-align:left;">
md5
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
tissue-cd8.h5ad
</td>
<td style="text-align:left;">
Counts, metadata, and TCR sequences for Tissue CD8 T cells
</td>
<td style="text-align:left;">
25,341
</td>
<td style="text-align:left;">
28,165
</td>
<td style="text-align:left;">
189.4 M
</td>
<td style="text-align:left;">
b8fa7e186ec7880113d8302c5f85995e
</td>
</tr>
<tr>
<td style="text-align:left;">
tissue-cd4.h5ad
</td>
<td style="text-align:left;">
Counts, metadata, and TCR sequences for Tissue CD4 T cells
</td>
<td style="text-align:left;">
14,503
</td>
<td style="text-align:left;">
28,165
</td>
<td style="text-align:left;">
117.3 M
</td>
<td style="text-align:left;">
97d20672437de481e979e398cec67998
</td>
</tr>
<tr>
<td style="text-align:left;">
tissue-myeloid.h5ad
</td>
<td style="text-align:left;">
Counts and metadata for Tissue Myeloid cells
</td>
<td style="text-align:left;">
2,242
</td>
<td style="text-align:left;">
28,165
</td>
<td style="text-align:left;">
59.4 M
</td>
<td style="text-align:left;">
d1ed746170914c8dbde5758f90a6812e
</td>
</tr>
<tr>
<td style="text-align:left;">
tissue-b.h5ad
</td>
<td style="text-align:left;">
Counts, metadata, and BCR sequences for Tissue B cells
</td>
<td style="text-align:left;">
40,352
</td>
<td style="text-align:left;">
28,165
</td>
<td style="text-align:left;">
216.6 M
</td>
<td style="text-align:left;">
35bbbc0ccde0752bd610251011f1f503
</td>
</tr>
<tr>
<td style="text-align:left;">
tissue-epithelial.h5ad
</td>
<td style="text-align:left;">
Counts and metadata for Tissue epithelial and mesenchymal nuclei
</td>
<td style="text-align:left;">
81,707
</td>
<td style="text-align:left;">
33,538
</td>
<td style="text-align:left;">
665.8 M
</td>
<td style="text-align:left;">
dfdb05b933487629cca405aadfb9e299
</td>
</tr>
<tr>
<td style="text-align:left;">
blood-cd8.h5ad
</td>
<td style="text-align:left;">
Counts, metadata, and TCR sequences for Blood CD8 T cells
</td>
<td style="text-align:left;">
68,785
</td>
<td style="text-align:left;">
33,538
</td>
<td style="text-align:left;">
361.7 M
</td>
<td style="text-align:left;">
9512b1d4d04ec4e921b7922b839fc7c7
</td>
</tr>
<tr>
<td style="text-align:left;">
blood-cd4.h5ad
</td>
<td style="text-align:left;">
Counts, metadata, and TCR sequences for Blood CD4 T cells
</td>
<td style="text-align:left;">
13,733
</td>
<td style="text-align:left;">
33,538
</td>
<td style="text-align:left;">
123.6 M
</td>
<td style="text-align:left;">
445226683e3be4a83663a5316bfabb68
</td>
</tr>
<tr>
<td style="text-align:left;">
blood-myeloid.h5ad
</td>
<td style="text-align:left;">
Counts and metadata for Blood Myeloid cells
</td>
<td style="text-align:left;">
32,338
</td>
<td style="text-align:left;">
33,538
</td>
<td style="text-align:left;">
225.8 M
</td>
<td style="text-align:left;">
549811219d0151155143ad091a86c743
</td>
</tr>
<tr>
<td style="text-align:left;">
blood-b.h5ad
</td>
<td style="text-align:left;">
Counts and metadata for Blood B cells
</td>
<td style="text-align:left;">
3,962
</td>
<td style="text-align:left;">
33,538
</td>
<td style="text-align:left;">
82.5 M
</td>
<td style="text-align:left;">
b067b97a1fa82e14896ee8264115a7a0
</td>
</tr>
</tbody>
</table>

# Analysis results

## Cell cluster differential abundance statistics for Case vs Control

Columns:

-   `analysis`: The dataset being tested.
-   `cluster`: The cluster being tested.
-   `value`: The name of the coefficient in the model.
-   `lrt_p`: The likelihoood ratio test p-value comparing two models.
-   `lrt_fdr`: The false discovery rate (adjusted p-value) for the likelihood ratio test.
-   `est`: Estimate of the coefficient.
-   `sd`: Coefficient standard deviation.
-   `z`: Coefficient Z-score.
-   `p`: Coefficient P-value.
-   `est_low`: Lower bound of the 95% confidence interval of the coefficient.
-   `est_high`: Upper bound of the 95% confidence interval of the coefficient.
-   `OR`: The logisitic regression odds ratio.
-   `OR_2.5`: The left side of the 95% confidence interval.
-   `OR_97.5`: The right side of the 95% confidence interval.

The likelihood ratio test compares two models:

    Model 0: in_cluster ~ 1 + sex + (1|donor)
    Model 1: in_cluster ~ 1 + sex + case + (1|donor)

The variable `case` has two levels `Case` and `Control`.

``` r
fread("cluster-abundance_case-vs-control.tsv")
```

    ##                      analysis cluster       value    lrt_p  lrt_fdr     est
    ##   1:       Tissue CD8 T cells       1 (Intercept) 1.71e-03 5.12e-03 -0.9113
    ##   2:       Tissue CD8 T cells       1        sexM 1.71e-03 5.12e-03 -0.0698
    ##   3:       Tissue CD8 T cells       1        Case 1.71e-03 5.12e-03 -0.5611
    ##   4:       Tissue CD8 T cells       2 (Intercept) 2.66e-01 2.91e-01 -1.8476
    ##   5:       Tissue CD8 T cells       2        sexM 2.66e-01 2.91e-01  0.2199
    ##  ---                                                                       
    ## 359: Luoma Tissue CD4 T cells       9        Case 1.22e-01 1.49e-01  0.2348
    ## 360: Luoma Tissue CD4 T cells      10 (Intercept) 1.31e-02 1.80e-02 -3.4681
    ## 361: Luoma Tissue CD4 T cells      10        Case 1.31e-02 1.80e-02  0.4143
    ## 362: Luoma Tissue CD4 T cells      11 (Intercept) 3.77e-06 4.15e-05 -5.5139
    ## 363: Luoma Tissue CD4 T cells      11        Case 3.77e-06 4.15e-05  1.2603
    ##          sd      z         p est_low est_high      OR  OR_2.5 OR_97.5
    ##   1: 0.1464  -6.23  4.79e-10 -1.1982   -0.624 0.40200 0.30173 0.53557
    ##   2: 0.1623  -0.43  6.67e-01 -0.3879    0.248 0.93256 0.67849 1.28176
    ##   3: 0.1615  -3.47  5.14e-04 -0.8777   -0.244 0.57059 0.41575 0.78311
    ##   4: 0.1606 -11.50  1.30e-30 -2.1624   -1.533 0.15762 0.11505 0.21595
    ##   5: 0.1763   1.25  2.12e-01 -0.1257    0.565 1.24591 0.88186 1.76024
    ##  ---                                                                 
    ## 359: 0.1464   1.60  1.09e-01 -0.0521    0.522 1.26471 0.94923 1.68504
    ## 360: 0.0985 -35.19 2.68e-271 -3.6612   -3.275 0.03118 0.02570 0.03782
    ## 361: 0.1484   2.79  5.25e-03  0.1234    0.705 1.51338 1.13136 2.02439
    ## 362: 0.1965 -28.06 3.11e-173 -5.8990   -5.129 0.00403 0.00274 0.00592
    ## 363: 0.2345   5.38  7.64e-08  0.8008    1.720 3.52663 2.22731 5.58392

## Cell cluster differential abundance statistics for anti-PD-1/anti-CTLA-4 vs anti-PD-1 within Cases, not including Controls

Columns:

-   `analysis`: The dataset being tested.
-   `cluster`: The cluster being tested.
-   `value`: The name of the coefficient in the model.
-   `lrt_p`: The likelihoood ratio test p-value comparing two models.
-   `lrt_fdr`: The false discovery rate (adjusted p-value) for the likelihood ratio test.
-   `est`: Estimate of the coefficient.
-   `sd`: Coefficient standard deviation.
-   `z`: Coefficient Z-score.
-   `p`: Coefficient P-value.
-   `est_low`: Lower bound of the 95% confidence interval of the coefficient.
-   `est_high`: Upper bound of the 95% confidence interval of the coefficient.
-   `OR`: The logisitic regression odds ratio.
-   `OR_2.5`: The left side of the 95% confidence interval.
-   `OR_97.5`: The right side of the 95% confidence interval.

The likelihood ratio test compares two models:

    Model 0: in_cluster ~ 1 + sex + (1|donor)
    Model 1: in_cluster ~ 1 + sex + drug + (1|donor)

The variable `drug` has two levels `anti-PD-1/anti-CTLA-4` and `anti-PD-1`.

``` r
fread("cluster-abundance_drug-within-cases.tsv")
```

    ##                analysis cluster           value  lrt_p lrt_fdr    est    sd
    ##   1: Tissue CD8 T cells       1     (Intercept) 0.1804   0.309 -1.441 0.149
    ##   2: Tissue CD8 T cells       1 drugPD-1/CTLA-4 0.1804   0.309 -0.331 0.235
    ##   3: Tissue CD8 T cells       1            sexM 0.1804   0.309  0.133 0.257
    ##   4: Tissue CD8 T cells      10     (Intercept) 0.0521   0.208 -4.147 0.155
    ##   5: Tissue CD8 T cells      10 drugPD-1/CTLA-4 0.0521   0.208  0.526 0.229
    ##  ---                                                                       
    ## 311:  Blood CD8 T cells       8 drugPD-1/CTLA-4 0.3179   0.445  0.481 0.464
    ## 312:  Blood CD8 T cells       8            sexM 0.3179   0.445  0.602 0.447
    ## 313:  Blood CD8 T cells       9     (Intercept) 0.1762   0.335 -3.091 0.262
    ## 314:  Blood CD8 T cells       9 drugPD-1/CTLA-4 0.1762   0.335 -0.680 0.512
    ## 315:  Blood CD8 T cells       9            sexM 0.1762   0.335  0.849 0.464
    ##            z         p est_low est_high
    ##   1:  -9.660  4.45e-22 -1.7334   -1.149
    ##   2:  -1.409  1.59e-01 -0.7919    0.130
    ##   3:   0.519  6.04e-01 -0.3697    0.636
    ##   4: -26.799 3.36e-158 -4.4508   -3.844
    ##   5:   2.298  2.16e-02  0.0774    0.975
    ##  ---                                   
    ## 311:   1.035  3.01e-01 -0.4296    1.391
    ## 312:   1.347  1.78e-01 -0.2739    1.478
    ## 313: -11.793  4.23e-32 -3.6053   -2.578
    ## 314:  -1.328  1.84e-01 -1.6840    0.324
    ## 315:   1.829  6.74e-02 -0.0607    1.759

## Cell cluster differential gene expression statistics

Columns:

-   `analysis`: The dataset being tested.
-   `contrast`: The cluster being tested (one versus all).
-   `gene`: The gene symbol.
-   `auc`: The area under the receiver operator curve for this gene predicting cell membership to this cluster.
-   `ensembl_id`: The Ensembl gene ID.
-   `log_fc`: The log2 fold-change.
-   `ci_l`: The left side of the 95% confidence interval.
-   `ci_r`: The right side of the 95% confidence interval.
-   `ave_expr`: The mean expression of the gene in log2CPM.
-   `t`: The t-statistic.
-   `p_value`: The limma p-value.
-   `adj_p_value`: The false discovery rate (adjusted p-value).
-   `b`: The log odds ratio.

``` r
fread("de-ova.tsv.gz")
```

    ##                    analysis contrast       gene   auc      ensembl_id    log_fc
    ##       1: Tissue CD8 T cells 1 vs all      ODAPH 0.502 ENSG00000174792  1.63e-01
    ##       2: Tissue CD8 T cells 1 vs all AC068987.3 0.500 ENSG00000260473  5.31e-02
    ##       3: Tissue CD8 T cells 1 vs all AC120193.1 0.510 ENSG00000253535  4.79e-01
    ##       4: Tissue CD8 T cells 1 vs all      HPSE2 0.501 ENSG00000172987  7.29e-02
    ##       5: Tissue CD8 T cells 1 vs all     OR52I1 0.500 ENSG00000232268  9.95e-02
    ##      ---                                                                       
    ## 2459151:  Blood CD8 T cells 9 vs all     MBTPS2 0.497 ENSG00000012174 -6.66e-05
    ## 2459152:  Blood CD8 T cells 9 vs all  LINC00685 0.499 ENSG00000226179  4.13e-05
    ## 2459153:  Blood CD8 T cells 9 vs all     NLGN4Y 0.500 ENSG00000165246  2.33e-05
    ## 2459154:  Blood CD8 T cells 9 vs all AL513320.1 0.500 ENSG00000238260 -1.54e-05
    ## 2459155:  Blood CD8 T cells 9 vs all    CCDC189 0.499 ENSG00000196118  1.63e-05
    ##             ci_l   ci_r ave_expr         t  p_value adj_p_val     b
    ##       1:  0.1240 0.2020  0.03080  8.170000 7.00e-15  1.61e-10 23.20
    ##       2:  0.0396 0.0666  0.00442  7.740000 1.33e-13  1.53e-09 20.40
    ##       3:  0.3530 0.6060  0.14700  7.440000 9.23e-13  7.09e-09 18.50
    ##       4:  0.0516 0.0942  0.01090  6.740000 7.28e-11  4.19e-07 14.30
    ##       5:  0.0683 0.1310  0.01630  6.280000 1.08e-09  4.99e-06 11.70
    ##      ---                                                           
    ## 2459151: -0.2940 0.2940  0.70800 -0.000446 1.00e+00  1.00e+00 -6.41
    ## 2459152: -0.2160 0.2160  0.35600  0.000376 1.00e+00  1.00e+00 -6.41
    ## 2459153: -0.1440 0.1440  0.13200  0.000318 1.00e+00  1.00e+00 -6.41
    ## 2459154: -0.1200 0.1200  0.10600 -0.000253 1.00e+00  1.00e+00 -6.41
    ## 2459155: -0.1540 0.1540  0.21700  0.000209 1.00e+00  1.00e+00 -6.41

## Differential gene expression statistics for three contrasts

Three contrasts:

-   Case vs Control
-   CasePD1CTLA4 vs CasePD1
-   ControlPD1 vs ControlNone

Columns:

-   `analysis`: The dataset being tested.
-   `cluster`: The name of the cell cluster. All cells if it ends with `-all`.
-   `contrast`: One of the three contrasts (Case vs Control, CasePD1CTLA4 vs CasePD1, ControlPD1 vs ControlNone).
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

``` r
fread("de-contrasts.tsv.gz")
```

    ##                    analysis   cluster       gene percent    log_fc   ci_l  ci_r
    ##       1: Tissue CD8 T cells T-CD8-all       IL26   3.910  4.69e+00  3.490 5.880
    ##       2: Tissue CD8 T cells T-CD8-all     DEPDC1   2.030  2.89e+00  2.080 3.700
    ##       3: Tissue CD8 T cells T-CD8-all      IL17A   1.030  4.85e+00  3.470 6.230
    ##       4: Tissue CD8 T cells T-CD8-all AC010424.3   0.564  2.86e+00  2.040 3.680
    ##       5: Tissue CD8 T cells T-CD8-all     DLGAP5   3.220  3.54e+00  2.520 4.570
    ##      ---                                                                       
    ## 3805056:  Blood CD8 T cells   B-CD8-9       RPS3      NA  6.57e-04 -0.905 0.907
    ## 3805057:  Blood CD8 T cells   B-CD8-9     ANP32B      NA  1.01e-03 -1.500 1.500
    ## 3805058:  Blood CD8 T cells   B-CD8-9     ZNF365      NA -6.52e-04 -1.190 1.180
    ## 3805059:  Blood CD8 T cells   B-CD8-9     RPL18A      NA  2.69e-04 -0.901 0.901
    ## 3805060:  Blood CD8 T cells   B-CD8-9     HIGD1A      NA  1.18e-05 -1.480 1.480
    ##          ave_expr         t  p_value adj_p_val     b      ensembl_id
    ##       1:    2.690  8.06e+00 1.25e-08  0.000184  9.70 ENSG00000111536
    ##       2:    3.400  7.29e+00 8.11e-08  0.000316  7.98 ENSG00000024526
    ##       3:    2.670  7.22e+00 9.76e-08  0.000316  7.81 ENSG00000112115
    ##       4:    1.260  7.14e+00 1.19e-07  0.000316  7.63 ENSG00000282925
    ##       5:    4.310  7.09e+00 1.35e-07  0.000316  7.51 ENSG00000126787
    ##      ---                                                            
    ## 3805056:    7.440  1.52e-03 9.99e-01  0.999000 -5.89 ENSG00000149273
    ## 3805057:    3.650  1.41e-03 9.99e-01  0.999000 -5.89 ENSG00000136938
    ## 3805058:    0.761 -1.15e-03 9.99e-01  0.999000 -5.89 ENSG00000138311
    ## 3805059:    7.070  6.24e-04 1.00e+00  1.000000 -5.89 ENSG00000105640
    ## 3805060:    1.130  1.67e-05 1.00e+00  1.000000 -5.89 ENSG00000181061
    ##                         contrast
    ##       1:         Case vs Control
    ##       2:         Case vs Control
    ##       3:         Case vs Control
    ##       4:         Case vs Control
    ##       5:         Case vs Control
    ##      ---                        
    ## 3805056: CasePD1CTLA4 vs CasePD1
    ## 3805057: CasePD1CTLA4 vs CasePD1
    ## 3805058: CasePD1CTLA4 vs CasePD1
    ## 3805059: CasePD1CTLA4 vs CasePD1
    ## 3805060: CasePD1CTLA4 vs CasePD1

## Cell-cell communication analysis with means of gene pairs

This table has Case vs Control statistics for means of 1610 gene pairs, for each pair of cell types.

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
-   `fdr`: The false discovery rate (adjusted p-value) for gene pairs within each cluster pair.
-   `epi`: Does the interaction involve epithelial nuclei?
-   `inter`: Does the interaction involve different cell types?

``` r
fread("ccc-means.tsv.gz")
```

    ##         c1              e1 donors1        g1  c2              e2 donors2
    ##     1:   M ENSG00000243646      20    IL10RB CD8 ENSG00000111536      17
    ##     2: CD8 ENSG00000184451      21     CCR10 CD8 ENSG00000156234      14
    ##     3:   E ENSG00000243646      26    IL10RB CD8 ENSG00000111536      17
    ##     4:   E ENSG00000116329      25     OPRD1 CD8 ENSG00000156234      14
    ##     5:   M ENSG00000163702      21    IL17RC CD8 ENSG00000112115      15
    ##    ---                                                                  
    ## 13235:   B ENSG00000138448      27     ITGAV   B ENSG00000108821      12
    ## 13236:   B ENSG00000240505      27 TNFRSF13B   B ENSG00000102524      25
    ## 13237:   B ENSG00000124145      11      SDC4   E ENSG00000137801      26
    ## 13238:   B ENSG00000155760      17      FZD7   E ENSG00000085741      26
    ## 13239:   B ENSG00000106799      27    TGFBR1 CD8 ENSG00000164404      15
    ##              g2     logFC   CI.L  CI.R AveExpr         t  P.Value adj.P.Val
    ##     1:     IL26  3.614530  2.825 4.404    2.71  9.274530 3.64e-11  9.09e-10
    ##     2:   CXCL13  3.011380  2.310 3.713    1.90  8.747690 5.50e-10  1.38e-08
    ##     3:     IL26  2.403400  1.794 3.013    2.44  7.988070 1.49e-09  1.87e-08
    ##     4:   CXCL13  2.930270  2.216 3.645    1.84  8.360580 1.88e-09  4.69e-08
    ##     5:    IL17A  3.565190  2.685 4.446    2.84  8.260210 2.65e-09  4.64e-08
    ##    ---                                                                     
    ## 13235:   COL1A1  0.000380 -0.437 0.438    1.31  0.001770 9.99e-01  9.99e-01
    ## 13236: TNFSF13B -0.000556 -1.030 1.029    4.27 -0.001093 9.99e-01  9.99e-01
    ## 13237:    THBS1 -0.000494 -0.931 0.930    3.48 -0.001082 9.99e-01  9.99e-01
    ## 13238:    WNT11 -0.000162 -0.409 0.409    1.33 -0.000810 9.99e-01  1.00e+00
    ## 13239:     GDF9 -0.000127 -0.574 0.574    1.70 -0.000444 1.00e+00  1.00e+00
    ##            B                     id percent1 percent2      fdr   epi inter
    ##     1: 15.37       CT IL26 M IL10RB  0.04683 0.039067 1.27e-08 FALSE  TRUE
    ##     2: 12.77     CT CCR10 CT CXCL13  0.00312 0.016456 2.40e-07 FALSE FALSE
    ##     3: 11.68       CT IL26 E IL10RB  0.00974 0.039067 5.64e-07  TRUE  TRUE
    ##     4: 11.56      CT CXCL13 E OPRD1  0.00286 0.016456 5.64e-07  TRUE  TRUE
    ##     5: 11.23      CT IL17A M IL17RC  0.04951 0.010299 4.62e-07 FALSE  TRUE
    ##    ---                                                                    
    ## 13235: -6.04       B ITGAV B COL1A1  0.01279 0.000768 9.99e-01 FALSE FALSE
    ## 13236: -6.17 B TNFRSF13B B TNFSF13B  0.23912 0.011871 9.99e-01 FALSE FALSE
    ## 13237: -6.65         B SDC4 E THBS1  0.00052 0.173278 9.99e-01  TRUE  TRUE
    ## 13238: -5.56         B FZD7 E WNT11  0.00245 0.006829 9.99e-01  TRUE  TRUE
    ## 13239: -4.64       B TGFBR1 CT GDF9  0.02163 0.001026 1.00e+00 FALSE  TRUE

## Cell-cell communication analysis with Spearman correlation of gene percentages

This table has Spearman correlation coefficients for 1441 gene pairs, for each pair of cell types.

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
-   `fdr`: The false discovery rate (adjusted p-value) for gene pairs within each cluster pair.
-   `epi`: Does the interaction involve epithelial nuclei?
-   `inter`: Does the interaction involve different cell types?

``` r
fread("ccc-spearman.tsv.gz")
```

    ##         c1      g1 donors1  c2    g2 donors2 estimate  p.value percent1
    ##     1: CD4   PDIA3      26 CD8  CALR      27    0.866 1.10e-08 0.336965
    ##     2: CD4   KCNA3      26 CD8  IL16      27    0.840 7.74e-08 0.074123
    ##     3: CD8   ITGAL      27 CD4 ICAM3      26    0.832 1.33e-07 0.283375
    ##     4: CD8    CD8A      27   B HLA-A      27    0.822 1.42e-07 0.647686
    ##     5: CD4   ITGB2      26 CD8 ICAM3      27    0.809 5.65e-07 0.352410
    ##    ---                                                                 
    ## 10097: CD8   SCN8A      11   B FGF13      11   -1.000 1.00e+00 0.000829
    ## 10098:   E   LRRC4      15 CD8 NTNG2      16    0.000 1.00e+00 0.000783
    ## 10099:   B    CCR1      25   E CCL13      12    0.000 1.00e+00 0.011325
    ## 10100:   E IL22RA2      13 CD4  IL22      12    1.000 1.00e+00 0.000281
    ## 10101: CD8   IL2RG      27 CD4  IL21      23    0.000 1.00e+00 0.265893
    ##           mean1 percent2    mean2      fdr   epi inter               id
    ##     1: 0.364139 0.499941 0.655255 4.37e-06 FALSE  TRUE  CT CALR T PDIA3
    ##     2: 0.069246 0.332465 0.348667 1.54e-05 FALSE  TRUE  CT IL16 T KCNA3
    ##     3: 0.290258 0.403503 0.444684 4.94e-05 FALSE  TRUE CT ITGAL T ICAM3
    ##     4: 1.263167 0.834382 1.796378 7.64e-05 FALSE  TRUE  B HLA-A CT CD8A
    ##     5: 0.398149 0.395249 0.424643 7.47e-05 FALSE  TRUE CT ICAM3 T ITGB2
    ##    ---                                                                 
    ## 10097: 0.000555 0.000421 0.000443 1.00e+00 FALSE  TRUE B FGF13 CT SCN8A
    ## 10098: 0.000628 0.001263 0.001155 1.00e+00  TRUE  TRUE CT NTNG2 E LRRC4
    ## 10099: 0.013007 0.000379 0.000415 1.00e+00  TRUE  TRUE   B CCR1 E CCL13
    ## 10100: 0.000259 0.005171 0.006594 1.00e+00  TRUE  TRUE E IL22RA2 T IL22
    ## 10101: 0.281055 0.033235 0.037720 1.00e+00 FALSE  TRUE  CT IL2RG T IL21

## Correlation analysis of cluster abundances

This table has Spearman correlation coefficients for 2145 cell cluster pairs.

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
-   `fdr`: The false discovery rate (adjusted p-value) for gene pairs within each cluster pair.
-   `epi`: Does the interaction involve epithelial nuclei?
-   `inter`: Does the interaction involve different cell types?

``` r
fread("spearman-cluster-abundance.tsv.gz")
```

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

# Pseudobulk expression data

We provide log2CPM values in a matrix with one row per gene and one column per sample, along with a corresponding table that describes each sample.

## Pseudobulk at the level of cell lineage

The expression data in log2CPM:

``` r
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
```

    ## [1] 33538   249

``` r
log2cpm[1:5,1:5]
```

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##                 Tissue CD8 T cells|MC_1 Tissue CD8 T cells|MC_2
    ## ENSG00000243485                    .                      .    
    ## ENSG00000238009                    .                      0.546
    ## ENSG00000241599                    .                      .    
    ## ENSG00000235146                    .                      .    
    ## ENSG00000237491                    5.36                   5.870
    ##                 Tissue CD8 T cells|MC_9 Tissue CD8 T cells|SIC_100
    ## ENSG00000243485                    .                         .    
    ## ENSG00000238009                    .                         0.549
    ## ENSG00000241599                    .                         .    
    ## ENSG00000235146                    .                         .    
    ## ENSG00000237491                    5.13                      5.401
    ##                 Tissue CD8 T cells|SIC_109
    ## ENSG00000243485                       .   
    ## ENSG00000238009                       1.42
    ## ENSG00000241599                       .   
    ## ENSG00000235146                       .   
    ## ENSG00000237491                       5.24

The table with one row per sample:

``` r
obs <- h5read("pseudobulk_donor.h5", "obs")
str(obs)
```

    ## 'data.frame':    249 obs. of  7 variables:
    ##  $ donor  : chr [1:249(1d)] "MC_1" "MC_2" "MC_9" "SIC_100" ...
    ##  $ case   : chr [1:249(1d)] "Control" "Control" "Control" "Case" ...
    ##  $ drug   : chr [1:249(1d)] "None" "None" "None" "PD-1" ...
    ##  $ sex    : chr [1:249(1d)] "F" "M" "F" "F" ...
    ##  $ dataset: chr [1:249(1d)] "Tissue CD8 T cells" "Tissue CD8 T cells" "Tissue CD8 T cells" "Tissue CD8 T cells" ...
    ##  $ tissue : chr [1:249(1d)] "Tissue" "Tissue" "Tissue" "Tissue" ...
    ##  $ id     : chr [1:249(1d)] "Tissue CD8 T cells|MC_1" "Tissue CD8 T cells|MC_2" "Tissue CD8 T cells|MC_9" "Tissue CD8 T cells|SIC_100" ...

Complete contents:

``` r
h5ls("pseudobulk_donor.h5")
```

    ##     group     name       otype   dclass     dim
    ## 0       /   matrix   H5I_GROUP                 
    ## 1 /matrix barcodes H5I_DATASET   STRING     249
    ## 2 /matrix     data H5I_DATASET    FLOAT 3567136
    ## 3 /matrix features H5I_DATASET   STRING   33538
    ## 4 /matrix  indices H5I_DATASET  INTEGER 3567136
    ## 5 /matrix   indptr H5I_DATASET  INTEGER     250
    ## 6 /matrix    shape H5I_DATASET  INTEGER       2
    ## 7       /      obs H5I_DATASET COMPOUND     249

## Pseudobulk at the level of cell cluster

The expression data in log2CPM:

``` r
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
```

    ## [1] 33538  2958

``` r
log2cpm[1:5,1:5]
```

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##                 Tissue CD8 T cells|1|MC_1 Tissue CD8 T cells|2|MC_1
    ## ENSG00000243485                      .                         .   
    ## ENSG00000238009                      .                         .   
    ## ENSG00000241599                      .                         .   
    ## ENSG00000235146                      .                         .   
    ## ENSG00000237491                      1.25                      1.64
    ##                 Tissue CD8 T cells|3|MC_1 Tissue CD8 T cells|4|MC_1
    ## ENSG00000243485                         .                      .   
    ## ENSG00000238009                         .                      .   
    ## ENSG00000241599                         .                      .   
    ## ENSG00000235146                         .                      .   
    ## ENSG00000237491                         .                      2.35
    ##                 Tissue CD8 T cells|5|MC_1
    ## ENSG00000243485                      .   
    ## ENSG00000238009                      .   
    ## ENSG00000241599                      .   
    ## ENSG00000235146                      .   
    ## ENSG00000237491                      1.57

The table with one row per sample:

``` r
obs <- h5read("pseudobulk_donor_cluster.h5", "obs")
str(obs)
```

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

``` r
h5ls("pseudobulk_donor_cluster.h5")
```

    ##     group     name       otype   dclass      dim
    ## 0       /   matrix   H5I_GROUP                  
    ## 1 /matrix barcodes H5I_DATASET   STRING     2958
    ## 2 /matrix     data H5I_DATASET    FLOAT 24693797
    ## 3 /matrix features H5I_DATASET   STRING    33538
    ## 4 /matrix  indices H5I_DATASET  INTEGER 24693797
    ## 5 /matrix   indptr H5I_DATASET  INTEGER     2959
    ## 6 /matrix    shape H5I_DATASET  INTEGER        2
    ## 7       /      obs H5I_DATASET COMPOUND     2958

# Single-cell expression data

We provide raw counts and log2CPM values in each of the `.h5` files, along with all of the relevant metadata, as well as the TCR and BCR sequencing data.

## Tissue CD8 T cells

We provide data in the [anndata](https://anndata.readthedocs.io/en/latest/) file format.

Read the data in R:

``` r
a1 <- anndata::read_h5ad("tissue-cd8.h5ad")
a1
```

    ## AnnData object with n_obs × n_vars = 25341 × 28165
    ##     obs: 'cell', 'channel', 'case', 'class', 'class2', 'class_short', 'donor', 'facs_sorting', 'facs_gate', 'n_counts', 'n_features', 'mito_counts', 'mito_pct', 'drug', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19', 'PC20', 'PC21', 'PC22', 'PC23', 'PC24', 'UMAP1', 'UMAP2', 'leiden0.5', 'leiden0.644', 'leiden0.789', 'leiden0.933', 'leiden1.08', 'leiden1.22', 'leiden1.37', 'leiden1.51', 'leiden1.66', 'leiden1.8', 'leiden', 'barcode', 'TRAV', 'TRAD', 'TRAJ', 'TRAC', 'TRA_cdr3', 'TRA_cdr3_nt', 'TRA_length', 'TRA_reads', 'TRA_umis', 'TRBV', 'TRBD', 'TRBJ', 'TRBC', 'TRB_cdr3', 'TRB_cdr3_nt', 'TRB_length', 'TRB_reads', 'TRB_umis', 'TRB_cdr3_trim', 'TRA_cdr3_trim', 'TRAV2', 'TRBV2', 'TRAV_cdr1', 'TRAV_cdr2', 'TRAV_cdr25', 'TRBV_cdr1', 'TRBV_cdr2', 'TRBV_cdr25', 'has_tcr', 'IGHV', 'IGHD', 'IGHJ', 'IGHC', 'IGH_cdr3', 'IGH_cdr3_nt', 'IGH_length', 'IGH_reads', 'IGH_umis', 'IGKV', 'IGKD', 'IGKJ', 'IGKC', 'IGK_cdr3', 'IGK_cdr3_nt', 'IGK_length', 'IGK_reads', 'IGK_umis', 'IGLV', 'IGLD', 'IGLJ', 'IGLC', 'IGL_cdr3', 'IGL_cdr3_nt', 'IGL_length', 'IGL_reads', 'IGL_umis', 'has_bcr', 'cluster'
    ##     var: 'mean', 'sd', 'percent', 'gene', 'exclude', 'include', 'fitted', 'residuals', 'rank'
    ##     uns: 'de', 'knn', 'mcv', 'pca', 'pca_h'

Read the data in Python:

``` python
import anndata
a1 = anndata.read_h5ad("tissue-cd8.h5ad")
a1
```

    ## AnnData object with n_obs × n_vars = 25341 × 28165
    ##     obs: 'cell', 'channel', 'case', 'class', 'class2', 'class_short', 'donor', 'facs_sorting', 'facs_gate', 'n_counts', 'n_features', 'mito_counts', 'mito_pct', 'drug', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18', 'PC19', 'PC20', 'PC21', 'PC22', 'PC23', 'PC24', 'UMAP1', 'UMAP2', 'leiden0.5', 'leiden0.644', 'leiden0.789', 'leiden0.933', 'leiden1.08', 'leiden1.22', 'leiden1.37', 'leiden1.51', 'leiden1.66', 'leiden1.8', 'leiden', 'barcode', 'TRAV', 'TRAD', 'TRAJ', 'TRAC', 'TRA_cdr3', 'TRA_cdr3_nt', 'TRA_length', 'TRA_reads', 'TRA_umis', 'TRBV', 'TRBD', 'TRBJ', 'TRBC', 'TRB_cdr3', 'TRB_cdr3_nt', 'TRB_length', 'TRB_reads', 'TRB_umis', 'TRB_cdr3_trim', 'TRA_cdr3_trim', 'TRAV2', 'TRBV2', 'TRAV_cdr1', 'TRAV_cdr2', 'TRAV_cdr25', 'TRBV_cdr1', 'TRBV_cdr2', 'TRBV_cdr25', 'has_tcr', 'IGHV', 'IGHD', 'IGHJ', 'IGHC', 'IGH_cdr3', 'IGH_cdr3_nt', 'IGH_length', 'IGH_reads', 'IGH_umis', 'IGKV', 'IGKD', 'IGKJ', 'IGKC', 'IGK_cdr3', 'IGK_cdr3_nt', 'IGK_length', 'IGK_reads', 'IGK_umis', 'IGLV', 'IGLD', 'IGLJ', 'IGLC', 'IGL_cdr3', 'IGL_cdr3_nt', 'IGL_length', 'IGL_reads', 'IGL_umis', 'has_bcr', 'cluster'
    ##     var: 'mean', 'sd', 'percent', 'gene', 'exclude', 'include', 'fitted', 'residuals', 'rank'
    ##     uns: 'de', 'knn', 'mcv', 'pca', 'pca_h'
