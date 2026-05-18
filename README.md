# DisCO.VG

Integrative tool to predict genetic variant-target gene links in a specific disease and cell type/tissue.

NOTE: 
- The first release is pending final adjustments to features.
- This repository contains as yet unpublished work.

## Description

DisCO-VG (Disease-specific CMMC Optimization of Variant-Gene pairs) is a computational method to link disease-associated variants to their target genes with cell type or tissue specificity. DisCO-VG uses a novel statistical learning method called *coupled matrix-matrix completion* (CMMC) ([Bagherian et al. (2021)](10.1093/bib/bbaa025)) to generate the variant-gene linking scores. By integrating various sources of functional (epi)genomic information such as chromatin accessibility, quantitative trait loci, and gene pathway information, DisCO-VG computes a score that represents the strength of variant-gene targeting given a specific disease and cell type/tissue.

The goal of this tool is to nominate variant-gene pairs in a disease-specific and cell-type-/tissue-specific manner that can be further validated with experimental approaches. By systematically mapping GWAS variants to their target genes, we come closer and closer to building a future in which our genomes inform our respective personalized healthcare.

![DisCO-VG overview](https://github.com/sartorlab/DisCO.VG/blob/main/fig/DisCO-VG_overview.png?raw=true)

**Inputs:** *GWAS variants*, *ATAC-Seq* for cell type or tissue of interest (bulk or single cell), and *eQTL* data for cell type or tissue of interest (bulk or single cell).

**Outputs:** DisCO-VG scores for all possible variant-gene pairs based on input data. The score ranges from [0, 1] and reflects the strength of disease-specific, cell-type- or tissue-specific targeting. (Note: preprocessed DisCO-VG scores and analytics of tool performance may also be available in a later version.)

## References

Bagherian M, Kim RB, Jiang C, Sartor MA, Derksen H, Najarian K. Coupled matrix-matrix and coupled tensor-matrix completion methods for predicting drug-target interactions. Brief Bioinform. 2021 Mar 22;22(2):2161-2171. doi: 10.1093/bib/bbaa025. PMID: 32186716; PMCID: PMC7986629.

(In preparation) E.Chou*, A. Guerra*, Z. Zhang, T. Qin, S. Li, H. Zhang, K. Wang, R. Sherpa, A.P. Boyle, J.T. Elder, L.C. Tsoi, M.A. Sartor. Prioritizing noncoding variant-gene pairs in psoriasis using a data fusion approach.
