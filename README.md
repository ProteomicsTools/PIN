# PIN

### Installation of development version from GitHub

```sh
$ library(devtools)
$ install_github("shaowenguang/PIN")
```

### Installation of release version from CRAN

Needs to wait for the official release...

```sh
$ install.packages("PIN")
```

### Description

This document describes a computational method, namely PIN (proteome integrity number), to assess a quantitative measure of protein stability directly from bottom-up proteomic datasets.

By quantifying the relative abundance of semi-tryptic peptides (i.e. the likely products of protein degradation), an individual protein integrity score (iPIS) was first calculated for each measured protein. A per-sample PIN was then derived to summarize a proteome-wide measurement, indicating the level of protein degradation.

This method has been applied to several clinical cohorts of prostate tissue samples, breast cancer tissue samples, gastric cancer samples and human plasma samples, respectively in which degradation is sometimes inevitable during the sample preparation. The datasets analyzed also represent the types of most commonly used clinical specimens, fresh frozen tissue, FFPE tissue and plasma and the most common mass spectrometric acquisition methods. Overall, the data indicate that the PIN algorithm is broadly applicable for the assessment for the proteome integrity in clinical studies.  


### Usage



First: to load the package.
```sh
$ library(PIN)
```

Second: to extract the peptide list, together with their quantities, from OpenSWATH outputs.
```sh
$ peptides <- generate_peptide_table(search_results="./openswath_search_results.tsv", sample_annotation="./sample_annotation_table", sptxt="./spectral_library.sptxt")
```

Third: to perform PIN analysis.
```sh
$ perform_PIN_analysis(peptides)
```
#### Outputs
| File Name | Description |
| ------ | ------ |
| PIN.tsv | PINs and their p-values, for each sample |
| iPIS.tsv | individual Protein Integrity Score, for each individual proteins |
| lib_peptide.tsv | a temp file that records the assay library information |





### Todos

 -  
 -  

License
----
