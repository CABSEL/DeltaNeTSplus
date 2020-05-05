

# DeltaNeTS+

- DeltaNeTS+ is a major improvement to our previous method DeltaNet.
- DeltaNeTS+ is a method for inferring direct gene targets of drug compounds and diseases from steady-state and/or __time-series__ transcriptional profiles. 
- DeltaNeTS+ incorporates gene regulatory information (if avaialbe) during the inference. 
- DeltaNeTS+ generates a perturbation score for each gene in every sample. The score magnitude reflects the confidence that the transcription process of this gene was directly affected by the external stimuli. 
- The score sign indicates the nature of the perturbation: positive for gene induction, negative for gene repression.

Please refer to [DeltaNeTS+ manuscript](https://www.biorxiv.org/content/10.1101/788968v1) for more detailed information.

### DeltaNeTS+ installation
To install `DeltaNeTS+` directly from github repository, `devtools` R package is required. 

1. Install and load `devtools` package in R.
2. Install the package called `deltanetsPlus`, using `devtools::install_github("CABSEL/DeltaNeTSplus/deltanetsPlus")`. Your package will e inatalled in R library directory.
3. Load the package, using `library(deltanetsPlus)`, and now you're ready to use!


### DeltaNeTS+ example codes using *C.elegans* expression data

`deltanetsPlus` package includes `lfc`, `pval`, `tp`, `experiment`, and `grn` data of *C.elegans*, which were processed from the orignal data of Baugh et al. 2005 (GSE2180).
- `lfc`: log2Fc of differential gene expression values. Rows are genes and Columna are time-series samples from three different experiments (10 time points per each experiment).  
- `pval`:
- `tp`:
- `experiment`:
- `grn`:


```{r warning=FALSE,eval=FALSE,echo=TRUE}
pgn <- generatePGN(glist = glist, tftg = tftg, ppi = ppi, tftg_thre = 0, ptf_thre = 0, 
                   ppi_thre = 500)
```
d.obj = createDeltanetsPlusObj(lfc=lfc, pval=pval, tp=tp, experiment=experiment, p.thres=0.05)

## in case of combining multiple datasets for the same organism or cell line, one can use combin2() for deltaents object. e.g. d.obj=cbind2(d.obj1,d.obj2)


grn = read.delim(file="./celegans_data/grn_celegans.txt")

### Acknowledgements
This work has been supported by the ETH Zurich Research Grant.



