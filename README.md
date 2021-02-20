

# DeltaNeTS+

- DeltaNeTS+ is a major improvement to our previous method DeltaNet.
- DeltaNeTS+ is a method for inferring direct gene targets of drug compounds and diseases from steady-state and/or __time-series__ transcriptional profiles. 
- DeltaNeTS+ incorporates gene regulatory information (if avaialbe) during the inference. 
- DeltaNeTS+ generates a perturbation score for each gene in every sample. The score magnitude reflects the confidence that the transcription process of this gene was directly affected by the external stimuli. 
- The score sign indicates the nature of the perturbation: positive for gene induction, negative for gene repression.

Please refer to [DeltaNeTS+ manuscript](https://www.biorxiv.org/content/10.1101/788968v1) for more detailed information.

### DeltaNeTS+ installation
To install `DeltaNeTS+` directly from github repository, `devtools` R package is required. 

1. Install and load `devtools` package in R using the follwoing commands:

```{r warning=FALSE,eval=FALSE,echo=TRUE}
## installing devtools from CRAN
install.packages("devtools")
## Loading devtools library
library(devtools)
```
2. Install the package called `deltanetsPlus`, using `devtools::install_github("CABSEL/DeltaNeTSplus/deltanetsPlus")`. When you are asked to install/update the denpendencies, type 1 for 'All'. Your package will be inatalled in R library directory.
3. Load the package, using `library(deltanetsPlus)`, and now you're ready to use!


### DeltaNeTS+ example codes using *Caenorhabditis elegans* expression data

`deltanetsPlus` package includes example data of `lfc`, `pval`, `tp`,and `experiment`, which were processed from *C. elegans* dataset of Baugh et al. 2005 (GSE2180), as well as the `grn` structure (tf-gene interactions) of *C. elegans*. This dataset consists of 30 time-series samples of three genetic perturbation experiments in *C. elegans* (10 time points in each experiment), and the gene knock-downs were mex-3 for Experiment A, pie-1 for Experiment B, and pie-1 and pal-1 for Experiment C.

- `lfc`: an nx30 matrix of log base2 Fold Change values of differential gene expressions. Rows are genes (n) and Columna are time-series samples from three different experiments (10 time points per each experiment).  
- `pval`: an nx30 matrix of statistical p-values of the differential gene expressions.
- `tp`: a 1X30 vector of sample time points
- `experiment`: a 1x30 vector of group/experiment indication.
- `grn`: TF-gene interaction list.

#### 1. Generate a deltanetsPlus object.

In this example, we are creating a DeltanetsPlus object using `createDeltanetsPlusObj()`, which will filter out unsignificant gene expressions and compute slopes of gene expression changes, given `pval` and `tp`, respectively. 
 
```{r warning=FALSE,eval=FALSE,echo=TRUE}
d.obj = createDeltanetsPlusObj(lfc=lfc, pval=pval, tp=tp, experiment=experiment, p.thres=0.05)
```

One can use `combin2()` to combine two DeltanetsPlus objects (e.g. `d.obj=cbind2(d.obj1,d.obj2)`). Note that the genes should match between the two datasets.

#### 2. Compute gene perturbation scores using deltanetPlus().

In this example, we will compute the perturbations for only a select set of genes.

```{r warning=FALSE,eval=FALSE,echo=TRUE}
gset = c("mex-3","pie-1","pal-1",sample(rownames(lfc),10))
```
In the above, `gset` includes the actual perturbation targets (mex-3, pie-1, and pal-1) as well as 10 random genes, which were not perturbed. In the manuscript describing DeltaNeTS+, we applied the method for the complete set of genes. If `gset` is not provided, `deltanetsPlus()` will compute the perturbation scores for the complete set of genes. The example below implements parallel computing (`par=TRUE`) with 2 nodes (`numClusters=2`). To turn off parallel option, users can set `par=FALSE`.

```{r warning=FALSE,eval=FALSE,echo=TRUE}
dts.res <- deltanetsPlus(d.obj, 
                         grn=grn, 
                         kfolds=10,
                         cv.method="cv",
                         perturbation="group",
                         gset=gset,
                         group=NULL,
                         lambda=10^seq(-2,5,length.out=100),
                         cv.opt = "lambda.1se",
                         par=TRUE, numClusters=2)
```

Finally, we can check the perturbation scores by DeltaNeTS+. An example of the result from the analysis above is shown below. The result expectedly shows large negative perturbation values for mex-3, pie-1, and pal-1/pie-1 in exp. A (1st Column), exp. B (2nd Column), and exp. C (3rd column), respectively, since these genes were knocked down in the experiments. 


```{r warning=FALSE,eval=FALSE,echo=TRUE}
gset2 = intersect(gset, rownames(dts.res$P))
print(dts.res$P[gset2,])
```

Example results:
##### __12 x 3 sparse Matrix of class "dgCMatrix"__


genes    |             1|           2|           3
--- | --- | --- | --- 
mex-3    |  **-1.776621450**| -0.59719729| -0.62503056
pie-1    |  -1.182894234| **-2.05782894**| **-2.40294302**
pal-1    |   .          | -0.02866807| **-3.01030870**
ife-3    |  -0.002458885| -0.00233180| -0.00266895
frm-5.1  |   0.044570707| -0.08326815|  .         
dnj-17   |  -0.138966291| -0.30727911| -0.26977825
F29D10.2 |  -0.187006724| -0.17008810| -0.18143735
pgrn-1   |  -0.307675932|  .         |  0.10154587
F25E5.3  |  -0.093054526|  .         | -0.12077273
Y69A2AR.23|  0.058642704|  .         |  .         
ZK112.6   |  0.084273230|  .         |  0.16039465
ncs-2     | -0.104642724|  .         |  .         



### Acknowledgements
This work was supported by the ETH Zurich Research Grant.



