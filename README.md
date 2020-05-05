

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
2. Install the package called deltanetPlus, using `devtools::install_github("CABSEL/DeltaNeTSplus/deltanetsPlus")`. Your package is inatalled in R library directory.
3. Load the package, using `library(deltanetsPlus)`.



### Acknowledgements
This work has been supported by the ETH Zurich Research Grant.



