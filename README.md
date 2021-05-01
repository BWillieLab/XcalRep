# XcalRep

Implementation of cross-calibration and precision analysis for HR-pQCT. 

Please cite: Mikolajewicz et al (2020) Bone. \nDOI: 10.1016/j.bone.2021.115880
https://pubmed.ncbi.nlm.nih.gov/33561589/ 

## Description

[XcalRep](https://nmikolajewicz.github.io/XcalRep/) (x-calibration & reproducibility) is an R implementation of methods used for cross-calibration and precision analysis of HR-pQCT It enables cross-calibration and precision analysis in R for single- and multi-centre studies, conducted at single or multiple timepoints. XcalRep was developed with intended use for HR-pQCT, however, it easily extends to other instruments which undergo routine calibration and reproducibility assessment. 

The XcalRep package provides the following utilities:
* single and multi-site precision error calculation (CV and SD-based errors)
* reference site selection
* scanner cross-calibration

Additionally, the XcalRep package provides multiple visualization tools to evaluate scanner reproducibility. 

Visit [XcalRep package site](https://nmikolajewicz.github.io/XcalRep/) for documentation and step-by-step tutorials/vignettes. 

## Installation

Firstly, please install or update the package devtools by running

```
install.packages("devtools")
```

Then the XcalRep can be installed via

```
library(devtools)
devtools::install_github(repo = "NMikolajewicz/XcalRep")
```
## Authors

* [Nicholas Mikolajewicz](https://scholar.google.ca/citations?user=LBWQMXsAAAAJ&hl=en&oi=ao)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Bettina Willie, Elizabeth Zimmermann
