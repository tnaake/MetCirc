# MetCirc

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![license](http://img.shields.io/badge/license-GPL%20%28%3E=%203%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![BioC checks](https://bioconductor.org/shields/build/devel/bioc/MetCirc.svg)](https://bioconductor.org/shields/build/devel/bioc/MetCirc.svg)

Navigating mass spectral similarity in high-resolution MS/MS metabolomics data

## Description
Please visit [MetCirc](https://bioconductor.org/packages/MetCirc) for further information. 

## Contact 

You are welcome to 

 * write a mail to <thomasnaake (at) googlemail (dot) com> 
 * submit suggestions and issues: <https://github.com/tnaake/MetCirc/issues>
 * send a pull request: <https://github.com/tnaake/MetCirc/issues> 

## Install
To install MetCirc, please use the stable version available via Bioconductor. 
To install, enter 

```r 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MetCirc")
``` 

to your console. The installation via BiocManager requires R version 3.6. 


If you would like to install the development version of MetCirc, you will first
have to install [devtools](http://cran.r-project.org/web/packages/devtools/index.html) package: 

```r
install.packages("devtools")
library("devtools")
install_github("tnaake/MetCirc")
```


