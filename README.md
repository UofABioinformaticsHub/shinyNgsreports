# shinyNgsreports

Graphical User Interface for inspecting ngs log files using the [ngsReports](https://github.com/UofABioinformaticsHub/ngsReports) R package. 
Currently only supports Fastqc reports

## Installation
To install required packages and the ngsReports base package follow the instructions below.

```
install.packages("BiocManager")
BiocManager::install('UofABioinformaticsHub/shinyNgsreports')
```
# Vignette
A simple user guide can be found [here](https://cmwbio.github.io/fastqcReports-vignette/index.html).

# Quick start quide
After loading both the base package `ngsReports` and the shiny app `ngsReportsShiny` the shiny app can simply be run using the 
function `fastqcShiny().`

```
library(ngsReports)
library(shinyNgsreports)

fastqcShiny()
```
Data can then be imported by clicking the "Choose Files" button and navigating to your fastqc files.
Plots will then be loaded automatically.

# Citation 

Please cite the ngsReports [preprint](https://www.biorxiv.org/content/early/2018/05/02/313148):

```
@article{ward2018ngsreports,
  title={ngsReports: An R Package for managing FastQC reports and other NGS related log files.},
  author={Ward, Christopher M and To, Hien and Pederson, Stephen M},
  journal={bioRxiv},
  pages={313148},
  year={2018},
  publisher={Cold Spring Harbor Laboratory}
}
```
