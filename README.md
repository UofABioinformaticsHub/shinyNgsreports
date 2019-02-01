# fastqcRShiny

Graphical User Interface for inspecting FastQC reports using the [ngsReports](https://github.com/UofABioinformaticsHub/ngsReports) R package. 

## Installation
To install required packages and the ngsReports base package follow the instructions below.
Currently you need to install the fastqcTheoreticalGC package separately.

```
source("https://bioconductor.org/biocLite.R")
pkgs <- c("BiocGenerics", "BiocStyle", "BSgenome", "checkmate", "devtools", "ggdendro",  "plotly", "reshape2", "Rsamtools", "scales", "ShortRead", "tidyverse",  "viridis", "viridisLite", "zoo", "mikelove/fastqcTheoreticalGC", "UofABioinformaticsHub/ngsReports")
BiocManager::install(pkgs)
BiocManager::install('UofABioinformaticsHub/fastqcRShiny')

library(ngsReports)
library(fastqcRShiny)
```

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
