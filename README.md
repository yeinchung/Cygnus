<p align="center"><img src="figures/image.png" alt="" width="350"></a></p>
<hr>

### Cygnus Version 0.1.0

**Single EV Imaging Data Analysis and Visualization Pipeline**

Cygnus offers data analysis tool for ___. 

<hr>

# Quick Installation of Cygnus

**First, install devtools (for installing GitHub packages) if it isn't already installed:**
``` r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
```

**Then, install BiocManager (for installing bioconductor packages) if it isn't already installed:**
``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```

**Then, install Cygnus:**
``` r
devtools::install_github("yeinchung/Cygnus", ref="master", repos = BiocManager::repositories())
```

# Issues using Cygnus?

Cygnus is currently in __beta__. We expect there to be bumps in the road. If you think you have found a bug, please first install the latest version of Cygnus via
``` r
devtools::install_github("yeinchung/Cygnus", ref="master", repos = BiocManager::repositories())
```
If this does not fix your problem, please [report an issue on Github](https://github.com/yeinchung/Cygnus/issues) with the __Bug Report__ form.
