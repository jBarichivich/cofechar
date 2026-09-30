# cofechar

R port of Holmes' COFECHA and EDT modules from the Dendrochronology Program
Library (DPL). Provides quality control and crossdating for tree-ring
measurement series, with a complete dplR integration layer.

## Installation

```r
remotes::install_github("jbarichivich/cofechar")
```

## Example

```r
library(cofechar)

# Load the included Fitzroya cupressoides dataset
# (Cerro Mirador, Alerce Costero, Barichivich 2005)
mir_file <- system.file("extdata", "CL-MIR.rwl", package = "cofechar")
rwl <- dpl_read_dec(mir_file, stop_val = -9999L, unit = "0.001mm")

# COFECHA quality control
cof <- dpl_cof(rwl)
cof$stats
subset(cof$segments, flag != "")
```

## Reference

Holmes, R. L. (1983). Computer-assisted quality control in tree-ring dating
and measurement. *Tree-Ring Bulletin* 43:69--78.
