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

## Licence and funding

`cofechar` is released under the GNU General Public License v3 (or later).

`cofechar` was developed by Prof. Jonathan Barichivich (CNRS-LSCE) as a
reimplementation in R of the original Fortran 77 source of COFECHA and EDT
written by Richard L. Holmes (Laboratory of Tree-Ring Research, University of
Arizona). Development was funded by the European Research Council (ERC) under the European Union's
Horizon Europe programme, Starting Grant **CATES** (*Long-term consequences of
altered tree growth and physiology in the Earth System*, grant agreement
no. [101043214](https://cordis.europa.eu/project/id/101043214)), hosted by CNRS
at the Laboratoire des Sciences du Climat et de l'Environnement (LSCE).

`cofechar` is part of **xDPL**, the CATES suite of open-source tools for
modern, reproducible and scriptable tree-ring research workflows.

## Citation

Run `citation("cofechar")` in R, or use the "Cite this repository" button on
GitHub. Please also cite Holmes, R. L. (1983), Tree-Ring Bulletin 43:69--78,
whose COFECHA program this package reimplements.
