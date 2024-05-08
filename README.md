
<!-- README.md is generated from README.Rmd. Please edit that file -->

# AMEND 2.0 Manuscript

<!-- badges: start -->
<!-- badges: end -->

This repository contains all of the necessary code and data to recreate
the analysis results from the manuscript titled “AMEND 2.0: Active
Module Identification and Multi-Omic Data Integration with
Multiplex-Heterogeneous Networks”.

AMEND (Active Module identification with Experimental data and Network
Diffusion) is an algorithm designed to find a subset of connected nodes
in a molecular interaction network that have large experimental values.
It makes use of random walk with restart (RWR) to create node weights,
and a heuristic approach for solving the Maximum-weight Connected
Subgraph problem using these weights. This is performed iteratively
until an optimal subnetwork (i.e., module) is found. AMEND can now
accommodate multiplex and/or heterogeneous networks, making it a widely
applicable tool. These complex networks can include several node types,
edge types, and/or data types, which increases the types of biological
questions one can address with this method.

## Organization

There are five main analyses that were performed for this study. Random
Walk with Restart for Multiplex/Heterogeneous Networks (RWR-MH) was
assessed on it’s ability to rank nodes related to the seed nodes using
KEGG pathways. Biased Random Walk (BRW) was assessed on it’s ability to
rank nodes related to the seed nodes and on it’s ability to affect the
final module in the context of AMEND. Various degree bias adjustment
methods and transition matrix types were assessed on their ability to
rank nodes related to seed nodes and on their ability to reduce the
correlation between diffusion scores and node degree. Finally, AMEND was
applied to two real-world multi-omic datasets (an O-GlcNAc Transferase
knockout study and the TCGA-KIRC project) to assess it’s ability to
return biologically relevant results.

## Authors

Sam Boyd, Chad Slawson, and Jeff Thompson.

<!-- ## Example -->
<!-- AMEND contains three objects -->
<!-- ```{r example} -->
<!-- library(AMEND) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
