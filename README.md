# Ingredients Analysis Package

This repository contains the scripts and data used to analyze the ingredients of semilocal density functional approximations for fcc metals and semiconductors, as presented in

> Vega Bazantes, J., *et al.* "Rationalizing defect formation energies in metals and semiconductors with semilocal density functionals." *arXiv* (2026), arXiv:2604.05385.

The corresponding preprint is available at: https://arxiv.org/pdf/2604.05385

The core implementation was originally developed by **Timo Lebeda** and subsequently extended by **Jorge Vega Bazantes** for the analysis presented in the paper.

This package includes 2 main folders:

* Si_diamond
* fcc_metals

Each folder contains Python codes for generating the ingredients and then plotting them in 1D and 2D. Each folder has its own README and also the scripts explain what the inputs and ouputs are. All the scripts were run on Python 3.12.2. The VASP outputs were generated with VASP 6.5.1. In this repository only the Si_diamond folder has the ouputs from VASP to generate the figures. For the fcc metals' ingredients contact the authors. 
