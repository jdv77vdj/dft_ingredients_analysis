In this folder we have several files to analyse the ingredients of the pristine Si-diamond and defective Si-diamond interstitials. We only analyzed the T interstitial.

## Folders with ingredients

* `Si`: Contains the Silicon diamond pristine outputs.
* `T`: Contains the Tetrahedral interstitial T.

All folders contain:

* `AECCAR0`, `AECCAR2`, `CHGCAR`, `ELFCAR`, `POSCAR` files coming from VASP.
* `data_2D.dat`: contains the computed ingredients in 2D from the generator file `2d_ingredients_generator_silicon.py`.
* `1d_data_path_ae.dat`: contains computed ingredients along one line (or path) from the `1d_ingredients_generator_silicon_given_2_positions.py` file.
* Note that for running the generator `.py` files, `generate_Fx.py` must be present.

## Scripts

We have one set of scripts for generating the DFAs ingredients and another set of scripts for plotting them. These scripts are in 2 folders: `1d` and `2d`.

The workflow to generate ingredients is the next:

1. Generate ingredients either in 1d or 2d using the appropriate `.py` script.
2. Plot ingredients with the `.py` script accordingly.

## Excel file with defect formation energies

It contains the defect formation energies obtained with all DFAs. While in the scripts for plotting, we don't need this Excel file, if you'd like to add a bar of defect formation energies, having the Excel file might be useful.
