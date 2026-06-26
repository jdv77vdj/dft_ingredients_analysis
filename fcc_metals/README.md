In this folder we have several files to analyse the ingredients of the pristine and defective (monovacant) fcc metals.

## Folders with ingredients

Since the folders are heavy we did not upload them, please send an email to [jvega5@tulane.edu](mailto:jvega5@tulane.edu) for requesting the computed ingredient files. Also you can see the files for Si-diamond, which are available in this repository.

The metals analyzed are Al, Ni, Pd, Pt, Ag, Au, Pt, Pb.

For each metal we have the next structure:

```text
metal_name/
├── metal_name
└── metal_name_v
```

The subfolder named `metal_name` has the pristine structure and the subfolder `metal_name_v` has the defective structure with a monovacancy.

All folders contain:

* `AECCAR0`, `AECCAR2`, `CHGCAR`, `ELFCAR`, `POSCAR` files coming from VASP.
* `data_2D.dat`: contains the computed ingredients in 2D from the generator file `2d_ingredients_generator_metals.py`.
* `1d_data_path_ae.dat`: contains computed ingredients along one line (or path) from the `1d_ingredients_generator_metals_given_2_positions.py` file or `1d_ingredients_generator_given_3_positions.py`.
* Note that for running the generator `.py` files, `generate_Fx.py` must be present.
* There is also one script named `1d_ingredients_generator_metals_original_version.py` which works for the 1d case, although this is the first version of the code and only works for 2 atoms.

## Scripts

We have one set of scripts for generating the DFAs ingredients and another set of scripts for plotting them. These scripts are in 2 folders: `1d` and `2d`.

The workflow to generate ingredients is the next:

1. Generate ingredients either in 1d or 2d using the appropriate `.py` script.
2. Plot ingredients with `.py` script accordingly.

## Excel file with defect formation energies

It contains the defect formation energies obtained with all DFAs for the fcc metals. This file is needed to generate the bar plots shown in the `\Delta R_{xc}` figures.
