#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd

from matplotlib.lines import Line2D

from matplotlib.patches import Patch



# =========  Global style controls  =========
#FIGSIZE = (6.2, 8.5)
FIGSIZE = (11, 4.5)
BASE_FONT = 16

plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": BASE_FONT,
    "axes.linewidth": 1.2,
    "xtick.labelsize": BASE_FONT,
    "ytick.labelsize": BASE_FONT,
})

### T, H, X ### Theoretical
#formation_data = {
#    "LDA":   [3.42, 3.39, 3.33],
#    "PBE":   [3.76, 3.60, 3.55],
#    "SCAN":  [4.56, 4.30, 4.21],
#    "r2SCAN":[4.50, 4.16, 4.21],
#    "LAK":   [5.33, 4.79, 4.88]
#}

### T, H, X ###
formation_data = {
    "LDA":   [3.70, 4.09, 21.97],
    "PBE":   [4.26, 4.66, 22.98],
    "SCAN":  [4.96, 5.18, 23.54],
    "r2SCAN":[4.89, 5.05, 23.58],
    "LAK":   [5.73, 5.72, 24.62]
}

methods = list(formation_data.keys())
defect_index = {
    "T":0,   # T
    "H":1,   # H
    "X":2    # X
}




# ======================================================
# Master plotting function
# ======================================================
def plot_defect(metal_name, defect_name, limits):
    ylim_left = limits["ylim_left"]
    ylim_right = limits["ylim_right"]
    xlim_left = limits["xlim_left"]
    xlim_right = limits["xlim_right"]
    
    Ef_range = limits["Ef_range"]
    x1 = limits["reg_x1"]
    x2 = limits["reg_x2"]
    x3 = limits["reg_x3"]
    x4 = limits["reg_x4"]
    
    
    x_limit = 41
    fig, (ax_left, ax_right) = plt.subplots(
    1, 2, figsize=FIGSIZE, sharex=True
    )

    # ---------------------------------------------------
    # Load data
    # ---------------------------------------------------
    
    base_path = f"/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_March7/Silicon"
    file_pr = f"{base_path}/{metal_name}/data_path_ae.dat"
    file_def = f"{base_path}/{defect_name}/data_path_ae.dat"

    data_pr = np.loadtxt(file_pr, comments="#")[:x_limit]
    data_def = np.loadtxt(file_def, comments="#")[:x_limit]

    x = data_pr[:, 1]
    n_pr, n_def = data_pr[:, 4], data_def[:, 4]
    ones_pr, ones_def = np.ones_like(x), np.ones_like(x)


    Fx_pr = [data_pr[:, 9], data_pr[:, 10], data_pr[:, 11], data_pr[:, 12], ones_pr]
    Fc_pr = [data_pr[:, 13], data_pr[:, 14], data_pr[:, 15], data_pr[:, 16], data_pr[:, 17]]
    Fx_def = [data_def[:, 9], data_def[:, 10], data_def[:, 11], data_def[:, 12], ones_def]
    Fc_def = [data_def[:, 13], data_def[:, 14], data_def[:, 15], data_def[:, 16], data_def[:, 17]]
    
    # xc energy density
    exc_pr = [n_pr**(4/3) * (Fx + Fc) for Fx, Fc in zip(Fx_pr, Fc_pr)]
    exc_def = [n_def**(4/3) * (Fx + Fc) for Fx, Fc in zip(Fx_def, Fc_def)]

    # Ratios
   # ratio_pbe, ratio_scan, ratio_r2scan, ratio_lak, ratio_lda = [d/p for d, p in zip(exc_def, exc_pr)]
    ratio_pbe, ratio_scan, ratio_r2scan, ratio_lak, ratio_lda = [p/d for d, p in zip(exc_def, exc_pr)]
    rs_pr, s_pr, a_pr = data_pr[:, 5], data_pr[:, 7], data_pr[:, 8]
    rs_def, s_def, a_def = data_def[:, 5], data_def[:, 7], data_def[:, 8]

    # Differences of the ingredients
    ds  = s_def  - s_pr
    da  = a_def  - a_pr
    drs = rs_def - rs_pr

    
    ############################################
    # Left plot: Delta Ratios
    ############################################

    # --- core / semicore shading ---
    color_base = "orange"
    ax_left.axvspan(x1[0], x1[1], color=color_base, alpha=0.40, zorder=0)
    ax_left.axvspan(x2[0], x2[1], color=color_base, alpha=0.15, zorder=0)
    ax_left.axvspan(x3[0], x3[1], color="#E65100", alpha=0.20, zorder=0)
    ax_left.axvspan(x4[0], x4[1], color=color_base, alpha=0.20, zorder=0)

    reference = ratio_lda
    lwidth=3;
    ax_left.plot(x, -(ratio_pbe    - reference), color="tab:blue",   lw=lwidth, ls="-", label="PBE")
    ax_left.plot(x, -(ratio_scan   - reference), color="tab:orange", lw=lwidth, ls="-", label="SCAN")
    ax_left.plot(x, -(ratio_r2scan - reference), color="tab:green",  lw=lwidth, ls="-", label="r2SCAN")
    ax_left.plot(x, -(ratio_lak    - reference), color="tab:cyan",   lw=lwidth, ls="-", label="LAK")
    ax_left.plot(x, -(ratio_lda    - reference), color="tab:purple", lw=lwidth, ls="-", label="LDA")
    
    ax_left.axhline(0.0, color="k", lw=0.8, ls="--")
    ax_left.set_ylabel(r"$\Delta(1/R_{xc})$", fontsize=20)
    ax_left.set_ylim(*ylim_left)
    ax_left.set_xlim(*xlim_left)
    ax_left.set_title(defect_name, fontsize=34, loc="left", y=0.9, pad=-14)
    ax_left.set_xlabel(r"$r (\mathrm{\AA})$", fontsize=20)
    #ax_left.set_title(defect_name, fontsize=30, loc="left", pad=-10)
    
    # --- Dummy patch for DMC ---
    dmc_patch = Patch(facecolor="gold", edgecolor="k", alpha=0.65, label="DMC")
    
    handles, labels = ax_left.get_legend_handles_labels()
    ax_left.legend(handles + [dmc_patch], labels + ["DMC"],
               frameon=True,
               ncol=3,
               loc="lower center",
               bbox_to_anchor=(0.5, 1.02),
               fontsize=BASE_FONT)

    
  
    # ===== Inset for formation energies =====

    axins3 = inset_axes(ax_left, width="18%", height="30%", loc="upper right")

    idx = defect_index[defect_name]

    energies = [formation_data[m][idx] for m in methods]

    colors = {
    "LDA":"tab:purple",
    "PBE":"tab:blue",
    "SCAN":"tab:orange",
    "r2SCAN":"tab:green",
    "LAK":"tab:cyan",
    "HSE06":"#8B4513",
    "RPA":"orange",
    "CCSD":"silver",
    "CCSD(T)":"teal"
    }

    x_bar = np.arange(len(methods))
    axins3.bar(x_bar, energies, color=[colors[m] for m in methods], width=0.7)

    # ---- QMC band ----
    qmc_min = [5.1,4.7,4.4]
    qmc_max = [5.4,5.13,4.96]

    axins3.axhspan(qmc_min[idx], qmc_max[idx], color="gold", alpha=0.65)

    axins3.set_xticks([])
    axins3.set_ylabel(r"$E_f$ (eV)", fontsize=11)
    axins3.tick_params(axis='y', labelsize=9)
    



   ############################################
    # Right panel : Delta Ingredients
    ############################################
    # --- core / semicore shading ---
    color_base = "orange"
    ax_right.axvspan(x1[0], x1[1], color=color_base, alpha=0.40, zorder=0)
    ax_right.axvspan(x2[0], x2[1], color=color_base, alpha=0.15, zorder=0)
    ax_right.axvspan(x3[0], x3[1], color="#E65100", alpha=0.20, zorder=0)
    ax_right.axvspan(x4[0], x4[1], color=color_base, alpha=0.20, zorder=0)

    ax_right.plot(x, ds,  color='tab:blue',  lw=lwidth, label=r"$\Delta s$")
    ax_right.plot(x, da,  color='tab:green', lw=lwidth, label=r"$\Delta \alpha$")
    ax_right.plot(x, drs, color='tab:red',   lw=lwidth, label=r"$\Delta r_s$")

    ax_right.axhline(0, color="k", lw=0.8, ls="--")

    ax_right.set_ylim(*ylim_right)
    ax_right.set_xlim(*xlim_right)
    
    ax_right.legend(frameon=True, fontsize=14)
    
    ax_right.legend(frameon=True,
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        fontsize=14)
    
    ax_right.set_xlabel(r"$r (\mathrm{\AA})$", fontsize=20)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(f"{metal_name}_{defect_name}.pdf", dpi=300, bbox_inches="tight")
    plt.show()
    print(f"{metal_name}_{defect_name}.pdf")
    plt.close()
    

#plot_limits = {
#    "T": {"ylim_left": (-0.0070, 0.020), "ylim_right": (-2.3, 1.0), "xlim_left": (0.3, 2.5), "xlim_right": (0.51, 2.4), "Ef_range": (0.0, 5.0),  "reg_x1":(0.0, 0.16), "reg_x2":(0.16, 0.82), "reg_x3":(0.82, 1.33),"reg_x4":(1.33, 2.53)}
    

plot_limits = {
    "T": {"ylim_left": (-0.0070, 0.020), "ylim_right": (-2.3, 1.0), "xlim_left": (0.3, 2.5), "xlim_right": (0.51, 2.4), "Ef_range": (0.0, 5.0),  "reg_x1":(0.0, 0.16), "reg_x2":(0.16, 0.50), "reg_x3":(0.50, 2.21),"reg_x4":(2.21, 2.53)}
#    "X": {"ylim_left": (-0.005, 0.002), "ylim_right": (-0.11, 0.15), "xlim_left": (0.3, 2.5), "xlim_right": (0.51, 2.8), "Ef_range": (0.0, 5.0),"reg_x1":(0.0, 0.20), "reg_x2":(0.20, 0.48)},
#    "H": {"ylim_left": (-0.013, 0.017), "ylim_right": (-8.5, 0.6), "xlim_left": (0.3, 2.5), "xlim_right": (0.51, 2.8), "Ef_range": (0.0, 5.0),"reg_x1":(0.0, 0.29), "reg_x2":(0.29, 0.62)}
}

metal = "Si"
#DEFECTS = ["T", "X", "H"]
DEFECTS = ["T"]
for defect in DEFECTS:
    plot_defect(metal, defect, plot_limits[defect])


#
#for metal in plot_limits.keys():
#    for defect in DEFECTS:
#        plot_defect(metal, defect, plot_limits[defect])
        
        
#main_metals = ["Al","Ni","Pd","Pt"]
#supp_metals = ["Cu","Ag","Au","Pb"]
#plot_group(main_metals, "main_metals.pdf")
#plot_group(supp_metals, "supplement_metals.pdf")
