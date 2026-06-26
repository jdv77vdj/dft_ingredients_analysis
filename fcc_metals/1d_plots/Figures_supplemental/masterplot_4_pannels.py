#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
from matplotlib.lines import Line2D


### Script to generate 1D ingredients, the original version for the fcc metals only.
### Example for executing it: python masterplot_4_pannels.py

### Required input files: data.dat file, written for all systems in their respective folder.
### Outputs files: a folder with a PDF grid plot for each metal with the ingredients.


# =========  Global style controls  =========
BASE_FONT = 26

plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": BASE_FONT,
    "axes.linewidth": 1.2,
    "xtick.labelsize": BASE_FONT,
    "ytick.labelsize": BASE_FONT,
    "legend.fontsize": BASE_FONT,
})

# ======================================================
# Master plotting function
# ======================================================
def plot_metal(metal_name, limits, ax_left_top, ax_left_bot, ax_right_top, ax_right_bot):
    ylim_right_top = limits["ylim_right_top"]
    ylim_bot = limits["ylim_bot"]
    Ef_range = limits["Ef_range"]
    x1 = limits["reg_x1"]
    x2 = limits["reg_x2"]
    x3 = limits["reg_x3"]
    x4 = limits["reg_x4"]
    bbox_right_bot = limits["bbox_right_bot"]
    inset_loc = limits["inset_loc"]
    
    x_limit = 41

    # ---------------------------------------------------
    # Load data
    # ---------------------------------------------------
    
    base_path = f"/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_March6"
    file_pr = f"{base_path}/{metal_name}/data_path_ae.dat"
    file_def = f"{base_path}/{metal_name}v/data_path_ae.dat"

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
    ratio_pbe, ratio_scan, ratio_r2scan, ratio_lak, ratio_lda = [d/p for d, p in zip(exc_def, exc_pr)]
    rs_pr, s_pr, a_pr = data_pr[:, 5], data_pr[:, 7], data_pr[:, 8]
    rs_def, s_def, a_def = data_def[:, 5], data_def[:, 7], data_def[:, 8]

    # Differences of the ingredients
    ds  = s_def  - s_pr
    da  = a_def  - a_pr
    drs = rs_def - rs_pr
    
    dR_pbe    = np.gradient(ratio_pbe, x)
    dR_scan   = np.gradient(ratio_scan, x)
    dR_r2scan = np.gradient(ratio_r2scan, x)
    dR_lak    = np.gradient(ratio_lak, x)
    dR_lda    = np.gradient(ratio_lda, x)

    ############################################
    # Left Top plot: Delta Ratios
    ############################################
    # --- core / semicore shading ---
    color_base = "orange"
    ax_left_top.axvspan(x1[0], x1[1], color=color_base, alpha=0.40, zorder=0)
    ax_left_top.axvspan(x2[0], x2[1], color=color_base, alpha=0.20, zorder=0)
    ax_left_top.axvspan(x3[0], x3[1], color="#E65100", alpha=0.20, zorder=0)
    ax_left_top.axvspan(x4[0], x4[1], color=color_base, alpha=0.20, zorder=0)

    reference = ratio_lda
    lwidth=3;
    
    ax_left_top.plot(x, -(ratio_pbe    - reference), color="tab:blue",   lw=lwidth, ls="-", label="PBE")
    ax_left_top.plot(x, -(ratio_scan   - reference), color="tab:orange", lw=lwidth, ls="-", label="SCAN")
    ax_left_top.plot(x, -(ratio_r2scan - reference), color="tab:green",  lw=lwidth, ls="-", label="r2SCAN")
    ax_left_top.plot(x, -(ratio_lak    - reference), color="tab:cyan",   lw=lwidth, ls="-", label="LAK")
    ax_left_top.plot(x, -(ratio_lda    - reference), color="tab:purple", lw=lwidth, ls="-", label="LDA")
    
    ax_left_top.axhline(0.0, color="k", lw=0.8, ls="--")
    ax_left_top.set_ylabel(r"$\Delta R_{xc}$", fontsize=28)
    #ax_left.set_xlabel(r"$(\AA)$", fontsize=20)
    ax_left_top.set_ylim(-0.030, 0.05)
    ax_left_top.set_xlim(0.1,1.75)
      
    # ===== Inset for formation energies =====
    df_main = pd.read_excel("convergence_report_experimental_formation_energies.xlsx", sheet_name="table_plot")
    row = df_main[df_main["Solids"] == metal_name].iloc[0]
    functionals = ["LDA", "PBE", "SCAN", "r2SCAN", "LAK"]
    colors_func = ["tab:purple", "tab:blue", "tab:orange", "tab:green", "tab:cyan"]
    energies = [row[f] for f in functionals]
    rec_exp = row["Recommended Expt"]

    #axins3 = inset_axes(ax_left_top, width="20%", height="30%", loc=inset_loc)
    inset_pos = limits.get("inset_pos", (0.75, 0.05))
    axins3 = inset_axes(
        ax_left_top,
        width="20%",
        height="30%",
        loc="lower left",
        bbox_to_anchor=(*inset_pos, 1, 1),
        bbox_transform=ax_left_top.transAxes,
        borderpad=0
        )
    
    x_f = np.arange(len(functionals))
    axins3.bar(x_f, energies, color=colors_func, alpha=1, width=0.7)
    axins3.set_zorder(200)
    axins3.axhline(0, color="k", lw=1.2) # zero reference line
    if not pd.isna(rec_exp):
        axins3.axhline(rec_exp, color="crimson", lw=2.0, ls="--")
    axins3.set_ylabel(r"$E_{f}$ (eV)", fontsize=15)
    #axins3.set_ylim(*limitsEf["Ef_range"])
    axins3.set_ylim(*Ef_range)
    axins3.set_xticks([])
    axins3.tick_params(axis='both', labelsize=15)
  

   ############################################
    # Right top : Delta Ingredients
    ############################################
    color_base = "orange"
    ax_right_top.axvspan(x1[0], x1[1], color=color_base, alpha=0.40, zorder=0)
    ax_right_top.axvspan(x2[0], x2[1], color=color_base, alpha=0.15, zorder=0)
    ax_right_top.axvspan(x3[0], x3[1], color="#E65100", alpha=0.20, zorder=0)
    ax_right_top.axvspan(x4[0], x4[1], color=color_base, alpha=0.20, zorder=0)
    ax_right_top.plot(x, 2*ds,  color='tab:blue',  lw=lwidth, label=r"$2\Delta s$")
    ax_right_top.plot(x, da,  color='tab:green', lw=lwidth, label=r"$\Delta \alpha$")
    ax_right_top.plot(x, drs, color='tab:red',   lw=lwidth, label=r"$\Delta r_s$")
    ax_right_top.axhline(0, color="k", lw=0.8, ls="--")
    #ax_right_top.set_ylim(*ylim_right_top) if using custom limits
    ax_right_top.set_ylim(-1.0, 2.5)
    ax_right_top.set_xlim(0.30, 1.50)
    
    
    
    ############################################
    # Left bottom panel : Pure Ratios R_{xc}
    ############################################
    # --- core / semicore shading ---
    color_base = "orange"
    ax_left_bot.axvspan(x1[0], x1[1], color=color_base, alpha=0.40, zorder=0)
    ax_left_bot.axvspan(x2[0], x2[1], color=color_base, alpha=0.20, zorder=0)
    ax_left_bot.axvspan(x3[0], x3[1], color="#E65100", alpha=0.20, zorder=0)
    ax_left_bot.axvspan(x4[0], x4[1], color=color_base, alpha=0.20, zorder=0)
    

    lwidth=3;
    
    ax_left_bot.plot(x, ratio_pbe, color="tab:blue",   lw=lwidth, ls="-", label="PBE")
    ax_left_bot.plot(x, ratio_scan, color="tab:orange", lw=lwidth, ls="-", label="SCAN")
    ax_left_bot.plot(x, ratio_r2scan , color="tab:green",  lw=lwidth, ls="-", label="r2SCAN")
    ax_left_bot.plot(x, ratio_lak, color="tab:cyan",   lw=lwidth, ls="-", label="LAK")
    ax_left_bot.plot(x, ratio_lda, color="tab:purple", lw=lwidth, ls="-", label="LDA")
    
    ax_left_bot.set_ylabel(r"$R_{xc}$", fontsize=28)
    
    
    ############################################
    # Right bottom panel : Pure Ingredients
    ############################################
    color_base = "orange"
    ax_right_bot.axvspan(x1[0], x1[1], color=color_base, alpha=0.40, zorder=0)
    ax_right_bot.axvspan(x2[0], x2[1], color=color_base, alpha=0.20, zorder=0)
    ax_right_bot.axvspan(x3[0], x3[1], color="#E65100", alpha=0.20, zorder=0)
    ax_right_bot.axvspan(x4[0], x4[1], color=color_base, alpha=0.20, zorder=0)
    
    ax_right_bot.plot(x, s_pr,  color='tab:blue',  lw=lwidth, ls="--", label=r"$s_{\mathrm{pris}}$")
    ax_right_bot.plot(x, a_pr,  color='tab:green', lw=lwidth, ls="--", label=r"$\alpha_{\mathrm{pris}}$")
    ax_right_bot.plot(x, rs_pr, color='tab:red',   lw=lwidth, ls="--", label=r"$r_{s,\mathrm{pris}}$")
    
    ax_right_bot.plot(x, s_def,  color='tab:blue',  lw=lwidth, label=r"$s_{\mathrm{def}}$")
    ax_right_bot.plot(x, a_def,  color='tab:green', lw=lwidth, label=r"$\alpha_{\mathrm{def}}$")
    ax_right_bot.plot(x, rs_def, color='tab:red',   lw=lwidth, label=r"$r_{s,\mathrm{def}}$")

    ax_right_bot.axhline(0, color="k", lw=0.8, ls="--")
    #ax_right_bot.set_ylim(0, 3)
    ax_right_bot.set_ylim(*ylim_bot)
    ax_right_bot.set_xlim(0.30, 1.50)
    

    ax_right_bot.axhline(1.0, color="k", ls="--", lw=1) # plot one horizontal line at y=1 in right plot


def plot_group(metals, filename):
    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(14, 16),
        sharex=True,
        gridspec_kw={"width_ratios": [1, 1.05]}
    )

    ax_left_top  = axes[0,0]
    ax_right_top = axes[0,1]
    ax_left_bot  = axes[1,0]
    ax_right_bot = axes[1,1]
    

    metal_name = metals
    plot_metal(
        metal,
        plot_limits[metal],
        ax_left_top,
        ax_left_bot,
        ax_right_top,
        ax_right_bot
    )
        
    handles_left, labels_left = ax_left_top.get_legend_handles_labels()
    handles_right, labels_right = ax_right_top.get_legend_handles_labels()
    handles_right_bot, labels_right_bot = ax_right_bot.get_legend_handles_labels()
    
        # Add experimental line
    rec_exp_legend = Line2D(
        [0], [0],
        color='crimson',
        lw=2.0,
        ls='--',
        label="Exp. Value"
    )

    handles_left += [rec_exp_legend]
    labels_left  += ["Exp. Value"]

    alpha_legends = 0.95
    bbox_right_bot = plot_limits[metal]["bbox_right_bot"]

    fig.legend(
        handles_left, labels_left,
        loc="upper center",
        bbox_to_anchor=(0.30, 0.94),
        ncol=2,
        frameon=True,
        framealpha=alpha_legends,
        fontsize=20
    )

    fig.legend(
        handles_right, labels_right,
        loc="upper center",
        bbox_to_anchor=(0.775, 0.94),
        ncol=3,
        frameon=True,
        framealpha=alpha_legends,
        fontsize=20
    )
    
    fig.legend(
        handles_right_bot, labels_right_bot,
        loc="upper center",
        #bbox_to_anchor=(0.78, 0.46),
        bbox_to_anchor=bbox_right_bot,
        ncol=2,
        frameon=True,
        framealpha=alpha_legends,
        fontsize=20
    )


    # --- Titles ---
    ax_left_top.set_title(metal_name, fontsize=48, loc="left")

    # --- Labels ---
    ax_left_bot.set_xlabel(r"$r\;\mathrm{in}\;\AA$")
    ax_right_bot.set_xlabel(r"$r\;\mathrm{in}\;\AA$")
    
 
    #ax_left_top.set_ylabel(r"$\Delta R_{xc}$")
    
    ax_right_top.set_ylabel(r"$\Delta$ Ingredients")
    ax_right_bot.set_ylabel(r"Ingredients")
    
    
    text_font_size = 28
    
    ax_left_top.text(0.02, 0.02, "(a)", transform=ax_left_top.transAxes,
                 fontsize=text_font_size, va='bottom', ha='left')

    ax_right_top.text(0.02, 0.02, "(b)", transform=ax_right_top.transAxes,
                      fontsize=text_font_size, va='bottom', ha='left')

    ax_left_bot.text(0.02, 0.02, "(c)", transform=ax_left_bot.transAxes,
                     fontsize=text_font_size, va='bottom', ha='left')

    ax_right_bot.text(0.02, 0.02, "(d)", transform=ax_right_bot.transAxes,
                      fontsize=text_font_size, va='bottom', ha='left')

    # --- Layout ---
    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.savefig(filename, dpi=300)
 
 

inset_x = 0.79
inset_y = 0.04

plot_limits = {
    "Al": {"ylim_right_top": (-1.5, 3.0), "ylim_bot": (0.0, 5.5),"Ef_range": (0.0, 1.3), "reg_x1":(0.0, 0.16), "reg_x2":(0.16, 0.42),"reg_x3":(0.42, 1.00),"reg_x4":(1.00, 1.58),"bbox_right_bot":(0.78, 0.49),"inset_loc": "lower right","inset_pos":(inset_x, inset_y)},
    "Ni": {"ylim_right_top": (-1.5, 3.0), "ylim_bot": (0.0, 3.0), "Ef_range": (0.0, 2.3),"reg_x1":(0.0, 0.20), "reg_x2":(0.20, 0.48), "reg_x3":(0.48, 0.88),"reg_x4":(0.88, 1.25),"bbox_right_bot":(0.78, 0.49),"inset_loc": "lower right","inset_pos":(inset_x, inset_y)},
    "Pd": {"ylim_right_top": (-1.5, 3.0), "ylim_bot": (0.0, 3.0), "Ef_range": (0.0, 2.0),"reg_x1":(0.0, 0.29), "reg_x2":(0.29, 0.62),"reg_x3":(0.62, 0.97),"reg_x4":(0.97, 1.3),"bbox_right_bot":(0.78, 0.49),"inset_loc": "lower right","inset_pos":(inset_x, inset_y)},
    "Pt": {"ylim_right_top": (-1.5, 3.0), "ylim_bot": (0.0, 3.0), "Ef_range": (0.0, 1.5),"reg_x1":(0.0, 0.34), "reg_x2":(0.34, 0.63),"reg_x3":(0.63, 0.98),"reg_x4":(0.98, 1.32),"bbox_right_bot":(0.78, 0.49),"inset_loc": "lower right","inset_pos":(inset_x, inset_y)},
    ## Supplementary metals below
    "Cu": {"ylim_right_top": (-1.5, 3.0), "ylim_bot": (0.0, 3.0), "Ef_range": (0.0, 1.6),"reg_x1":(0.0, 0.18), "reg_x2":(0.18, 0.51),"reg_x3":(0.51, 0.89),"reg_x4":(0.89, 1.28),"bbox_right_bot":(0.78, 0.49),"inset_loc": "lower right","inset_pos":(inset_x, inset_y)},
    "Ag": {"ylim_right_top": (-1.5, 3.0), "ylim_bot": (0.0, 3.0), "Ef_range": (0.0, 1.30),"reg_x1":(0.0, 0.30), "reg_x2":(0.30, 0.68),"reg_x3":(0.68, 1.02),"reg_x4":(1.02, 1.36),"bbox_right_bot":(0.78, 0.49),"inset_loc": "lower right","inset_pos":(inset_x, inset_y)},
    "Au": {"ylim_right_top": (-1.5, 3.0), "ylim_bot": (0.0, 3.0), "Ef_range": (-0.4, 1.0),"reg_x1":(0.0, 0.35), "reg_x2":(0.35, 0.66),"reg_x3":(0.66, 1.02),"reg_x4":(1.02, 1.36),"bbox_right_bot":(0.78, 0.49),"inset_loc": "center","inset_pos":(0.23, 0.45)},
    "Pb": {"ylim_right_top": (-1.5, 3.0), "ylim_bot": (0.0, 3.0), "Ef_range": (-0.3, 0.65),"reg_x1":(0.0, 0.32), "reg_x2":(0.32, 0.83),"reg_x3":(0.83, 1.23),"reg_x4":(1.23, 1.74),"bbox_right_bot":(0.78, 0.49),"inset_loc": "center","inset_pos":(0.23, 0.45)},
}


main_metals = ["Al","Ni","Pd","Pt"]
supp_metals = ["Cu","Ag","Au","Pb"]

for metal in main_metals:
    plot_group(metal, f"{metal}_4panel.pdf")
    
for metal in supp_metals:
    plot_group(metal, f"{metal}_4panel.pdf")
