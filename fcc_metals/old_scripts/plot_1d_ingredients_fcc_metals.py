#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
from matplotlib.lines import Line2D


# =========  Global style controls  =========
BASE_FONT = 24
LEGENDS = 22

plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": BASE_FONT,
    "axes.linewidth": 1.2,
    "xtick.labelsize": BASE_FONT,
    "ytick.labelsize": BASE_FONT,
})


# ======================================================
# Master plotting function
# ======================================================
def plot_metal(metal_name, limits, ax_left, ax_mid, ax_right):
    ylim_mid = limits["ylim_mid"]
    ylim_bot = limits["ylim_bot"]
    Ef_range = limits["Ef_range"]
    x1 = limits["reg_x1"]
    x2 = limits["reg_x2"]
    x3 = limits["reg_x3"]
    x4 = limits["reg_x4"]

    # ---------------------------------------------------
    # Load data
    # ---------------------------------------------------
    
    base_path = f"/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_May19_metals_only"
    #file_pr = f"{base_path}/{metal_name}/1d_data_path_ae.dat"
    file_pr = f"{base_path}/{metal_name}/1d_data_path_ae.dat"
    
    #file_def = f"{base_path}/{metal_name}v/1d_data_path_ae.dat"
    file_def = f"{base_path}/{metal_name}v/1d_data_path_ae.dat"

    data_pr = np.loadtxt(file_pr, comments="#")
    data_def = np.loadtxt(file_def, comments="#")

    #x = data_pr[:, 1]
    #x = data_pr[:, 1]
    #y = data_pr[:, 2]
    y = data_pr[:, 1] # this is actually x
    n_pr, n_def = data_pr[:, 4], data_def[:, 4]
    #ones_pr, ones_def = np.ones_like(x), np.ones_like(x)
    ones_pr, ones_def = np.ones_like(y), np.ones_like(y)


    Fx_pr = [data_pr[:, 9], data_pr[:, 10], data_pr[:, 11], data_pr[:, 12], ones_pr]
    Fc_pr = [data_pr[:, 13], data_pr[:, 14], data_pr[:, 15], data_pr[:, 16], data_pr[:, 17]]
    Fx_def = [data_def[:, 9], data_def[:, 10], data_def[:, 11], data_def[:, 12], ones_def]
    #Fc_def = [data_def[:, 13], data_def[:, 14], data_def[:, 15], data_def[:, 16], data_pr[:, 17]] # had a typo in data_pr[:, 17]], it should be data_def[:, 17]]
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

    ############################################
    # Left plot: Delta Ratios
    ############################################
    # --- core / semicore shading ---
    color_base = "orange"
#    ax_left.axvspan(x1[0], x1[1], color=color_base, alpha=0.40, zorder=0)
#    ax_left.axvspan(x2[0], x2[1], color=color_base, alpha=0.20, zorder=0)
#    ax_left.axvspan(x3[0], x3[1], color="#E65100", alpha=0.20, zorder=0)
#    ax_left.axvspan(x4[0], x4[1], color=color_base, alpha=0.20, zorder=0)
    

    reference = ratio_lda
    lwidth=3;
    
    ax_left.plot(y, -(ratio_pbe    - reference), color="tab:blue",   lw=lwidth, ls="-", label="PBE")
    ax_left.plot(y, -(ratio_scan   - reference), color="tab:orange", lw=lwidth, ls="-", label="SCAN")
    ax_left.plot(y, -(ratio_r2scan - reference), color="tab:green",  lw=lwidth, ls="-", label="r2SCAN")
    ax_left.plot(y, -(ratio_lak    - reference), color="tab:cyan",   lw=lwidth, ls="-", label="LAK")
    ax_left.plot(y, -(ratio_lda    - reference), color="tab:purple", lw=lwidth, ls="-", label="LDA")
    
    ax_left.axhline(0.0, color="k", lw=0.8, ls="--")
    ax_left.set_ylabel(r"$\Delta R_{xc}$", fontsize=24)
    ax_left.set_xlabel(r"$(\AA)$", fontsize=20)
#    ax_left.set_ylim(-0.0030, 0.0030)
#    ax_left.set_xlim(2.2,4.00)
      
    # ===== Inset for formation energies =====
    df_main = pd.read_excel("convergence_report_experimental_formation_energies.xlsx", sheet_name="table_plot")
    row = df_main[df_main["Solids"] == metal_name].iloc[0]
    functionals = ["LDA", "PBE", "SCAN", "r2SCAN", "LAK"]
    colors_func = ["tab:purple", "tab:blue", "tab:orange", "tab:green", "tab:cyan"]
    energies = [row[f] for f in functionals]
    rec_exp = row["Recommended Expt"]

    axins3 = inset_axes(ax_left, width="20%", height="30%", loc="lower right")
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
    # Mid panel : Delta Ingredients
    ############################################
    color_base = "orange"
#    ax_mid.axvspan(x1[0], x1[1], color=color_base, alpha=0.40, zorder=0)
#    ax_mid.axvspan(x2[0], x2[1], color=color_base, alpha=0.15, zorder=0)
#    ax_mid.axvspan(x3[0], x3[1], color="#E65100", alpha=0.20, zorder=0)
#    ax_mid.axvspan(x4[0], x4[1], color=color_base, alpha=0.20, zorder=0)
    ax_mid.plot(y, 3*ds,  color='tab:blue',  lw=lwidth, label=r"$3\Delta s$")
    ax_mid.plot(y, da,  color='tab:green', lw=lwidth, label=r"$\Delta \alpha$")
    ax_mid.plot(y, drs, color='tab:red',   lw=lwidth, label=r"$\Delta r_s$")
    ax_mid.axhline(0, color="k", lw=0.8, ls="--")
    #ax_mid.set_ylim(*ylim_mid)
    #ax_mid.set_xlim(0.30, 1.60)
    
    ############################################
    # Right panel : Pure Ingredients
    ############################################
    color_base = "orange"
#    ax_right.axvspan(x1[0], x1[1], color=color_base, alpha=0.40, zorder=0)
#    ax_right.axvspan(x2[0], x2[1], color=color_base, alpha=0.20, zorder=0)
#    ax_right.axvspan(x3[0], x3[1], color="#E65100", alpha=0.20, zorder=0)
#    ax_right.axvspan(x4[0], x4[1], color=color_base, alpha=0.20, zorder=0)
    
    ax_right.plot(y, s_pr,  color='tab:blue',  lw=lwidth, ls="--", label=r"$s_{\mathrm{pris}}$")
    ax_right.plot(y, a_pr,  color='tab:green', lw=lwidth, ls="--", label=r"$\alpha_{\mathrm{pris}}$")
    ax_right.plot(y, rs_pr, color='tab:red',   lw=lwidth, ls="--", label=r"$r_{s,\mathrm{pris}}$")
    
    ax_right.plot(y, s_def,  color='tab:blue',  lw=lwidth, label=r"$s_{\mathrm{def}}$")
    ax_right.plot(y, a_def,  color='tab:green', lw=lwidth, label=r"$\alpha_{\mathrm{def}}$")
    ax_right.plot(y, rs_def, color='tab:red',   lw=lwidth, label=r"$r_{s,\mathrm{def}}$")

    ax_right.axhline(0, color="k", lw=0.8, ls="--")
    #ax_right.set_ylim(*ylim_mid)
    #ax_right.set_xlim(0.30, 1.7)
    ax_right.axhline(1.0, color="k", ls="--", lw=1) # plot one horizontal line at y=1 in right plot
#    reference = ratio_lda
#    lwidth=3;
#    ax_right.plot(x, ratio_pbe, color="tab:blue",   lw=lwidth, ls="-", label="PBE")
#    ax_right.plot(x, ratio_scan, color="tab:orange", lw=lwidth, ls="-", label="SCAN")
#    ax_right.plot(x, ratio_r2scan, color="tab:green",  lw=lwidth, ls="-", label="r2SCAN")
#    ax_right.plot(x, ratio_lak, color="tab:cyan",   lw=lwidth, ls="-", label="LAK")
#    ax_right.plot(x, ratio_lda, color="tab:purple", lw=lwidth, ls="-", label="LDA")
#
#    ax_right.axhline(0.0, color="k", lw=0.8, ls="--")
#    #ax_right.set_ylabel(r"$R_{xc}$", fontsize=20)
#    #ax_right.set_xlabel(r"$(\AA)$", fontsize=20)
#    ax_right.set_ylim(-0.030, 0.04)
#    ax_right.set_xlim(0.1,1.75)

 
def plot_group(metals, filename):
    
    fig, axes = plt.subplots(
    nrows=4,
    ncols=3,
    #figsize=(16, 14),   # wider figure
    figsize=(14, 18),
    sharex='col',
    gridspec_kw={"width_ratios":[1,1,1]}
    )

    for i, metal in enumerate(metals):

        ax_left  = axes[i,0]
        ax_mid = axes[i,1]
        ax_right = axes[i,2]
        
#        plot_metal(
#            metal,
#            plot_limits[metal],
#            ax_left,
#            ax_mid
#        )
        plot_metal(metal, plot_limits[metal], ax_left, ax_mid, ax_right)
        
        if i == 0:
            handles_left, labels_left = ax_left.get_legend_handles_labels()

            rec_exp_legend = Line2D(
                [0], [0],
                color='crimson',
                lw=2.0,
                ls='--',
                label="Exp. Value"
            )

            handles_left = handles_left + [rec_exp_legend]
            labels_left  = labels_left  + ["Exp. Value"]
            handles_right, labels_right = ax_mid.get_legend_handles_labels()
            handles_third, labels_third = ax_right.get_legend_handles_labels()
        
            fig.legend(
                handles_left, labels_left,
                loc="upper center",
                bbox_to_anchor=(0.20, 1.0),
                ncol=2,
                frameon=True,
                fontsize=LEGENDS
            )

            fig.legend(
                handles_right, labels_right,
                loc="upper center",
                bbox_to_anchor=(0.53, 1.0),
                ncol=1,
                frameon=True,
                fontsize=LEGENDS
            )
            
            fig.legend(
            handles_third, labels_third,
            loc="upper center",
            bbox_to_anchor=(0.85, 1.0),
            ncol=2,
            fontsize=LEGENDS
            )


        ax_left.set_title(metal, fontsize=42, loc="left", y=0.9, pad=-14)

    axes[-1,0].set_xlabel(r"$r (\AA)$")
    axes[-1,1].set_xlabel(r"$r (\AA)$")
    axes[-1,2].set_xlabel(r"$r (\AA)$")

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    #plt.tight_layout()
    #plt.tight_layout(rect=[0,0,1,0.96])
    plt.savefig(filename, dpi=300)
    #plt.show()
    #plt.close()
    
    
    for ax in axes[:,0]:
        ax.set_ylabel(r"$\Delta R_{xc}$")

    for ax in axes[:,1]:
        ax.set_ylabel(r"$\Delta$ Ingredients")
        
    for ax in axes[:,2]:
        ax.set_ylabel(r"$R_{xc}$")



plot_limits = {
    "Al": {"ylim_mid": (-1.5, 4.5), "ylim_bot": (0.0, 4.5),"Ef_range": (0.0, 1.3), "reg_x1":(0.0, 0.16), "reg_x2":(0.16, 0.42),"reg_x3":(0.42, 1.00),"reg_x4":(1.00, 1.58)},
    "Ni": {"ylim_mid": (-1.5, 4.5), "ylim_bot": (0.0, 4.5), "Ef_range": (0.0, 2.3),"reg_x1":(0.0, 0.20), "reg_x2":(0.20, 0.48), "reg_x3":(0.48, 0.88),"reg_x4":(0.88, 1.25)},
    "Pd": {"ylim_mid": (-1.5, 4.5), "ylim_bot": (0.0, 4.5), "Ef_range": (0.0, 2.0),"reg_x1":(0.0, 0.29), "reg_x2":(0.29, 0.62),"reg_x3":(0.62, 0.97),"reg_x4":(0.97, 1.3)},
    "Pt": {"ylim_mid": (-1.5, 4.5), "ylim_bot": (0.0, 4.5), "Ef_range": (0.0, 1.5),"reg_x1":(0.0, 0.34), "reg_x2":(0.34, 0.63),"reg_x3":(0.63, 0.98),"reg_x4":(0.98, 1.32)},
    ## Supplementary metals below
    "Cu": {"ylim_mid": (-1.5, 4.5), "ylim_bot": (0.0, 4.5), "Ef_range": (0.0, 1.6),"reg_x1":(0.0, 0.18), "reg_x2":(0.18, 0.51),"reg_x3":(0.51, 0.89),"reg_x4":(0.89, 1.28)},
    "Ag": {"ylim_mid": (-1.5, 4.5), "ylim_bot": (0.0, 4.5), "Ef_range": (0.0, 1.30),"reg_x1":(0.0, 0.30), "reg_x2":(0.30, 0.68),"reg_x3":(0.68, 1.02),"reg_x4":(1.02, 1.36)},
    "Au": {"ylim_mid": (-1.5, 4.5), "ylim_bot": (0.0, 4.5), "Ef_range": (-0.4, 1.0),"reg_x1":(0.0, 0.35), "reg_x2":(0.35, 0.66),"reg_x3":(0.66, 1.02),"reg_x4":(1.02, 1.36)},
    "Pb": {"ylim_mid": (-1.5, 4.5), "ylim_bot": (0.0, 4.5), "Ef_range": (-0.3, 0.65),"reg_x1":(0.0, 0.32), "reg_x2":(0.32, 0.83),"reg_x3":(0.83, 1.23),"reg_x4":(1.23, 1.74)},
}


main_metals = ["Al","Ni","Pd","Pt"]
#supp_metals = ["Cu","Ag","Au","Pb"]

plot_group(main_metals, "main_metals.pdf")
#plot_group(supp_metals, "supplement_metals.pdf")
