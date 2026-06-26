#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt


### Script to plot the ratio 1/R_{xc} for Si-diamond interstitials
### Make sure you are in the same folder with the folders with the generated ingredients.

### Inputs: .dat files present in all folders for all metals pristine and defective, T_def_image.png (which has no background) to be added in the background of the plots.
### Outputs: a PDF(s) with a 1x4 grid with the inverse ratio \Delta(1/R_{xc}) for PBE, LAK, SCAN, and r2SCAN.



### Plotting function
def plot_2d_ingredients(defect):
    #=========================================================
    # Load data
    #=========================================================
    base_path = "/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_May19_silicon_only"
    
    
    file_pr = f"{base_path}/Si/data_2D.dat"
    file_def = f"{base_path}/{defect}/data_2D.dat"

    data_pr  = np.loadtxt(file_pr, comments="#")
    data_def = np.loadtxt(file_def, comments="#")

    # =========================================================
    # Extract columns
    # =========================================================

    x_pr = data_pr[:,1]
    y_pr = data_pr[:,2]

    x_def = data_def[:,1]
    y_def = data_def[:,2]
    
    n_pr, n_def = data_pr[:, 4], data_def[:, 4]
    ones_pr, ones_def = np.ones_like(n_pr), np.ones_like(n_def)

    rs_pr    = data_pr[:,5]
    s_pr     = data_pr[:,7]
    alpha_pr = data_pr[:,8]

    rs_def    = data_def[:,5]
    s_def     = data_def[:,7]
    alpha_def = data_def[:,8]
    
    # Differences of the ingredients
    ds  = s_def  - s_pr
    da  = alpha_def  - alpha_pr
    drs = rs_def - rs_pr
    print("drs", drs)
    
    # Exractions of Fxc
    Fx_pr = [data_pr[:, 9], data_pr[:, 10], data_pr[:, 11], data_pr[:, 12], ones_pr]
    Fc_pr = [data_pr[:, 13], data_pr[:, 14], data_pr[:, 15], data_pr[:, 16], data_pr[:, 17]]
    Fx_def = [data_def[:, 9], data_def[:, 10], data_def[:, 11], data_def[:, 12], ones_def]
    Fc_def = [data_def[:, 13], data_def[:, 14], data_def[:, 15], data_def[:, 16], data_def[:, 17]]
    
    # xc energy density
    exc_pr = [n_pr**(4/3) * (Fx + Fc) for Fx, Fc in zip(Fx_pr, Fc_pr)]
    exc_def = [n_def**(4/3) * (Fx + Fc) for Fx, Fc in zip(Fx_def, Fc_def)]
    
    
   # Ratios, for Si is inverse: p/d, i.e., pristine over defective
   # ratio_pbe, ratio_scan, ratio_r2scan, ratio_lak, ratio_lda = [d/p for d, p in zip(exc_def, exc_pr)]
    
    ratio_pbe, ratio_scan, ratio_r2scan, ratio_lak, ratio_lda = [p/d for d, p in zip(exc_def, exc_pr)]
    
    reference = ratio_lda
    r_pbe = -(ratio_pbe - reference)
    r_scan = -(ratio_scan - reference)
    r_r2scan = -(ratio_r2scan - reference)
    r_lak = -(ratio_lak - reference)
    
    # =========================================================
    # Grid dimensions
    # =========================================================

    NGX = 160
    NGY = 160

    # =========================================================
    # Reshape into 2D grids
    # =========================================================

    X_pr = x_pr.reshape(NGY, NGX)
    Y_pr = y_pr.reshape(NGY, NGX)

    X_def = x_def.reshape(NGY, NGX)
    Y_def = y_def.reshape(NGY, NGX)

    RS_pr = rs_pr.reshape(NGY, NGX)
    S_pr  = s_pr.reshape(NGY, NGX)
    A_pr  = alpha_pr.reshape(NGY, NGX)

    RS_def = rs_def.reshape(NGY, NGX)
    S_def  = s_def.reshape(NGY, NGX)
    A_def  = alpha_def.reshape(NGY, NGX)
    
    DS = ds.reshape(NGY, NGX)
    DA = da.reshape(NGY, NGX)
    DRS = drs.reshape(NGY, NGX)
    
    R_PBE = r_pbe.reshape(NGY, NGX)
    R_SCAN = r_scan.reshape(NGY, NGX)
    R_R2SCAN = r_r2scan.reshape(NGY, NGX)
    R_LAK = r_lak.reshape(NGY, NGX)

    # =========================================================
    # Figure
    # =========================================================

    fig, axes = plt.subplots(4, 1, figsize=(10, 12))
    
    fig.suptitle(
    r"$\Delta(1/R_{xc})$",
    fontsize=20,
    x=0.75,
    y=0.92
    )
    #=\left(\frac{e_{xc}^{pris}}{e_{xc}^{def}}\right)^{LDA} -  \left(\frac{e_{xc}^{pris}}{e_{xc}^{def}}\right)^{DFA}
    
    xmin = X_pr.min()
    xmax = X_pr.max()

    ymin = Y_pr.min()
    ymax = Y_pr.max()
    
    img = plt.imread("T_def_image.png")
    for i in range(0,4):
        axes[i].imshow(
        img,
        extent=[xmin,xmax,ymin,ymax],
        alpha=0.6,
        zorder=2
        )

    # =========================================================
    # r_s
    # =========================================================
    z_min = -0.02
    z_max = 0.02
    z_ord = 1
    
    levels_rs = np.linspace(z_min, z_max, 100)
    im1 = axes[0].contourf(
    X_pr, Y_pr, R_PBE,
    levels=levels_rs,
    cmap='viridis',
    extend='max',
    zorder=z_ord
    )
    im1.set_rasterized(True)
    
    levels_s = np.linspace(z_min, z_max, 100)
    im2 = axes[1].contourf(
        X_pr, Y_pr, R_LAK,
        levels=levels_rs,
        cmap='viridis',
        extend='max',
        zorder=z_ord
    )
    im2.set_rasterized(True)


    levels_a = np.linspace(z_min, z_max, 100)
    im3 = axes[2].contourf(
    X_pr, Y_pr, R_SCAN,
    levels=levels_rs,
    cmap='viridis',
    extend='max',
    zorder=z_ord
    )
    im3.set_rasterized(True)
    
    levels_a = np.linspace(z_min, z_max, 100)
    im4 = axes[3].contourf(
    X_pr, Y_pr, R_R2SCAN,
    levels=levels_rs,
    cmap='viridis',
    extend='max',
    zorder=z_ord
    )
    im4.set_rasterized(True)
    
    # =========================================================
    # Titles
    # =========================================================
    #axes.set_title(r"$\Delta (1/R_{xc})$", fontsize=18)
#    axes[0].set_title(r"$\delta r_s$", fontsize=18)
#    axes[1].set_title(r"$\delta s$", fontsize=18)
#    axes[2].set_title(r"$\delta \alpha$", fontsize=18)


    # =========================================================
    # Labels and formatting
    # =========================================================

    #row_labels = [r"$r_s$", r"$s$", r"$\alpha$"]
    row_labels = ["PBE", "LAK", "SCAN", "R2SCAN"]

    for i in range(4):
        axes[i].set_ylabel("y (Å)", fontsize=14)
        axes[i].set_xlabel("x (Å)", fontsize=14)
        #axes[i].axis("equal")
        axes[i].set_aspect('equal')

        # ingredient label on left side
        axes[i].text(
            -0.55, 0.5,
            row_labels[i],
            transform=axes[i].transAxes,
            fontsize=20,
            rotation=90,
            va="center"
        )

    # =========================================================
    # Colorbars
    # =========================================================

    cbar1 = fig.colorbar(
        im1,
        ax=axes[0],
        location="right",
        fraction=0.046,
        pad=0.04
    )

    #cbar1.set_label(r"$r_s$", fontsize=16)


    cbar2 = fig.colorbar(
        im2,
        ax=axes[1],
        location="right",
        fraction=0.046,
        pad=0.04
    )

    #cbar2.set_label(r"$s$", fontsize=16)


    cbar3 = fig.colorbar(
        im3,
        ax=axes[2],
        location="right",
        fraction=0.046,
        pad=0.04
    )

    #cbar3.set_label(r"$\alpha$", fontsize=16)
    
    cbar4 = fig.colorbar(
        im4,
        ax=axes[3],
        location="right",
        fraction=0.046,
        pad=0.04
    )

    #cbar3.set_label(r"$\alpha$", fontsize=16)

        # =========================================================
    # Layout + save
    # =========================================================

    plt.subplots_adjust(left=0.15)

    #plt.tight_layout()
    #
    #plt.savefig("Al_2d_ingredients.pdf", dpi=300)
    plt.subplots_adjust(
        left=0.15,
        right=0.88,
        hspace=0.4,
        wspace=0.15
    )

    plt.savefig(
        f"{defect}_delta_Rxc_2d_ingredients.pdf",
        dpi=600,
        bbox_inches="tight"
    )

    #plt.show()
    plt.close()

#defects = ["T", "X", "H"]
#defects = ["T"]
defects = ["T"]
for defect in defects:
    print(f"Plotting {defect}...")
    plot_2d_ingredients(defect)
