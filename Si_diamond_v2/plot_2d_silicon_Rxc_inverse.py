#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# =========================================================
# Load data
# =========================================================

#base_path = "/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_May19_silicon_only"

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
    #N_pr = n_pr.reshape(NGY, NGX)
    #N_def = n_def.reshape(NGY, NGX)
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
    
    #ratio_pbe, ratio_scan, ratio_r2scan, ratio_lak, ratio_lda = [p/d for d, p in zip(exc_def, exc_pr)]

    eps = 1e-12
    ratio_pbe, ratio_scan, ratio_r2scan, ratio_lak, ratio_lda = [
        p / (d + eps) for d, p in zip(exc_def, exc_pr)]
    
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
    
    R_PBE = ratio_pbe.reshape(NGY, NGX)
    R_SCAN = ratio_scan.reshape(NGY, NGX)
    R_R2SCAN = ratio_r2scan.reshape(NGY, NGX)
    R_LAK = ratio_lak.reshape(NGY, NGX)
    R_LDA = ratio_lda.reshape(NGY, NGX)
    reference = R_LDA
    
#    print(R_PBE)
   
    
    # =========================================================
    # Figure
    # =========================================================

    fig, axes = plt.subplots(3, 1, figsize=(8, 10))

    # =========================================================
    # r_s
    # =========================================================
    z_min = -1
    z_max = 1
    
    levels_rs = np.linspace(z_min, z_max, 100)
    im1 = axes[0].contourf(
    X_pr, Y_pr, R_LDA,
    levels=levels_rs,
    cmap='viridis',
    extend='max'
    )
    im1.set_rasterized(True)
    
    levels_s = np.linspace(z_min, z_max, 100)
    im2 = axes[1].contourf(
        X_pr, Y_pr, R_LAK,
        levels=levels_s,
        cmap='viridis',
        extend='max'
    )
    im2.set_rasterized(True)

    # =========================================================
    # s
    # =========================================================

    levels_a = np.linspace(z_min, z_max, 100)
    im3 = axes[2].contourf(
    X_pr, Y_pr, R_PBE,
    levels=levels_a,
    cmap='viridis',
    extend='max'
    )
    im3.set_rasterized(True)
    
    # =========================================================
    # Titles
    # =========================================================

#    axes[0].set_title(r"$\delta r_s$", fontsize=18)
#    axes[1].set_title(r"$\delta s$", fontsize=18)
#    axes[2].set_title(r"$\delta \alpha$", fontsize=18)


    # =========================================================
    # Labels and formatting
    # =========================================================

    #row_labels = [r"$r_s$", r"$s$", r"$\alpha$"]
    row_labels = ["LDA", "LAK", "PBE"]

    for i in range(3):
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
        f"{defect}_ratios_2d_ingredients.pdf",
        dpi=600,
        bbox_inches="tight"
    )

    #plt.show()
    plt.close()

#defects = ["Al","Ni","Pd","Pt","Cu","Ag","Au","Pb"]
#defects = ["T", "X", "H"]
#defects = ["T"]
defects = ["T"]
for defect in defects:
    print(f"Plotting {defect}...")
    plot_2d_ingredients(defect)
