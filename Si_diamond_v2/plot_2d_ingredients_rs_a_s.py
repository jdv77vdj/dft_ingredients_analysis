#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# =========================================================
# Load data
# =========================================================

base_path = "/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_April19_2d_ingredients"

### Plotting function
def plot_2d_ingredients(defect):
    #=========================================================
    # Load data
    #=========================================================
    base_path = "/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_April19_2d_ingredients/Silicon"
    
    
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

    rs_pr    = data_pr[:,5]
    s_pr     = data_pr[:,7]
    alpha_pr = data_pr[:,8]

    rs_def    = data_def[:,5]
    s_def     = data_def[:,7]
    alpha_def = data_def[:,8]

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

    # =========================================================
    # Shared normalization
    # =========================================================

#    vmin_rs = min(RS_pr.min(), RS_def.min())
#    vmax_rs = max(RS_pr.max(), RS_def.max())
#
#    vmin_s = min(S_pr.min(), S_def.min())
#    vmax_s = max(S_pr.max(), S_def.max())
#
#    vmin_a = min(A_pr.min(), A_def.min())
#    vmax_a = max(A_pr.max(), A_def.max())

    # =========================================================
    # Figure
    # =========================================================

    fig, axes = plt.subplots(3, 2, figsize=(14, 16))

    # =========================================================
    # r_s
    # =========================================================

#    im1 = axes[0,0].contourf(
#        X_pr, Y_pr, RS_pr,
#        levels=100,
#        vmin=vmin_rs,
#        vmax=vmax_rs
#    )
#    im1.set_rasterized(True)
#
#    im2 = axes[0,1].contourf(
#        X_def, Y_def, RS_def,
#        levels=100,
#        vmin=vmin_rs,
#        vmax=vmax_rs
#    )
#    im2.set_rasterized(True)

    levels_rs = np.linspace(0, 4.0, 100)
    im1 = axes[0,0].contourf(
    X_pr, Y_pr, RS_pr,
    levels=levels_rs,
    cmap='inferno',
    extend='max'
    )
    im1.set_rasterized(True)
    
    im2 = axes[0,1].contourf(
        X_def, Y_def, RS_def,
        levels=levels_rs,
        cmap='inferno',
        extend='max'
    )
    im2.set_rasterized(True)

    # =========================================================
    # s
    # =========================================================

#    im3 = axes[1,0].contourf(
#        X_pr, Y_pr, S_pr,
#        levels=100,
#        vmin=vmin_s,
#        vmax=vmax_s
#    )
#    im3.set_rasterized(True)
#
#    im4 = axes[1,1].contourf(
#        X_def, Y_def, S_def,
#        levels=100,
#        vmin=vmin_s,
#        vmax=vmax_s
#    )
#    im4.set_rasterized(True)

    levels_s = np.linspace(0, 1.5, 100)
    im3 = axes[1,0].contourf(
    X_pr, Y_pr, S_pr,
    levels=levels_s,
    cmap='magma',
    extend='max'
    )
    im3.set_rasterized(True)
    
    im4 = axes[1,1].contourf(
        X_def, Y_def, S_def,
        levels=levels_s,
        cmap='magma',
        extend='max'
    )
    im4.set_rasterized(True)

    # =========================================================
    # alpha
    # =========================================================

#    im5 = axes[2,0].contourf(
#        X_pr, Y_pr, A_pr,
#        levels=100,
#        vmin=vmin_a,
#        vmax=vmax_a
#    )
#    im5.set_rasterized(True)
#
#    im6 = axes[2,1].contourf(
#        X_def, Y_def, A_def,
#        levels=100,
#        vmin=vmin_a,
#        vmax=vmax_a
#    )
#    im6.set_rasterized(True)

    levels_a = np.linspace(0, 5.0, 100)

    im5 = axes[2,0].contourf(
        X_pr, Y_pr, A_pr,
        levels=levels_a,
        extend='max'
    )
    im5.set_rasterized(True)
    im6 = axes[2,1].contourf(
        X_def, Y_def, A_def,
        levels=levels_a,
        extend='max'
    )
    im6.set_rasterized(True)

    # =========================================================
    # Titles
    # =========================================================

    axes[0,0].set_title(f"Si pristine", fontsize=18)
    axes[0,1].set_title(f"{defect}", fontsize=18)


    # =========================================================
    # Labels and formatting
    # =========================================================

    row_labels = [r"$r_s$", r"$s$", r"$\alpha$"]

    for i in range(3):

        axes[i,0].set_ylabel("y (Å)", fontsize=14)
        axes[i,1].set_ylabel("y (Å)", fontsize=14)

        axes[i,0].set_xlabel("x (Å)", fontsize=14)
        axes[i,1].set_xlabel("x (Å)", fontsize=14)

        axes[i,0].axis("equal")
        axes[i,1].axis("equal")

        # ingredient label on left side
        axes[i,0].text(
            -0.18, 0.5,
            row_labels[i],
            transform=axes[i,0].transAxes,
            fontsize=20,
            rotation=90,
            va="center"
        )

    # =========================================================
    # Colorbars
    # =========================================================

    cbar1 = fig.colorbar(
        im1,
        ax=axes[0,:],
        location="right",
        fraction=0.046,
        pad=0.04
    )

    #cbar1.set_label(r"$r_s$", fontsize=16)


    cbar2 = fig.colorbar(
        im3,
        ax=axes[1,:],
        location="right",
        fraction=0.046,
        pad=0.04
    )

    #cbar2.set_label(r"$s$", fontsize=16)


    cbar3 = fig.colorbar(
        im5,
        ax=axes[2,:],
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
        hspace=0.25,
        wspace=0.15
    )

    plt.savefig(
        f"{defect}_2d_ingredients.pdf",
        dpi=600,
        bbox_inches="tight"
    )

    #plt.show()
    plt.close()

#defects = ["Al","Ni","Pd","Pt","Cu","Ag","Au","Pb"]
#defects = ["T", "X", "H"]
#defects = ["T"]
defects = ["H"]
for defect in defects:
    print(f"Plotting {defect}...")
    plot_2d_ingredients(defect)
