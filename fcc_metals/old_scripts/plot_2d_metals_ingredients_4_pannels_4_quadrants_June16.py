#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams.update({
    "font.size": 30,          # base font
    "axes.titlesize": 34,     # subplot titles
    "axes.labelsize": 32,     # x/y labels
    "xtick.labelsize": 24,    # x tick labels
    "ytick.labelsize": 24,    # y tick labels
    "legend.fontsize": 24,
})

# =========================================================
# Load data
# =========================================================
base_path = "/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_May19_metals_only"

### Plotting function
def plot_2d_ingredients(metal):
    metals = ["Al","Ni","Pd","Pt","Cu","Ag","Au","Pb"]
    #=========================================================
    # Load data
    #=========================================================
    base_path = "/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_May19_metals_only"
    
    file_pr  = f"{base_path}/{metal}/data_2D.dat"
    file_def = f"{base_path}/{metal}v/data_2D.dat"

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
    
    # Differences of the ingredients (Delta ingredients)
    ds  = s_def  - s_pr
    da  = alpha_def  - alpha_pr
    drs = rs_def - rs_pr

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
    #=========================================================
    # Shared normalization
    #=========================================================

    vmin_rs = min(RS_pr.min(), RS_def.min())
    vmax_rs = max(RS_pr.max(), RS_def.max())

    vmin_a = min(A_pr.min(), A_def.min())
    vmax_a = max(A_pr.max(), A_def.max())

    # =========================================================
    # Figure
    # =========================================================
    fig, axes = plt.subplots(3, 3, figsize=(27, 24))


    ### Plain ingredients BEGIN ###
    ## r_s ##
    all_rs = np.concatenate([
    RS_pr.flatten(),
    RS_def.flatten(),
    DRS.flatten()
    ])
    vmax_rs = np.percentile(all_rs, 98)
    vmin_rs = np.percentile(all_rs, 2)
    levels_rs = np.linspace(vmin_rs, vmax_rs, 100)
    
    ##  s  ##
    all_s = np.concatenate([
        S_pr.flatten(),
        S_def.flatten()
    ])
    vmin_s = np.percentile(all_s, 2)
    vmax_s = np.percentile(all_s, 98)
    levels_s = np.linspace(vmin_s, vmax_s, 100)
    
    ## alpha  ##
    all_a = np.concatenate([
        A_pr.flatten(),
        A_def.flatten()
    ])
    vmin_a = np.percentile(all_a, 2)
    vmax_a = np.percentile(all_a, 98)
    levels_a = np.linspace(vmin_a, vmax_a, 100)
    ### Plain ingredients END ###
    
    
    ### Delta ingredients BEGIN ###
    ## delta r_s ##
    vmin_drs = -1
    vmax_drs = 1
    levels_drs = np.linspace(vmin_drs, vmax_drs, 100)

    ## delta d_s ##
    vmin_ds = -0.5
    vmax_ds = 0.5
    levels_ds = np.linspace(vmin_ds, vmax_ds, 100)

    ## delta d_a ##
    vmin_da = -0.1
    vmax_da = 0.1
    levels_da = np.linspace(vmin_da, vmax_da, 100)
    
    ### Delta ingredients END ###

    ### Plotting countors BEGIN ###
    
    #=========================================================
    # r_s
    #=========================================================

    im1 = axes[0,0].contourf(X_pr, Y_pr, RS_pr, levels=levels_rs, extend='both')
    im2 = axes[0,1].contourf(X_def, Y_def, RS_def, levels=levels_rs, extend='both')
    im3 = axes[0,2].contourf(X_pr, Y_pr, DRS, levels=levels_drs, cmap='RdBu_r', extend='both')
    
    #=========================================================
    # s
    #=========================================================
    im4 = axes[1,0].contourf(X_pr, Y_pr, S_pr, levels=levels_s, extend='both')
    im5 = axes[1,1].contourf(X_def, Y_def, S_def, levels=levels_s, extend='both')
    im6 = axes[1,2].contourf(X_pr, Y_pr, DS, levels=levels_ds, cmap='RdBu_r', extend='both')
    
    #=========================================================
    # alpha
    #=========================================================
    im7 = axes[2,0].contourf(X_pr, Y_pr, A_pr, levels=levels_a, extend='both')
    im8 = axes[2,1].contourf(X_def, Y_def, A_def, levels=levels_a, extend='both')
    im9 = axes[2,2].contourf(X_pr, Y_pr, DA, levels=levels_da, cmap='RdBu_r', extend='both')
    
    ### Plotting countors END ###



    #Rasterization
    im1.set_rasterized(True)
    im2.set_rasterized(True)
    im3.set_rasterized(True)
    im4.set_rasterized(True)
    im5.set_rasterized(True)
    im6.set_rasterized(True)
    im7.set_rasterized(True)
    im8.set_rasterized(True)
    im9.set_rasterized(True)
    
    
    # =========================================================
    # Titles
    # =========================================================

    FS = 40
    axes[0,0].set_title(f"{metal} pristine", fontsize=FS)
    axes[0,1].set_title(f"{metal} vacancy", fontsize=FS)
    axes[0,2].set_title(f"{metal} $\Delta$ Ingredient $def - pris$", fontsize=FS)
    
    # =========================================================
    # Labels and formatting
    # =========================================================

    row_labels = [r"$r_s$", r"$s$", r"$\alpha$"]

    for i in range(3):
        axes[i,0].axis("equal")
        axes[i,1].axis("equal")

        axes[i,0].text(
            -0.18, 0.5,
            row_labels[i],
            transform=axes[i,0].transAxes,
            fontsize=50,
            rotation=90,
            va="center"
        )

    # =========================================================
    # Set Colorbars
    # =========================================================

    # ---- Main ingredient colorbars, only appearing in middle column ----
    cbar1 = fig.colorbar(
        im1,
        ax=axes[0,0:2],
        fraction=0.025,
        pad=0.02
    )

    cbar2 = fig.colorbar(
        im4,
        ax=axes[1,0:2],
        fraction=0.025,
        pad=0.02
    )

    cbar3 = fig.colorbar(
        im7,
        ax=axes[2,0:2],
        fraction=0.025,
        pad=0.02
    )

    # ---- Delta ingredient colorbars ----
    cbar1_delta = fig.colorbar(
        im3,
        ax=axes[0,2],
        fraction=0.05,
        pad=0.04
    )

    cbar2_delta = fig.colorbar(
        im6,
        ax=axes[1,2],
        fraction=0.05,
        pad=0.04
    )

    cbar3_delta = fig.colorbar(
        im9,
        ax=axes[2,2],
        fraction=0.05,
        pad=0.04
    )

    # =========================================================
    # Font sizes
    # =========================================================
    for cbar in [cbar1, cbar2, cbar3,
                 cbar1_delta, cbar2_delta, cbar3_delta]:
        cbar.ax.tick_params(labelsize=28)
                
    
    #========================================================
    # Layout + save
    #=========================================================

    plt.subplots_adjust(left=0.15)
    
    plt.subplots_adjust(
    left=0.12,
    right=0.90,
    hspace=0.25,
    wspace=0.40
    )
    
    
    folder = "2d_ingredients_4all_quadrants"
    os.makedirs(folder, exist_ok=True)

    plt.savefig(
        f"{folder}/{metal}.pdf",
        dpi=600,
        bbox_inches="tight"
    )

    #plt.show()
    plt.close()


### Running code with function ###

#metals = ["Al","Ni","Pd","Pt","Cu","Ag","Au","Pb"]
#metals  = ["Ni","Pd","Pt","Cu","Ag","Au","Pb"]
metals = ["Pb"]
#metals = ["Al"]
for metal in metals:
    print(f"Plotting {metal}...")
    plot_2d_ingredients(metal)
