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
    base_path = "/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_May19_metals_only"
    
    
    file_pr = f"{base_path}/{defect}/data_2D.dat"
    file_def = f"{base_path}/{defect}v/data_2D.dat"

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
    
    # Differences of the ingredients (Delta Ingredients)
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
    
    
   # Note Ratios, for Si is inverse: p/d, i.e., pristine over defective
   # Note for metals we do: d/p, i.e., defective over pristine.
   # ratio_pbe, ratio_scan, ratio_r2scan, ratio_lak, ratio_lda = [d/p for d, p in zip(exc_def, exc_pr)]
    ratio_pbe, ratio_scan, ratio_r2scan, ratio_lak, ratio_lda = [d/p for d, p in zip(exc_def, exc_pr)]
    
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
    
    all_data = np.concatenate([
    R_PBE.flatten(),
    R_SCAN.flatten(),
    R_R2SCAN.flatten(),
    R_LAK.flatten()
    ])

    ### Set limits for contours
    z_min = -0.01
    z_max = 0.01
    
    R_PBE_plot = R_PBE
    R_SCAN_plot = R_SCAN
    R_R2SCAN_plot = R_R2SCAN
    R_LAK_plot = R_LAK
    
    # Results using the function do not change much
#    R_PBE_plot = np.clip(R_PBE, z_min, z_max)
#    R_SCAN_plot = np.clip(R_SCAN, z_min, z_max)
#    R_R2SCAN_plot = np.clip(R_R2SCAN, z_min, z_max)
#    R_LAK_plot = np.clip(R_LAK, z_min, z_max)

    # =========================================================
    # Building Figure
    # =========================================================
    fig, axes = plt.subplots(4, 1, figsize=(10, 12))

    # =========================================================
    # r_s
    # =========================================================
    levels_rs = np.linspace(z_min, z_max, 100)
    im1 = axes[0].contourf(
    X_pr, Y_pr, R_PBE_plot,
    levels=levels_rs,
    cmap='viridis',
    extend='both'
    )
    im1.set_rasterized(True)
    
    im2 = axes[1].contourf(
        X_pr, Y_pr, R_LAK_plot,
        levels=levels_rs,
        cmap='viridis',
        extend='both'
    )
    im2.set_rasterized(True)

    im3 = axes[2].contourf(
    X_pr, Y_pr, R_SCAN_plot,
    levels=levels_rs,
    cmap='viridis',
    extend='both'
    )
    im3.set_rasterized(True)
    
    im4 = axes[3].contourf(
    X_pr, Y_pr, R_R2SCAN_plot,
    levels=levels_rs,
    cmap='viridis',
    extend='both'
    )
    im4.set_rasterized(True)
    
    # =========================================================
    # Titles
    # =========================================================
    axes[0].set_title(f"{defect}", fontsize=24)
    
    # =========================================================
    # Labels and formatting
    # =========================================================

    #row_labels = [r"$r_s$", r"$s$", r"$\alpha$"]
    row_labels = ["PBE", "LAK", "SCAN", "R2SCAN"]

    for i in range(4):
        axes[i].set_ylabel("y (Å)", fontsize=14)
        axes[i].set_xlabel("x (Å)", fontsize=14)
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

    cbar2 = fig.colorbar(
        im2,
        ax=axes[1],
        location="right",
        fraction=0.046,
        pad=0.04
    )
    
    cbar3 = fig.colorbar(
        im3,
        ax=axes[2],
        location="right",
        fraction=0.046,
        pad=0.04
    )
    
    cbar4 = fig.colorbar(
        im4,
        ax=axes[3],
        location="right",
        fraction=0.046,
        pad=0.04
    )

    # =========================================================
    # Layout + save
    # =========================================================

    plt.subplots_adjust(left=0.15)
    plt.subplots_adjust(
        left=0.15,
        right=0.88,
        hspace=0.4,
        wspace=0.15
    )

    plt.savefig(
        f"{defect}_delta_ratios_2d_ingredients.pdf",
        dpi=600,
        bbox_inches="tight"
    )

    #plt.show()
    plt.close()

defects = ["Al","Ni","Pd","Pt","Cu","Ag","Au","Pb"]

### Plotting ONLY one graph per metal including 4 DFAs BEGIN
#for defect in defects:
#    print(f"Plotting {defect}...")
#    plot_2d_ingredients(defect)
### Plotting ONLY one graph per metal including 4 DFAs END





#### Plotting in groups of 4 metals BEGIN ####
plt.rcParams.update({
    "axes.titlesize": 40,
    "xtick.labelsize": 28,
    "ytick.labelsize": 28,
})

def plot_page(defects):
    fig, axes = plt.subplots(4,4, figsize=(25,25))
    for j, defect in enumerate(defects):
        base_path = "/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_May19_metals_only"
        
        file_pr = f"{base_path}/{defect}/data_2D.dat"
        file_def = f"{base_path}/{defect}v/data_2D.dat"

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
        ratio_pbe, ratio_scan, ratio_r2scan, ratio_lak, ratio_lda = [d/p for d, p in zip(exc_def, exc_pr)]
        
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
        
        all_data = np.concatenate([
        R_PBE.flatten(),
        R_SCAN.flatten(),
        R_R2SCAN.flatten(),
        R_LAK.flatten()
        ])
        
        ## Limits for all plots ##
        z_min = -0.01
        z_max = 0.01
        
        R_PBE_plot = np.clip(R_PBE, z_min, z_max)
        R_SCAN_plot = np.clip(R_SCAN, z_min, z_max)
        R_R2SCAN_plot = np.clip(R_R2SCAN, z_min, z_max)
        R_LAK_plot = np.clip(R_LAK, z_min, z_max)
        
        levels_rs = np.linspace(z_min, z_max, 100)
        im1 = axes[0,j].contourf(
        X_pr, Y_pr, R_PBE_plot,
        levels=levels_rs,
        cmap='viridis',
        extend='both'
        )
        im1.set_rasterized(True)
        
        im2 = axes[1,j].contourf(
            X_pr, Y_pr, R_LAK_plot,
            levels=levels_rs,
            cmap='viridis',
            extend='both'
        )
        im2.set_rasterized(True)

        levels_a = np.linspace(z_min, z_max, 100)
        im3 = axes[2,j].contourf(
        X_pr, Y_pr, R_SCAN_plot,
        levels=levels_rs,
        cmap='viridis',
        extend='both'
        )
        im3.set_rasterized(True)
        
        im4 = axes[3,j].contourf(
        X_pr, Y_pr, R_R2SCAN_plot,
        levels=levels_rs,
        cmap='viridis',
        extend='both'
        )
        im4.set_rasterized(True)
        
        
    # =========================================================
    # Titles
    # =========================================================
        axes[0,j].set_title(f"{defect}", fontsize=40)


    # =========================================================
    # Labels and formatting
    # =========================================================
    row_labels = ["PBE","LAK","SCAN","R2SCAN"]
    
    for i in range(4):
        for j in range(4):
            axes[i,j].set_aspect('equal')

    for i in range(4):

        axes[i,0].text(
            -0.4,
            0.5,
            row_labels[i],
            transform=axes[i,0].transAxes,
            rotation=90,
            fontsize=40,
            va="center"
        )
            
    # =========================================================
    # Colorbars
    # =========================================================

    cbar1 = fig.colorbar(
            im1,
            ax=axes[0],
            location="right",
            fraction=0.024,
            pad=0.04
        )
        
    cbar2 = fig.colorbar(
        im2,
        ax=axes[1],
        location="right",
        fraction=0.024,
        pad=0.04
    )

    cbar3 = fig.colorbar(
        im3,
        ax=axes[2],
        location="right",
        fraction=0.024,
        pad=0.04
    )
    
    cbar4 = fig.colorbar(
        im4,
        ax=axes[3],
        location="right",
        fraction=0.024,
        pad=0.04
    )
    
    for cbar in [cbar1, cbar2, cbar3, cbar4]:
        cbar.ax.tick_params(labelsize=28)


    plt.subplots_adjust(
        left=0.15,
        right=0.88,
        hspace=0.4,
        wspace=0.15
    )

    plt.savefig(
        f"Delta_Rxc_ratios.pdf",
        dpi=600,
        bbox_inches="tight"
    )

    plt.close()
        
        

    
plot_page(["Al","Ni","Pd","Pt"])
#plot_page(["Cu","Ag","Au","Pb"])


#### Plotting in groups of 4 metals END ####
