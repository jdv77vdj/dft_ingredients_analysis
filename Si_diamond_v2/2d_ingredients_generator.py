
import os
import argparse
import shutil
import subprocess
import re
import math
import sys
import numpy as np
#sys.path.append('/Users/tlebeda/Documents/scripts')
#BASE_PATH = '/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_March6/Silicon'

BASE_PATH = '/Users/jdv/Desktop/Tulane/Tulane-Research/Projects/Defects_SCAN/fcc_metals_vac_formation_energy/updated_ingredients_April19_2d_ingredients/Silicon'

sys.path.append(BASE_PATH)
import generate_Fx as Fx

#testset = 'BG'   # set of 14 band gaps from TASK-prr paper

systems = dict()
  # koff is the offset required to plot a binding region path that includes the vacancy, since the vacyncy is at different positions for different systems.
  # LCs in Angstrom, but twice the  (important for calculating the gradient!)

# The next ones are the same pristine system
#systems['Sit']     = {'lattice_const': 10.860000, 'kxoff': .375, 'kyoff': .375, 'kzoff': .375, 'yp': 0.0, 'z': 4.072}    # needs NO /2 in ELFCAR due to spin
#
#systems['Six']     = {'lattice_const': 10.860000, 'kxoff': .375, 'kyoff': .375, 'kzoff': .375, 'yp': 0.0, 'z': 6.787}    # needs NO /2 in ELFCAR due to spin

systems['Sih']     = {'lattice_const': 10.860000, 'kxoff': .375, 'kyoff': .375, 'kzoff': .375, 'yp': 0.0, 'z': 3.2}    # needs NO /2 in ELFCAR due to spin


#systems['T']   = {'lattice_const': 10.860000, 'kxoff': .375, 'kyoff': .375, 'kzoff': .375, 'yp': 0.0, 'z': 4.072}    # needs NO /2 in ELFCAR due to spin
#systems['X']   = {'lattice_const': 10.860000, 'kxoff': .375, 'kyoff': .375, 'kzoff': .375, 'yp': 0.0, 'z': 6.787}    # needs NO /2 in ELFCAR due to spin
systems['H']   = {'lattice_const': 10.860000, 'kxoff': .375, 'kyoff': .375, 'kzoff': .375, 'yp': 0.0, 'z': 3.2}    # needs NO /2 in ELFCAR due to spin


# Have this system on to break the loop, otherwise goes in an infinite loop
systems['Ir']    = {'lattice_const': 7.662000, 'kxoff': .5, 'kyoff': .5, 'kzoff': .5, 'yp': 0.25}

xyz_dir = '_xyz'

Eh2eV = 27.211385
bohr2A = 0.529177249
FxSCAN=0.0
FxTASK=0.0
FxPBE =0.0
FcSCAN=0.0
FcPBE =0.0

def rewrite_ELFCAR(nuclei=False, path=False, twodim=False, xc_func='lda'):    #TODO: nuclei=False, when remove nuclei is enabled
    for system, properties in systems.items():
        dir = os.path.join(BASE_PATH,system)
        # read ELFCAR
        ELFCAR = os.path.join(dir, 'ELFCAR')
        print(f"Checking file exists: {ELFCAR}")
        print("Exists?", os.path.isfile(ELFCAR))
        assert (os.path.isfile(ELFCAR)), ELFCAR + " not existing!"
        input = open(ELFCAR, 'r')
        NGX = int(re.findall(' ([0-9]+)[ \t]+[0-9]+[ \t]+[0-9]+', input.read())[0])
        NGY = NGX
        NGZ = NGX
        print(f'NGX: {NGX}')
        input = open(ELFCAR, 'r')
        ELF = re.findall(r' [0-9].[0-9]{5}[ E][\-\+]?[0-9]?[0-9]?', input.read())#jv
        #ELF = re.findall(' [0-9].[0-9]{5}[ E][\-\+]?[0-9]?[0-9]?', input.read())
        #len_ELF = int(len(ELF)/2)    # Divide by two because both spins included. #com
        len_ELF = int(len(ELF))   # use this line for spin-UNpolarized systems
        if NGX*NGX*NGX!=len_ELF:
            print(f'ERROR:  NGX*NGY*NGZ != len(ELF) for system: {system} !!!')
        print(f'len(ELF): {len_ELF}')
        
        # read CHGCAR or AECCCAR2
        CHGCAR = os.path.join(dir, 'AECCAR2')#'CHGCAR')
        assert (os.path.isfile(CHGCAR)), CHGCAR + " not existing!"
        input = open(CHGCAR, 'r')
        RhoV = re.findall(r' [0-9\-].[0-9]{11}[ E][\-\+]?[0-9]?[0-9]?', input.read())
        #RhoV = re.findall(' [0-9\-].[0-9]{11}[ E][\-\+]?[0-9]?[0-9]?', input.read())
        len_RhoV = len(RhoV)
        if NGX*NGX*NGX!=len_RhoV:  # *8, if NGXF = 2*NGX and analogous for y, z
            print(f'ERROR:  NGXF*NGYF*NGZF != len(RhoV) for system: {system} !!!')
        print(f'len(RhoV): {len_RhoV}')
        
        # read AECCAR0
        COREDENS = os.path.join(dir, 'AECCAR0')#CHGCAR')
        assert (os.path.isfile(COREDENS)), COREDENS + " not existing!"
        input = open(COREDENS, 'r')
        RhoC = re.findall(r' [0-9\-].[0-9]{11}[ E][\-\+]?[0-9]?[0-9]?', input.read()) #jv
        
        #RhoC = re.findall(' [0-9\-].[0-9]{11}[ E][\-\+]?[0-9]?[0-9]?', input.read())
        len_RhoC = len(RhoC)
        if NGX*NGX*NGX!=len_RhoC:  # *8, if NGXF = 2*NGX and analogous for y, z
            print(f'ERROR:  NGXF*NGYF*NGZF != len(RhoC) for system: {system} !!!')
        print(f'len(RhoC): {len_RhoC}')
        
        # calculate x, y, z [k]
        lattice_const = properties['lattice_const']
        kxoff = properties['kxoff']
        kyoff = properties['kyoff']
        kzoff = properties['kzoff']
        yp = properties['yp']
        Vcell = (lattice_const/bohr2A)**3
        dx = lattice_const/float(NGX)
        x = [0.0]*len_ELF
        y = [0.0]*len_ELF
        z = [0.0]*len_ELF
        Rho = [0.0]*len_ELF  # Note: len_ELF=NGX*NGY*NGZ
        s = [0.0]*len_ELF
        alpha = [0.0]*len_ELF
        print(f'dx: {dx}')
        Rho_tot = 0.0
        Rho_val = 0.0
        for k in range(0,len_ELF):
            Rho[k] = float(RhoV[k])/Vcell + float(RhoC[k])/Vcell
     #       Rho_tot += (float(RhoV[k]) + float(RhoC[k]) ) / NGX**3.
            Rho_val += float(RhoV[k])/NGX**3.
     #   print(f'Gesamtdichte: {Rho_tot}')
        print(f'Valenzdichte: {Rho_val}')
        for k in range(0,len_ELF):
            x[k] = float(k%NGX) * dx
            y[k] = float((k//NGX)%NGY) * dx             # float(math.floor(k/NGX)%NGY) * dx             # k // N means floor(k/N)
            z[k] = float((k//(NGX*NGY))%NGZ) * dx       # float(math.floor((k/NGX)/NGY)%NGZ) * dx
     #       Rho[k] = float(RhoV[k])/Vcell
            if Rho[k] > 0.0:
                s[k] = calc_s_single_k(k, NGX, dx, Rho)
                alpha[k] = math.sqrt((1./float(ELF[k])-1.))
            else:
                s[k] = -1.0; alpha[k] = -1.0
        
        if path==False:
            if nuclei==False:
                if twodim==False:   # 3D
                    # write to file
                    plot_dir = os.path.join(BASE_PATH,system)
                    outfile = os.path.join(plot_dir, 'data.dat')
                    output = open(outfile, 'w')
                    with open(output.name, 'w', encoding=output.encoding) as o:
                        o.write(f'# System: {system}\tNGX: {NGX}\tNumber of points: {len_ELF}\n')
                       # o.write(f'# Nr   x   y   z in A   rho   rs   t   s   alpha \n\n')
                        o.write(f'# Nr \t x \t\t y \t \t z in A \t rho \t \t rs \t \t t \t \t s \t \t alpha \t \t FxPBE \t \t FxSCAN \t Fxr2SCAN \t FxLAK \t \t FcPBE \t \t FcSCAN \t Fcr2SCAN \t FcLAK \t FcLDA\n\n')
                        for k in range(0,len_ELF):
                            if Rho[k] > 0.0:
                                s = calc_s_single_k(k, NGX, dx, Rho)
                                alpha = math.sqrt((1./float(ELF[k])-1.))
                            else:
                                s = -1.0; alpha = -1.0
                         #   o.write(f'{k}\t{x[k]:.6f}\t{y[k]:.6f}\t{z[k]:.6f}\t{Rho[k]:.6f}\t{rs[k]:.6f}\t{t[k]:.6f}\t{s[k]:.6f}\t{alpha[k]:.6f}\n')
                            outstring = calc_write_single_k(k, x[k], y[k], z[k], Rho[k], s, alpha)
                            o.write(outstring)
                            
                if twodim==True:
                    # write to file, but select only one z-value and print emptylines, whenever y changes (for 2D-projection)
                    #plot_dir = os.path.join(BASE_PATH, xc_func, system)
                    plot_dir = os.path.join(BASE_PATH,system)
                    print("Current directory", plot_dir)
                    os.makedirs(plot_dir, exist_ok=True)
                    outfile = os.path.join(plot_dir, 'data_2D.dat')
                    
                    output = open(outfile, 'w')
                    
                    with open(output.name, 'w', encoding=output.encoding) as o:
                        o.write(f'# System: {system}\tNGX: {NGX}\tNumber of points: {len_ELF}\n')
                        o.write(f'# Nr \t x \t\t y \t \t z in A \t rho \t \t rs \t \t t \t \t s \t \t alpha \t \t FxPBE \t \t FxSCAN \t Fxr2SCAN \t FxLAK \t \t FcPBE \t \t FcSCAN \t Fcr2SCAN \t FcLAK \t FcLDA \n\n')
                        yprev=0.0

                        z_target = properties['z']
                        z_plane_index = np.argmin(np.abs(np.array(z) - z_target))
                        z_offset = (z_plane_index // (NGX*NGY)) * NGX * NGY
                        zfix = z[z_offset]

                        for k in range(0,NGX*NGY):
                            kloc = k + z_offset
                            if z[kloc]!=zfix:
                                print(f'ERROR in twodim: z has changed! !!!')
                            if y[kloc]!=yprev:
                                o.write(f'\n')
                                yprev=y[kloc]
                            if Rho[kloc] > 0.0:
                                s = calc_s_single_k(kloc, NGX, dx, Rho)
                                alpha = math.sqrt((1./float(ELF[kloc])-1.))
                            else:
                                s = -1.0; alpha = -1.0
                            outstring = calc_write_single_k(kloc, x[kloc], y[kloc], z[kloc], Rho[kloc], s, alpha)
                            o.write(outstring)
                    print("Done writing")

                            

        if path==True:                  # PATH: (0,0,0) atomic center -> (0,0.5,0.5) face center
            # write only path to file
                plot_dir = os.path.join(BASE_PATH, system)
                outfile = os.path.join(plot_dir, 'data_path_ae.dat') #data_path.dat')
                output = open(outfile, 'w')
                with open(output.name, 'w', encoding=output.encoding) as o:
                    o.write(f'# System: {system}\tNGX: {NGX}\tNumber of points: {len_ELF}\n')
                    #o.write(f'# NrPATH   x   y   z in A   rho   rs   t   s   alpha   NrCELL \n\n')
                    o.write(f'# Nr \t x \t\t y \t \t z in A \t rho \t \t rs \t \t t \t \t s \t \t alpha \t \t FxPBE \t \t FxSCAN \t Fxr2SCAN \t FxLAK \t \t FcPBE \t \t FcSCAN \t Fcr2SCAN \t FcLAK \t FcLDA \n\n')
                    i=0    # count points along PATH
                    # (0,0,0) -> (0,0.5,0.5)
                    for k in range(0,math.ceil(NGX/2)):
                        kloc = k*(1+int(4*yp)*NGX) + math.ceil(NGX*(kxoff-0.25) + NGX*(NGY*(kyoff-yp)) + NGX*NGY*(NGZ*kzoff) )  # start 0.25 before vacyncy; vacancy at (kxoff,kyoff,kzoff).   #  WARNING: 4*yp only works for yp=0 and yp=0.25 !!!
                        if Rho[kloc] > 0.0:
                            s = calc_s_single_k(kloc, NGX, dx, Rho)
                            alpha = math.sqrt((1./float(ELF[kloc])-1.))
                        else:
                            s = -1.0; alpha = -1.0
                        xloc = x[kloc] - (kxoff-0.25)*lattice_const # shift x (and y) to start at zero for easier plots (shifts vacancy to (0.25,0.25,z) )
                        yloc = y[kloc] - yp*lattice_const
                        #o.write(f'{i}\t{x[kloc]:.6f}\t{y[kloc]:.6f}\t{z[kloc]:.6f}\t{Rho[kloc]:.6f}\t{rs[kloc]:.6f}\t{t[kloc]:.6f}\t{s[kloc]:.6f}\t{alpha[kloc]:.6f}\t{kloc}\n')
                        outstring = calc_write_single_k(i, xloc, yloc, z[kloc], Rho[kloc], s, alpha)
                        o.write(outstring)
                        i = i+1
                    # manually add value for grid point out of cell
                    kloc = math.ceil(NGX*(kxoff+0.25) + NGX*(NGY*(kyoff+yp)) + NGX*NGY*(NGZ*kzoff) ) # LAST value of the loop
                    if Rho[0] > 0.0:
                        s = calc_s_single_k(kloc, NGX, dx, Rho)
                        alpha = math.sqrt((1./float(ELF[0])-1.))
                    else:
                        s = -1.0; alpha = -1.0
                    xloc = x[kloc] - (kxoff-0.25)*lattice_const # shift x (and y) to start at zero for easier plots (shifts vacancy to (0.25,0.25,z) )
                    yloc = y[kloc] - yp*lattice_const
                    outstring = calc_write_single_k(i, xloc, yloc, z[kloc], Rho[kloc], s, alpha)
                    o.write(outstring)
                    i = i+1
        print(f'{system} finished!\n')
        
        
def calc_write_single_k(k, x, y, z, Rho, s, alpha):
    if Rho > 0.0:
        rs = (4.0/3.0*math.pi*Rho)**(-1.0/3.0)
        t = (3.0*math.pi*math.pi/16.0)**(1.0/3.0) * s/rs
        FxLAK = Fx.calc_FxLAK()(s,alpha)
        Fxr2SCAN = Fx.calc_Fxr2SCAN()(s,alpha)
        FxSCAN = Fx.calc_FxSCAN()[0](s,alpha)
        FxTASK = Fx.calc_FxTASK()[0](s,alpha)
        FxPBE  = Fx.calc_FxPBE()[0](s)
        FxTASKd170 = Fx.calc_FxTASKd170()[0](s,alpha)
        FcLAK = Fx.calc_FcLAK()(rs,s,alpha)
        Fcr2SCAN = Fx.calc_Fcr2SCAN()(rs,s,alpha)
        FcSCAN = Fx.calc_FcSCAN()[0](rs,s,alpha)
        FcPBE  = Fx.calc_FcPBE()(rs,s)
        FcPW92 = Fx.calc_FcPW92()(rs)
        da_FxSCAN = Fx.calc_FxSCAN()[2](s,alpha)
        da_FxTASK = Fx.calc_FxTASK()[2](s,alpha)
        da_FxTASKd170 = Fx.calc_FxTASKd170()[2](s,alpha)
        da_FcSCAN = Fx.calc_FcSCAN()[1](rs,s,alpha)
    else:
        rs = -1.0
        t = -1.0
        FxLAK = -1.0
        Fxr2SCAN = -1.0
        FxSCAN = -1.0
        FxPBE  = -1.0
        FxTASK = -1.0
        FxTASKd170 = -1.0
        FcLAK = -1.0
        Fcr2SCAN = -1.0
        FcSCAN = -1.0
        FcPBE  = -1.0
        FcPW92 = -1.0
        da_FxSCAN = 10.0
        da_FxTASK = 10.0
        da_FxTASKd170 = 10.0
        da_FcSCAN = 10.0
    return f'{k}\t{x:.6f}\t{y:.6f}\t{z:.6f}\t{Rho:.6f}\t{rs:.6f}\t{t:.6f}\t{s:.6f}\t{alpha:.6f}\t{FxPBE:.6f}\t{FxSCAN:.6f}\t{Fxr2SCAN:.6f}\t{FxLAK:.6f}\t{FcPBE:.6f}\t{FcSCAN:.6f}\t{Fcr2SCAN:.6f}\t{FcLAK:.6f}\t{FcPW92:.6f}\n'
 #   return f'{k}\t{x:.6f}\t{y:.6f}\t{z:.6f}\t{Rho:.6f}\t{rs:.6f}\t{t:.6f}\t{s:.6f}\t{alpha:.6f}\t{FxSCAN:.6f}\t{FxTASK:.6f}\t{FxTASKd170:.6f}\t{FcSCAN:.6f}\t{FcPW92:.6f}\t{da_FxSCAN:.6f}\t{da_FxTASK:.6f}\t{da_FxTASKd170:.6f}\t{da_FcSCAN:.6f}\n'
    
def calc_s_single_k(k, NGX, dx, Rho):
    NGY = NGX; NGZ = NGX
    dxbohr = dx/bohr2A # For the gradient, dx is required in Bohr (because the density is in 1/Bohr^3)!!
    # grad_x RHO
    if k%NGX==0:                                 # x on left border of cell; take value of "left cell"
        RHO_xp = Rho[k+1]           # Rho(x+dx)
        RHO_xm = Rho[k+NGX-1]           # Rho(x-dx)
    elif k%NGX==NGX-1:                           # x on right border of cell; take value of "right cell"
        RHO_xp = Rho[k-(NGX-1)]         # Rho(x+dx)
        RHO_xm = Rho[k-1]           # Rho(x-dx)
    else:
        RHO_xp = Rho[k+1]           # Rho(x+dx)
        RHO_xm = Rho[k-1]           # Rho(x-dx)
    # grad_y RHO
    if (k//NGX)%NGY==0:                          # y on left border of cell; take value of "left cell"     # <-> k%(NGX*NGY) < NGX
        RHO_yp = Rho[k+NGX]         # Rho(y+dy)
        RHO_ym = Rho[k+NGX*(NGY-1)]     # Rho(y-dy)
    elif (k//NGX)%NGY==NGX-1:                    # y on right border of cell; take value of "right cell"
        RHO_yp = Rho[k-NGX*(NGY-1)]     # Rho(y+dy)
        RHO_ym = Rho[k-NGX]         # Rho(y-dy)
    else:
        RHO_yp = Rho[k+NGX]         # Rho(y+dy)
        RHO_ym = Rho[k-NGX]         # Rho(y-dy)
    # grad_z RHO
    if (k//(NGX*NGY))%NGZ==0:                    # z on left border of cell; take value of "left cell"     # <-> k < NGX*NGY
        RHO_zp = Rho[k+NGX*NGY]     # Rho(z+dz)
        RHO_zm = Rho[k+NGX*NGY*(NGZ-1)]     # Rho(z-dz)
    elif (k//(NGX*NGY))%NGZ==NGX-1:              # z on right border of cell; take value of "right cell"
        RHO_zp = Rho[k-NGX*NGY*(NGZ-1)]     # Rho(z+dz)
        RHO_zm = Rho[k-NGX*NGY]     # Rho(z-dz)
    else:
        RHO_zp = Rho[k+NGX*NGY]     # Rho(z+dz)
        RHO_zm = Rho[k-NGX*NGY]     # Rho(z-dz)
    grad_x_RHO = (RHO_xp - RHO_xm)/(2.0*dxbohr)
    grad_y_RHO = (RHO_yp - RHO_ym)/(2.0*dxbohr)    # dy=dx
    grad_z_RHO = (RHO_zp - RHO_zm)/(2.0*dxbohr)    # dz=dx
 #   print(f'k={k}\n Rho[k]={Rho[k]}\n RHO_xp={RHO_xp}\n RHO_xm={RHO_xm}\n RHO_yp={RHO_yp}\n RHO_ym={RHO_ym}\n RHO_zp={RHO_zp}\n RHO_zm={RHO_zm}\n ')
    # s
    return math.sqrt( grad_x_RHO*grad_x_RHO + grad_y_RHO*grad_y_RHO + grad_z_RHO*grad_z_RHO ) / (2.0*(3.0*math.pi*math.pi)**(1.0/3.0) * Rho[k]**(4./3.) )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Rewrite CHGCAR and ELFCAR for plot of mGGA-input variables")
    parser.add_argument('-a', '--analyze', action='store_true',
                        help='analyze and rewrite CHGCAR and ELFCAR')
    parser.add_argument('-n', '--nuclei', action='store_true',
                        help='if TRUE: leave away region near nuclei')
    parser.add_argument('-p', '--path', action='store_true',
                        help='if TRUE: print along path in unit cell')
    parser.add_argument('-t', '--twodim', action='store_true',
                        help='if TRUE: blank line in output whenever z changes for 2D-plots')
    parser.add_argument('-f', '--functionals', default='lda',
                        metavar='F', nargs='+', help='functionals', type=str)
    args = parser.parse_args()
    if(args.analyze):
        for functional in args.functionals:
            rewrite_ELFCAR(args.nuclei, args.path, args.twodim, xc_func=functional)
            
