#!/usr/bin/env python3
#import os
#import argparse
#import shutil
#import subprocess
#import re
import math
import numpy as np
import matplotlib.pyplot as plt
#import sys

# EXCHANGE

# FxTASK
# FxLAK
# Fxr2SCAN, FxSCAN, FxPBE
# FxTASKd170, FxTL23v1, FxTASKmod

# DEF Fx TASK

def calc_FxTASK():
    hx0 = 1.174      # LO-bound for two-electron systems
    Na = 2
    Nb = 4
    a = np.array([0.938719, -0.076371, -0.0150899])   # GE4x und hx(0)=1
    b = np.array([-0.628591, -2.10315, -0.5, 0.103153, 0.128591])  # GE4x , fx(0)=1, fx(1)=0 und fx(alpha-infty) = -3
    c = 4.9479  # H-atom norm
    D = 10.  # for numerical stability
    
    gx = lambda s: calc_gx(c,s)
    hx1 = lambda s: ChebEv(a,Na,s*s)
    fx = lambda alpha: ChebEv(b,Nb,alpha)

    Fx = lambda s, alpha: hx0*gx(s) + (1.-fx(alpha)) * (hx1(s)-hx0) * gx(s)**D
#    print(f'Fx(0,1) = {Fx(0,1)}')
    
    ds2_gx = lambda s: calc_dgx(c,s)
    ds2_hx1 = lambda s: ChebDerEv(a,Na,s*s)
    ds2_Fx = lambda s, alpha: hx0*ds2_gx(s) + (1.-fx(alpha)) * ( ds2_hx1(s) * gx(s)**D + (hx1(s)-hx0) * D*gx(s)**(D-1)*ds2_gx(s) )
#    print(f'ds2_Fx(0,1) = {ds2_Fx(0,1)}')
    
    dalpha_fx = lambda alpha: ChebDerEv(b,Nb,alpha)
#    print(f'dalpha_fx(1) = {dalpha_fx(1)}')
    dalpha_Fx = lambda s, alpha: dalpha_fx(alpha) * (hx0 - hx1(s)) * gx(s)**D
#    print(f'dalpha_Fx(0,1) = {dalpha_Fx(0,1)}')
    
    mu_alpha = -0.209897      # in principle free, but fixed by separation of s and alpha in ansatz for Fx (due to satisfaction of GE4x)
    mu_s2 = (10. + 60.*mu_alpha) / 81.
    C_s4 = - (1606. - 50.*mu_alpha) / 18225.
    C_s2alpha = - (511. - 50.*mu_alpha) / 13500.
    C_alpha2 = - (73. - 50.*mu_alpha) / 5000.
    
#    print(f'mu_s2 = {mu_s2},  mu_alpha = {mu_alpha},  C_s4 = {C_s4}  C_s2alpha = {C_s2alpha}  C_alpha2 = {C_alpha2}')
    
    return Fx, ds2_Fx, dalpha_Fx
    
def calc_gx(c,x):
    return np.where(x>0., 1.0 - np.exp(-c*x**(-1.0/2.0+(x==0))), 1.0)  # the "x==0" is to avoid an error due to hypothetic division by zero...
    
def calc_dgx(c,x):
    return np.where(x>0., - c/4. * x**(-5./2.+3.*(x==0)) * np.exp(-c*x**(-1.0/2.0+(x==0))), 0.0)  # the "x==0" is to avoid an error due to hypothetic division by zero...
    
def ChebEv(coef, n, r):
    d = 0.
    dd = 0.
    x = lambda r: (r-1.)/(r+1.)
    
    for j in range(n):   # range(n) == range(0,n); for loop in C++: for(j=0;j<=n-1;j++)
        k = n-j
        sv = d
        d = 2.*x(r)*d - dd + coef[k]
        dd = sv
    return x(r)*d - dd + coef[0]
    
def ChDer(coef, n):
    coefder=np.zeros(n+1)
    coefder[n]   = 0.
    coefder[n-1] = 2.*n*coef[n]
    
    if n>=2 :
        for j in range(n-1):
            k = n-2 -j
            coefder[k] = coefder[k+2] + 2.*(k+1)*coef[k+1]
            
    coefder[0] = coefder[0]/2.
    return coefder
    
def ChebDerEv(coef, n, r):
    coefder=np.zeros(n)
    coefder = ChDer(coef,n)
    return ChebEv(coefder,n,r) * 2. / ((1.+r)*(1.+r))
    
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# DEF Fx LAK
    
def calc_FxLAK():
    hx0 = 1.174     # LO-bound for two-electron systems
    ax = 1.1        # Construction principle (p2)
    bx = 4.9479     # H-atom norm
    anum = 5.       # Numerical stability
    mu_ax = -0.209897      # in principle free, but fixed by separation of s and alpha in ansatz for Fx (due to satisfaction of GE4x)
    nu_ax = (73. - 50.*mu_ax) / 5000. # fixed by GE4x and mu_ax
    c1 = mu_ax / (hx0-1.)           # prefactor of (alpha-1) for GE2x
    c2 = (mu_ax+nu_ax) / (hx0-1.)   # prefactor of (alpha-1)^2 for GE4x
    
    gx = lambda s: calc_gx(bx,s)
    gnum = lambda s: calc_gnum(anum,s*s)
    hx1 = lambda s: calc_hx1LAK(mu_ax,hx0,ax,bx,s)
    fx = lambda alpha: calc_fxLAK(c1,c2,alpha)

    Fx = lambda s, alpha: hx0*gx(s) + (1.-fx(alpha)) * (hx1(s)-hx0) * gnum(s)
#    print(f'Fx(0,1) = {Fx(0,1)}')
    
    ds2_gx = lambda s: calc_dgx(bx,s)
    ds2_gnum = lambda s: calc_dgnum(anum,s*s)
  #  ds2_hx1 = lambda s: calc_dhx1LAK(mu_ax,hx0,ax,bx,s)
  #  ds2_Fx = lambda s, alpha: hx0*ds2_gx(s) + (1.-fx(alpha)) * ( ds2_hx1(s) * gnum(s) + (hx1(s)-hx0) * *ds2_gnum(s) )
#    print(f'ds2_Fx(0,1) = {ds2_Fx(0,1)}')
    
 #   dalpha_fx = lambda alpha: calc_dfxLAK(c1,c2,alpha)
#    print(f'dalpha_fx(1) = {dalpha_fx(1)}')
 #   dalpha_Fx = lambda s, alpha: dalpha_fx(alpha) * (hx0 - hx1(s)) * gx(s)**D
#    print(f'dalpha_Fx(0,1) = {dalpha_Fx(0,1)}')
    
    return Fx#, ds2_Fx, dalpha_Fx
    
def calc_hx1LAK(muax,hx0,ax,bx,s):
    s2 = s*s
    mu_sx = (10. + 60.*muax) / 81.
    nu_sx = - (1606. - 50.*muax) / 18225.
    gx = lambda s: calc_gx(bx,s)
    hxGE4 = lambda s: 1. + mu_sx*s2 + nu_sx*s2**2. + hx0*(1.-gx(s))
    kx = lambda s: np.exp(-ax*ax/ (s2*(1.+s2)+(s2==0)) )
    return hxGE4(s) + kx(s) * (ax - hxGE4(s))
    
def calc_fxLAK(c1,c2,a):
    fx = 2./np.pi * math.atan(np.pi/2. * (c1*(a-1)/(a+(a==0)) + c2*(a-1)**2.))
    return fx
    
def calc_gnum(c,x):
    return np.where(x>0., 1.0 - np.exp(-c*c*x**(-1.0+(x==0))), 1.0)  # the "x==0" is to avoid an error due to hypothetic division by zero...
    
def calc_dgnum(c,x):
    return np.where(x>0., - c*c * x**(-2.+(x==0)) * np.exp(-c*c*x**(-1.0+(x==0))), 0.0)  # the "x==0" is to avoid an error due to hypothetic division by zero...
    
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# DEF Fx r2SCAN
    
def calc_Fxr2SCAN():
    hx0 = 1.174
    c = 4.9479; c1x = 0.667; c2x = 0.8; dx = 1.24; k1 = 0.065 # parameters determined by appropriate norms
    eta = 0.001 # numerical parameters
    
    gx = lambda s: calc_gx(c,s)   # from SCAN
    hx1 = lambda s: calc_hx1r2SCAN(k1,eta,s)
    fx = lambda s, alpha: calc_fr2SCAN(c1x,c2x,dx,eta,s,alpha)
    
  #  Fx = lambda s, alpha: hx0*gx(s) + (1.-fx(alpha)) * (hx1(s,alpha)-hx0) * gx(s)      # TASK-notation! Change formula to SCANx AFTER CHECK!
    Fx = lambda s, alpha: (hx1(s) + fx(s,alpha) * (hx0 - hx1(s)) ) * gx(s)
#    print(f'FxSCAN(0,1) = {Fx(0,1)}')
    
    return Fx
    
def calc_hx1r2SCAN(k1,eta,s):
    s2 = s*s
    c0 = 1.; c1 = -0.667; c2 = -0.44455555; c3 = -0.663086601049; c4 = 1.45129704490; c5 = -0.887998041597; c6 = 0.234528941479; c7 = -0.023185843322
    dp2 = 0.361
    mu = 10./81.
    Ceta = 20./27. + 5./3.*eta
    C2x = (c1 + 2.*c2 + 3.*c3 + 4.*c4 + 5.*c5 + 6.*c6 + 7.*c7) * (1.174-1.)
    x = (Ceta*C2x*np.exp(-s2/dp2**4.)+mu)*s2
    
    return 1. + k1 - k1 / (1. + x/k1)
    
def calc_fr2SCAN(c1x,c2x,dx,eta,s,a):
    at = a / (1.+eta*5./3.*s*s)
    omat = 1.-at
    c0 = 1.; c1 = -0.667; c2 = -0.44455555; c3 = -0.663086601049; c4 = 1.45129704490; c5 = -0.887998041597; c6 = 0.234528941479; c7 = -0.023185843322
    fxsum = c0 + c1*at + c2*at**2. + c3*at**3. + c4*at**4. + c5*at**5. + c6*at**6. + c7*at**7.
    return np.where(a<1., np.exp(-c1x*at/(omat+(omat==0))), np.where(a>2.5, -dx*np.exp(c2x/(omat+(omat==0))), fxsum) )
    
    
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# DEF Fx SCAN
    
def calc_FxSCAN():
    hx0 = 1.174
    c = 4.9479; c1x = 0.667; c2x = 0.8; dx = 1.24; k1 = 0.065 # parameters determined by appropriate norms
    
    gx = lambda s: calc_gx(c,s)
    hx1 = lambda s, alpha: calc_hx1SCAN(k1,s,alpha)
    fx = lambda alpha: calc_fSCAN(c1x,c2x,dx,alpha)
    
  #  Fx = lambda s, alpha: hx0*gx(s) + (1.-fx(alpha)) * (hx1(s,alpha)-hx0) * gx(s)      # TASK-notation! Change formula to SCANx AFTER CHECK!
    Fx = lambda s, alpha: (hx1(s,alpha) + fx(alpha) * (hx0 - hx1(s,alpha)) ) * gx(s)
#    print(f'FxSCAN(0,1) = {Fx(0,1)}')
    
    ds2_gx = lambda s: calc_dgx(c,s)
    ds2_hx1 = lambda s, alpha: calc_ds2_hx1SCAN(k1,s,alpha)
    ds2_Fx = lambda s, alpha: hx0*ds2_gx(s) + (1.-fx(alpha)) * ( ds2_hx1(s,alpha) * gx(s) + (hx1(s,alpha)-hx0) * ds2_gx(s) )
#    print(f'ds2_FxSCAN(0,1) = {ds2_Fx(0,1)}')
    
    dalpha_hx1 = lambda s, alpha: calc_dalpha_hx1SCAN(k1,s,alpha)
    dalpha_fx = lambda alpha: calc_dalpha_fSCAN(c1x,c2x,dx,alpha)
    dalpha_Fx = lambda s, alpha: dalpha_fx(alpha) * (hx0 - hx1(s,alpha)) * gx(s)
#    print(f'dalpha_FxSCAN(0,1) = {dalpha_Fx(0,1)}')
    
    return Fx, ds2_Fx, dalpha_Fx
    
def calc_hx1SCAN(k1,s,a):
    s2 = s*s
    oma = 1.-a
    mu_AK = 10./81. # = mu_s2 for mu_alpha=0
    b2 = math.sqrt(5913./405000.); b1 = (511./13500.)/(2.*b2); b3 = 0.5; b4 = mu_AK*mu_AK/k1 - 1606./18225. - b1*b1
    x = mu_AK*s2 * (1. + b4*s2/mu_AK * np.exp(-np.abs(b4)*s2/mu_AK) ) + (b1*s2 + b2*oma * np.exp(-b3*oma*oma))**2
    
    return 1. + k1 - k1 / (1. + x/k1)

def calc_ds2_hx1SCAN(k1,s,a):
    s2 = s*s
    oma = 1.-a
    mu_AK = 10./81. # = mu_s2 for mu_alpha=0
    b2 = math.sqrt(5913./405000.); b1 = (511./13500.)/(2.*b2); b3 = 0.5; b4 = mu_AK*mu_AK/k1 - 1606./18225. - b1*b1
    x = mu_AK*s2 * (1. + b4*s2/mu_AK * np.exp(-np.abs(b4)*s2/mu_AK) ) + (b1*s2 + b2*oma * np.exp(-b3*oma*oma))**2
    ds2_x = mu_AK + b4*s2 * (2. - np.abs(b4)/mu_AK*s2) * np.exp(-np.abs(b4)*s2/mu_AK) + 2.*b1 * (b1*s2 + b2*oma * np.exp(-b3*oma*oma))

    return ds2_x / (1. + x/k1)**2
    
def calc_dalpha_hx1SCAN(k1,s,a):
    s2 = s*s
    oma = 1.-a
    mu_AK = 10./81. # = mu_s2 for mu_alpha=0
    b2 = math.sqrt(5913./405000.); b1 = (511./13500.)/(2.*b2); b3 = 0.5; b4 = mu_AK*mu_AK/k1 - 1606./18225. - b1*b1
    x = mu_AK*s2 * (1. + b4*s2/mu_AK * np.exp(-np.abs(b4)*s2/mu_AK) ) + (b1*s2 + b2*oma * np.exp(-b3*oma*oma))**2
    dalpha_x = 2.*b2 (2.*b3 * oma*oma - 1.)*np.exp(-b3*oma*oma) (b1*s2 + b2*oma * np.exp(-b3*oma*oma))
    
    return dalpha_x / (1. + x/k1)**2
    
def calc_fSCAN(c1,c2,d,a):
    oma = 1.-a    
    return np.where(a<1., np.exp(-c1*a/(oma+(oma==0))), np.where(a>1., -d*np.exp(c2/(oma+(oma==0))), 0.) )
    
def calc_dalpha_fSCAN(c1,c2,d,a):
    oma = 1.-a    
    return np.where(a<1., -c1/(oma*oma+(oma==0))*np.exp(-c1*a/(oma+(oma==0))), np.where(a>1., -d*c2/(oma*oma+(oma==0))*np.exp(c2/(oma+(oma==0))), 0.) )
    
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# DEF Fx PBE
    
def calc_FxPBE():
    k1 = 0.804; mu_PBE = 0.21951 # parameters determined by appropriate norms
    
    Fx = lambda s: 1. + k1 - k1 / (1. + mu_PBE*s*s/k1)
#    print(f'FxPBE(0) = {Fx(0)}')
    
    ds2_Fx = lambda s: mu_PBE / (1. + mu_PBE*s*s/k1)**2
#    print(f'ds2_FxPBE(0) = {ds2_Fx(0)}')

    return Fx, ds2_Fx
    
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# DEF Fx TASKd170

def calc_FxTASKd170():
    hx0 = 1.174      # LO-bound for two-electron systems
    Na = 2
    Nb = 4
    a = np.array([0.938719, -0.076371, -0.0150899])   # GE4x und hx(0)=1
    b = np.array([-0.628591, -2.10315, -0.5, 0.103153, 0.128591])  # GE4x , fx(0)=1, fx(1)=0 und fx(alpha-infty) = -3
    c = 4.9479  # H-atom norm
    D = 170.  # for numerical stability
    
    gx = lambda s: calc_gx(c,s)
    hx1 = lambda s: ChebEv(a,Na,s*s)
    fx = lambda alpha: ChebEv(b,Nb,alpha)

    Fx = lambda s, alpha: hx0*gx(s) + (1.-fx(alpha)) * (hx1(s)-hx0) * gx(s)**D
    
    ds2_gx = lambda s: calc_dgx(c,s)
    ds2_hx1 = lambda s: ChebDerEv(a,Na,s*s)
    ds2_Fx = lambda s, alpha: hx0*ds2_gx(s) + (1.-fx(alpha)) * ( ds2_hx1(s) * gx(s)**D + (hx1(s)-hx0) * D*gx(s)**(D-1)*ds2_gx(s) )
    
    dalpha_fx = lambda alpha: ChebDerEv(b,Nb,alpha)
    dalpha_Fx = lambda s, alpha: dalpha_fx(alpha) * (hx0 - hx1(s)) * gx(s)**D
    
    return Fx, ds2_Fx, dalpha_Fx
    
# DEF Fx TL23v1

def calc_FxTL23v1():#TODO
    hx0 = 1.174      # LO-bound for two-electron systems
    c = 4.9479  # H-atom norm
    D = 5.  # for numerical stability
    
    gx = lambda s: calc_gx(c,s)
    hx1 = lambda s: calc_hx1TL23(s*s)
    fx = lambda alpha: calc_fxTL23(alpha)
    hxinf = lambda s: calc_hxinfTL23(s*s)
    fxinf = lambda alpha: calc_fxinfTL23(alpha)

    Fx = lambda s, alpha: hx0*gx(s) + (1.-fx(alpha)) * (hx1(s)-hx0) * gx(s)**D + fxinf(alpha) * (hxinf(s) - (hx1(s)+fx(alpha)*(hx0-hx1(s)) ) ) *gx(s)**(2.*D)
#    print(f'Fx(0,1) = {Fx(0,1)}')
    
    ds2_gx = lambda s: calc_dgx(c,s)
    ds2_hx1 = lambda s: calc_ds2_hx1TL23(s*s)
    ds2_hxinf = lambda s: calc_ds2_hxinfTL23(s*s)
    ds2_Fx = lambda s, alpha: hx0*ds2_gx(s) + (1.-fx(alpha)) * ( ds2_hx1(s) * gx(s)**D + (hx1(s)-hx0) * D*gx(s)**(D-1)*ds2_gx(s) )
#    print(f'ds2_Fx(0,1) = {ds2_Fx(0,1)}')
    
    dalpha_fx = lambda alpha: ChebDerEv(b,Nb,alpha)
#    print(f'dalpha_fx(1) = {dalpha_fx(1)}')
    dalpha_Fx = lambda s, alpha: dalpha_fx(alpha) * (hx0 - hx1(s)) * gx(s)**D
#    print(f'dalpha_Fx(0,1) = {dalpha_Fx(0,1)}')
    
    mu_alpha = -0.209897      # in principle free, but fixed by separation of s and alpha in ansatz for Fx (due to satisfaction of GE4x)
    mu_s2 = (10. + 60.*mu_alpha) / 81.
    C_s4 = - (1606. - 50.*mu_alpha) / 18225.
    C_s2alpha = - (511. - 50.*mu_alpha) / 13500.
    C_alpha2 = - (73. - 50.*mu_alpha) / 5000.
    
#    print(f'mu_s2 = {mu_s2},  mu_alpha = {mu_alpha},  C_s4 = {C_s4}  C_s2alpha = {C_s2alpha}  C_alpha2 = {C_alpha2}')
    
    return Fx, ds2_Fx, dalpha_Fx
    
def calc_GE4xTL23v1(s):#TODO
    s2 = s*s
    mu_AK = 10./81. # = mu_s2 for mu_alpha=0
    b2 = math.sqrt(5913./405000.); b1 = (511./13500.)/(2.*b2); b3 = 0.5; b4 = mu_AK*mu_AK/k1 - 1606./18225. - b1*b1
    x = mu_AK*s2 * (1. + b4*s2/mu_AK * np.exp(-np.abs(b4)*s2/mu_AK) ) + (b1*s2 + b2*oma * np.exp(-b3*oma*oma))**2
    
    return 1. + k1 - k1 / (1. + x/k1)

def calc_ds2_GE4xTL23v1(s):#TODO
    s2 = s*s
    oma = 1.-a
    mu_AK = 10./81. # = mu_s2 for mu_alpha=0
    b2 = math.sqrt(5913./405000.); b1 = (511./13500.)/(2.*b2); b3 = 0.5; b4 = mu_AK*mu_AK/k1 - 1606./18225. - b1*b1
    x = mu_AK*s2 * (1. + b4*s2/mu_AK * np.exp(-np.abs(b4)*s2/mu_AK) ) + (b1*s2 + b2*oma * np.exp(-b3*oma*oma))**2
    ds2_x = mu_AK + b4*s2 * (2. - np.abs(b4)/mu_AK*s2) * np.exp(-np.abs(b4)*s2/mu_AK) + 2.*b1 * (b1*s2 + b2*oma * np.exp(-b3*oma*oma))

    return ds2_x / (1. + x/k1)**2
    
def calc_hx1TL23(s):#TODO
    s2 = s*s
    mu_AK = 10./81. # = mu_s2 for mu_alpha=0
    b2 = math.sqrt(5913./405000.); b1 = (511./13500.)/(2.*b2); b3 = 0.5; b4 = mu_AK*mu_AK/k1 - 1606./18225. - b1*b1
    x = mu_AK*s2 * (1. + b4*s2/mu_AK * np.exp(-np.abs(b4)*s2/mu_AK) ) + (b1*s2 + b2*oma * np.exp(-b3*oma*oma))**2
    
    return 1. + k1 - k1 / (1. + x/k1)

def calc_ds2_hx1TL23(s):#TODO
    s2 = s*s
    oma = 1.-a
    mu_AK = 10./81. # = mu_s2 for mu_alpha=0
    b2 = math.sqrt(5913./405000.); b1 = (511./13500.)/(2.*b2); b3 = 0.5; b4 = mu_AK*mu_AK/k1 - 1606./18225. - b1*b1
    x = mu_AK*s2 * (1. + b4*s2/mu_AK * np.exp(-np.abs(b4)*s2/mu_AK) ) + (b1*s2 + b2*oma * np.exp(-b3*oma*oma))**2
    ds2_x = mu_AK + b4*s2 * (2. - np.abs(b4)/mu_AK*s2) * np.exp(-np.abs(b4)*s2/mu_AK) + 2.*b1 * (b1*s2 + b2*oma * np.exp(-b3*oma*oma))

    return ds2_x / (1. + x/k1)**2
    
def calc_fxinfTL23(c1,c2,a):
    return np.where(a>0., c1* np.exp(-c2*x**(-4.0+(a==0))), 1.0)  # the "a==0" is to avoid an error due to hypothetic division by zero...
    
def calc_dalpha_fxinfTL23(c1,c2,a):#TODO
    return np.where(a>0., c1* np.exp(-c2*x**(-4.0+(a==0))), 1.0)  # the "a==0" is to avoid an error due to hypothetic division by zero...
    
    
# DEF Fx modTASK

def calc_FxmodTASK():
    hx0 = 1.174      # LO-bound for two-electron systems
    Na = 2
   # Nb = 4
    a = np.array([0.938719, -0.076371, -0.0150899])   # GE4x und hx(0)=1
   # b = np.array([-0.628591, -2.10315, -0.5, 0.103153, 0.128591])  # GE4x , fx(0)=1, fx(1)=0 und fx(alpha-infty) = -3
    c = 4.9479  # H-atom norm
    c1 = 3.341; c2 = 1.2062701            # constants in fx(alpha) for fx(1)=0 and fx'(1)=fxTASK'(1)   TODO: satisfy also fx''(1)=fxTASK''(1) for GE4x!
    D = 10.  # for numerical stability
    
    gx = lambda s: calc_gx(c,s)
    hx1 = lambda s: ChebEv(a,Na,s*s)
    fx = lambda alpha: calc_fxmodTASK(c1,c2,alpha)

    Fx = lambda s, alpha: hx0*gx(s) + (1.-fx(alpha)) * (hx1(s)-hx0) * gx(s)**D
#    print(f'Fx(0,1) = {Fx(0,1)}')
    
    ds2_gx = lambda s: calc_dgx(c,s)
    ds2_hx1 = lambda s: ChebDerEv(a,Na,s*s)
    ds2_Fx = lambda s, alpha: hx0*ds2_gx(s) + (1.-fx(alpha)) * ( ds2_hx1(s) * gx(s)**D + (hx1(s)-hx0) * D*gx(s)**(D-1)*ds2_gx(s) )
#    print(f'ds2_Fx(0,1) = {ds2_Fx(0,1)}')
    
    dalpha_fx = lambda alpha: calc_dalpha_fxmodTASK(c1,c2,alpha)
#    print(f'dalpha_fx(1) = {dalpha_fx(1)}')
    dalpha_Fx = lambda s, alpha: dalpha_fx(alpha) * (hx0 - hx1(s)) * gx(s)**D
#    print(f'dalpha_Fx(0,1) = {dalpha_Fx(0,1)}')
    
    mu_alpha = -0.209897      # in principle free, but fixed by separation of s and alpha in ansatz for Fx
    mu_s2 = (10. + 60.*mu_alpha) / 81.
    C_s4 = - (1606. - 50.*mu_alpha) / 18225.
    C_s2alpha = - (511. - 50.*mu_alpha) / 13500.
    C_alpha2 = - (73. - 50.*mu_alpha) / 5000.
    
#    print(f'mu_s2 = {mu_s2},  mu_alpha = {mu_alpha},  C_s4 = {C_s4}  C_s2alpha = {C_s2alpha}  C_alpha2 = {C_alpha2}')
    
    return Fx, ds2_Fx, dalpha_Fx
    
def calc_fxmodTASK(c1,c2,a):
    return np.where(a>0., 1. - c1 * np.exp(-c2/(a+(a==0))), 1.)
    
def calc_dalpha_fxmodTASK(c1,c2,a):
    return np.where(a>0., - c1*c2 /(a+(a==0))**2 * np.exp(-c2/(a+(a==0))), 0.)

    
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# DEF ex LDA
    
def calc_exLDA():
    Ax = -0.738558766382022405884230032680836  # -3/4 (3/pi)**1/3
    ex = lambda rs: Ax * (3./4. / np.pi)**(1./3.) / rs   # Ax*rho**(1/3)       Ex = int n*ex

    return ex
    
    
    
    
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# CORRELATION

# FcLAK
# Fcr2SCAN
# FcSCAN
# FcPBE
# FcPW92
    
    
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# DEF Fc LAK    (for zeta=0)

def calc_FcLAK():
    ecLAK = calc_ecLAK()
    exLDA = calc_exLDA()
    Fc = lambda rs, s, alpha: ecLAK(rs,s,alpha)/exLDA(rs)
#    ds2_Fc = lambda rs, s, alpha: ds2_ecLAK(rs,s,alpha)/exLDA(rs)
#    dalpha_Fc = lambda rs, s, alpha: dalpha_ecLAK(rs,s,alpha)/exLDA(rs)
    
    return Fc#, ds2_Fc, dalpha_Fc
    
def calc_ecLAK():
    muax = -0.209897
    mua = -muax/2.
    b1c = 0.0468; b2c = 0.205601; b3c = 2.85; chi0 = 1.55344; ac = 10 # parameters determined by appropriate norms and construction principles
    anum = 5. # parameter for numerical stability
    
    ecLDA0 = lambda rs: calc_ecLDA0(b1c,b2c, rs)
    w0 = lambda rs: np.exp(-ecLDA0(rs)/b1c) - 1.
    g0 = lambda s: (1. + 4.*chi0 * s*s)**(-1./4.)
    H0 = lambda rs, s: b1c * np.log(1. + w0(rs) * (1. - g0(s)) )
    ec0 = lambda rs, s: ecLDA0(rs) + H0(rs,s)   # Gc(zeta=0) = 1
    
    ecLSDA = lambda rs: calc_ecLSDA(rs)[0]
    H1 = lambda rs, s: calc_H1LAK(mua,muax, ac,b3c, rs,s)
    ec1 = lambda rs, s: ecLSDA(rs) + H1(rs,s)
    
    fc = lambda rs, alpha: calc_fcLAK(mua,muax, b1c,b2c, rs,alpha)
    gnum = lambda s: calc_gnum(anum,s*s)
    ec = lambda rs, s, alpha: ec0(rs,s) + (1.-fc(rs,alpha)) * (ec1(rs,s) - ec0(rs,s)) * gnum(s)
#    print(f'ecLAK(1,0,1) = {ec(1,0,1)}')
    
    return ec#, ds2_ec, dalpha_ec
        
def calc_H1LAK(mua,muax, ac,b3c, rs,s):
    ct = (3.*np.pi**2./16.)**(2./3.)
    t = lambda rs, s: np.sqrt(ct * s*s / rs)  # t(rs,s)
    Ax = -3./4. * (3./np.pi)**(1./3.); Ac = (4.*np.pi/3.)**(-1./3.) * Ax
    Cma = lambda rs: calc_CmaLAK(rs)   # depends on mua
    mu_sc = lambda rs: 10./81. * (Cma(rs)*(1.+6.*mua) - (1.+6.*muax) )
    bt = lambda rs: Ac/ct * mu_sc(rs)
    
    gamma = (1.-np.log(2.))/np.pi**2.
    w1 = lambda rs: np.exp(-calc_ecLSDA(rs)[0]/gamma) - 1.
    A = lambda rs: bt(rs) / (gamma * w1(rs))
    g1p = lambda rs, t: (1. + 4.*A(rs)*t*t)**(-1./4.)
    g2p = lambda rs, t: (1. + (A(rs)*t*t)**2. )**(-1.)
    g3p = lambda rs, t: (1. + ac*A(rs)*t*t)**(-1.)
    g1m = lambda rs, t: (1. - 4.*A(rs)*t*t)**(-1./4.)
    g2m = lambda rs, t: (1. + (A(rs)*t*t)**2. )**(-1.)
    g3m = lambda rs, t: (1. - (w1(rs)+b3c)*A(rs)*t*t )**(-1.)
    
    H1p = lambda rs, t: gamma * np.log(1. + w1(rs) * (1. - g1p(rs,t)) * (1. - g2p(rs,t) + g3p(rs,t)) )
    H1m = lambda rs, t: gamma * np.log(1. + w1(rs) * (1. - g1m(rs,t)) * (1. - g2m(rs,t) - g3m(rs,t)) )
    
    return np.where(bt(rs)>0., H1p(rs,t(rs,s)), np.where(bt(rs)<0., H1m(rs,t(rs,s)), 0.) )
    
def calc_ecLDA0(b1c,b2c, rs):
    ecLDA0 = -b1c / (1. + b2c*rs)
    return ecLDA0
    
def calc_fcLAK(mua,muax, b1c,b2c, rs,a):
    at = (a-1)/(rs*a)
    Cma = lambda rs: calc_CmaLAK(rs)   # depends on mua
    mu_ac = lambda rs: Cma(rs)*mua - muax
    Ax = -3./4. * (3./np.pi)**(1./3.); Ac = (4.*np.pi/3.)**(-1./3.) * Ax
    beta_at = lambda rs: -Ac*mu_ac(rs)
    ec0s0 = lambda rs: calc_ecLDA0(b1c,b2c, rs)
    ec1s0 = lambda rs: calc_ecLSDA(rs)[0]
    fcGE2 = lambda rs: beta_at(rs) / (ec1s0(rs) - ec0s0(rs))
    fc = 2./np.pi * math.atan(np.pi/2. * (fcGE2(rs)*at) )
    return fc
    
def calc_CmaLAK(rs):
    Cs0 = - 16.*np.pi*(3.*np.pi**2.)**(1./3.) * 2.568/(3000.*10./81.)
    CmaLAK = Cs0 * (1.+0.1*rs**0.65) / ( (1.+0.065*rs**0.9)*(1.+0.03*rs**1.2) )   # depends on mua
    return CmaLAK
    
    
# -----------------------------------------------------------------------------------------------------------------------------------------------------
# DEF Fc SCAN    (for zeta=0)

def calc_FcSCAN():
    ecSCAN, ds2_ecSCAN, dalpha_ecSCAN = calc_ecSCAN()
    exLDA = calc_exLDA()
    FcSCAN = lambda rs, s, alpha: ecSCAN(rs,s,alpha)/exLDA(rs)
    ds2_FcSCAN = lambda rs, s, alpha: ds2_ecSCAN(rs,s,alpha)/exLDA(rs)
    dalpha_FcSCAN = lambda rs, s, alpha: dalpha_ecSCAN(rs,s,alpha)/exLDA(rs)
    
    return FcSCAN, ds2_FcSCAN, dalpha_FcSCAN
    
def calc_ecSCAN():
    b3c = 0.125541   # from lower bound on xc-energies of 2-electron systems, Fxc <= 1.67082
    c1c = 0.64; c2c = 1.5; dc = 0.7; b1c = 0.0285764; b2c = 0.0889; # parameters determined by appropriate norms
    
    ecLDA0 = lambda rs: -b1c / (1. + b2c*np.sqrt(rs) + b3c*rs)
    w0 = lambda rs: np.exp(-ecLDA0(rs)/b1c) - 1.
    chi_inf = 0.128026
    ginf = lambda s: 1. / (1. + 4.*chi_inf * s*s)**(1./4.)
    ds2_ginf = lambda s: -chi_inf / (1. + 4.*chi_inf * s*s)**(5./4.)
    H0 = lambda rs, s: b1c * np.log(1. + w0(rs) * (1. - ginf(s)) )
    ds2_H0 = lambda rs, s: -b1c / (1. + w0(rs) * (1. - ginf(s)) ) * w0(rs)*ds2_ginf(s)
    ec0 = lambda rs, s: ecLDA0(rs) + H0(rs,s)   # Gc(zeta=0) = 1
    ds2_ec0 = lambda rs, s: ds2_H0(rs,s)   # Gc(zeta=0) = 1
    
    ecLSDA = lambda rs: calc_ecLSDA(rs)[0]
    H1 = lambda rs, s: calc_H1SCAN(rs,s)   # Note that s is in principle also zeta-dependent!
    ds2_H1 = lambda rs: calc_ds2_H1SCAN(rs)   # Note that s is in principle also zeta-dependent!
    ec1 = lambda rs, s: ecLSDA(rs) + H1(rs,s)
    ds2_ec1 = lambda rs: ds2_H1(rs)
    
    fc = lambda alpha: calc_fSCAN(c1c,c2c,dc,alpha)
    ec = lambda rs, s, alpha: ec1(rs,s) + fc(alpha) * (ec0(rs,s) - ec1(rs,s))
#    print(f'ecSCAN(1,0,1) = {ec(1,0,1)}')
    
    dalpha_fc = lambda alpha: calc_dalpha_fSCAN(c1c,c2c,dc,alpha)
    dalpha_ec = lambda rs, s, alpha: dalpha_fc(alpha) * (ec0(rs,s) - ec1(rs,s))
    ds2_ec = lambda rs, s, alpha: ds2_ec1(rs) + fc(alpha) * (ds2_ec0(rs,s) - ds2_ec1(rs))
#    print(f'dalpha_ecSCAN(1,0,1) = {dalpha_ec(1,0,1)}')
    
    return ec, ds2_ec, dalpha_ec
        
def calc_H1SCAN(rs,s):
    t = (3.*np.pi*np.pi/16.)**(1./3.) * s / np.sqrt(rs)  # t(rs,s)
    gamma = 0.031091
    beta0 = 0.066725
    beta = lambda rs: beta0 * (1. + 0.1*rs) / (1. + 0.1778*rs)
    w1 = lambda rs: np.exp(-calc_ecLSDA(rs)[0]/gamma) - 1.
    A = lambda rs: beta(rs) / (gamma * w1(rs))
    g = lambda rs, t: 1. / (1. + 4.*A(rs)*t*t)**(1./4.)
    
    return gamma * np.log(1. + w1(rs) * (1. - g(rs,t)) )
        
def calc_ds2_H1SCAN(rs):
    ds2_t2 = (3.*np.pi*np.pi/16.)**(2./3.) / rs  # dt²(rs,s)/ds²
    gamma = 0.031091
    beta0 = 0.066725
    beta = lambda rs: beta0 * (1. + 0.1*rs) / (1. + 0.1778*rs)
    w1 = lambda rs: np.exp(-calc_ecLSDA(rs)[0]/gamma) - 1.
    A = lambda rs: beta(rs) / (gamma * w1(rs))
    ds2_g = lambda rs: 1. / (1. + 4.*A(rs)*ds2_t2)**(1./4.)
    
    return gamma * np.log(1. + w1(rs) * (1. - ds2_g(rs)) )
    
def calc_fcSCAN(c1,c2,d,a):
    oma = 1.-a    
    return np.where(a<1., np.exp(-c1*a/(oma+(oma==0))), np.where(a>1., -d*np.exp(c2/(oma+(oma==0))), 0.) )
    
def calc_dalpha_fcSCAN(c1,c2,d,a):
    oma = 1.-a    
    return np.where(a<1., -c1/(oma*oma+(oma==0))*np.exp(-c1*a/(oma+(oma==0))), np.where(a>1., -d*c2/(oma*oma+(oma==0))*np.exp(c2/(oma+(oma==0))), 0.) )
    
    
#  -----------------------------------------------------------------------------------------------------------------------------------------------------
# DEF Fc r2SCAN    (for zeta=0)

def calc_Fcr2SCAN():
    ecr2SCAN = calc_ecr2SCAN()
    exLDA = calc_exLDA()
    Fcr2SCAN = lambda rs, s, alpha: ecr2SCAN(rs,s,alpha)/exLDA(rs)
    
    return Fcr2SCAN
    
def calc_ecr2SCAN():
    b3c = 0.125541   # from lower bound on xc-energies of 2-electron systems, Fxc <= 1.67082
    c1c = 0.64; c2c = 1.5; dc = 0.7; b1c = 0.0285764; b2c = 0.0889; # parameters determined by appropriate norms
    eta = 0.001 # numerical parameters
    
    # alpha=0
    ecLDA0 = lambda rs: -b1c / (1. + b2c*np.sqrt(rs) + b3c*rs)
    decLDA0drs = lambda rs: -b1c / (1. + b2c*np.sqrt(rs) + b3c*rs)**2. * (b3c + b2c/(2.*np.sqrt(rs)))
    w0 = lambda rs: np.exp(-ecLDA0(rs)/b1c) - 1.
    chi_inf = 0.128026
    ginf = lambda s: 1. / (1. + 4.*chi_inf * s*s)**(1./4.)
    H0 = lambda rs, s: b1c * np.log(1. + w0(rs) * (1. - ginf(s)) )
    ec0 = lambda rs, s: ecLDA0(rs) + H0(rs,s)   # Gc(zeta=0) = 1
    
    # alpha=1
    ecLSDA = lambda rs: calc_ecLSDA(rs)[0]
    decLSDAdrs = lambda rs: calc_ecLSDA(rs)[1]
    
    t = lambda rs, s: (3.*np.pi*np.pi/16.)**(1./3.) * s / np.sqrt(rs)  # t(rs,s)
    dp2 = 0.361 # numerical parameter
    gamma = 0.031090690869655
    beta0 = 0.066725
    beta = lambda rs: beta0 * (1. + 0.1*rs) / (1. + 0.1778*rs)
    w1 = lambda rs: np.exp(-calc_ecLSDA(rs)[0]/gamma) - 1.
    y = lambda rs, t: beta(rs) / (gamma * w1(rs)) * t*t
    
    c0 = 1.; c1 = -0.64; c2 = -0.4352; c3 = -1.535685604549; c4 = 3.061560252175; c5 = -1.915710236206; c6 = 0.516884468372; c7 = -0.051848879792
    dfc2 = (c1 + 2.*c2 + 3.*c3 + 4.*c4 + 5.*c5 + 6.*c6 + 7.*c7)
    Dc = lambda rs: dfc2 / (27.*gamma*w1(rs)) * ( 20.*rs * (decLDA0drs(rs) - decLSDAdrs(rs)) - 45.*eta* (ecLDA0(rs) - ecLSDA(rs)) )
    dy = lambda rs, s: Dc(rs) * s*s * np.exp(- s**4./dp2**4.)
    
    g = lambda rs, s: 1. / (1. + 4.* (y(rs,t(rs,s)) - dy(rs,s)) )**(1./4.)
    
    H1 = lambda rs, s: gamma * np.log(1. + w1(rs) * (1. - g(rs,s)) )
    ec1 = lambda rs, s: ecLSDA(rs) + H1(rs,s)
    
    # Combine
    fc = lambda s, alpha: calc_fcr2SCAN(c1c,c2c,dc,eta,s,alpha)
    ec = lambda rs, s, alpha: ec1(rs,s) + fc(s,alpha) * (ec0(rs,s) - ec1(rs,s))
#    print(f'ecSCAN(1,0,1) = {ec(1,0,1)}')
    
    return ec
    
def calc_fcr2SCAN(c1c,c2c,dc,eta,s,a):
    at = a / (1.+eta*5./3.*s*s)
    omat = 1.-at
    c0 = 1.; c1 = -0.64; c2 = -0.4352; c3 = -1.535685604549; c4 = 3.061560252175; c5 = -1.915710236206; c6 = 0.516884468372; c7 = -0.051848879792
    fcsum = c0 + c1*at + c2*at**2. + c3*at**3. + c4*at**4. + c5*at**5. + c6*at**6. + c7*at**7.
    return np.where(a<1., np.exp(-c1c*at/(omat+(omat==0))), np.where(a>2.5, -dc*np.exp(c2c/(omat+(omat==0))), fcsum) )
    
    
#  -----------------------------------------------------------------------------------------------------------------------------------------------------
# DEF Fc PBE    (for zeta=0)

def calc_FcPBE():

    ecLSDA = lambda rs: calc_ecLSDA(rs)[0]
    H1 = lambda rs, s: calc_H1PBE(rs,s)   # Note that H1 is in principle also zeta-dependent!
    ecPBE = lambda rs, s: ecLSDA(rs) + H1(rs,s)
#    print(f'ecSCAN(1,0) = {ec(1,0)}')
    exLDA = calc_exLDA()
    FcPBE = lambda rs, s: ecPBE(rs,s)/exLDA(rs)
    
    return FcPBE
        
def calc_H1PBE(rs,s):
    t = (3.*np.pi*np.pi/16.)**(1./3.) * s / np.sqrt(rs)  # t(rs,s)
    gamma = 0.031091
    beta = 0.066725
    w1 = lambda rs: np.exp(-calc_ecLSDA(rs)[0]/gamma) - 1.
    A = lambda rs: beta / (gamma * w1(rs))
    g = lambda rs, t: 1. / (1. + A(rs)*t*t + (A(rs)*t*t)**2. )
    
    return gamma * np.log(1. + w1(rs) * (1. - g(rs,t)) ) # This is equivalent to beta/gamma * t^2 * g(At^2) with g = lambda rs, t: (1. + A(rs)*t*t) / (1. + A(rs)*t*t + (A(rs)*t*t)**2. )





# -----------------------------------------------------------------------------------------------------------------------------------------------------
# DEF ecLSDA    (for zeta=0)

def calc_FcPW92():
    ecLSDA = lambda rs: calc_ecLSDA(rs)[0]
    exLDA = calc_exLDA()
    FcPW92 = lambda rs: ecLSDA(rs)/exLDA(rs)

    return FcPW92

def calc_ecLSDA(rs):
    fzz = 1.7099209341613656176   # f''(0)
    zeta = 0.0    # NO spin-dependence implemented, so far!
    
    ec0=0.0; drs_ec0=0.0; ec1=0.0; drs_ec1=0.0; ac=0.0; drs_ac=0.0
    ec0, drs_ec0 = pw92_g(rs,0)
    ec1, drs_ec1 = pw92_g(rs,1)
    ac, drs_ac = pw92_g(rs,2)

    f = ((1.0+zeta)**(4./3.)+(1.0-zeta)**(4./3.)-2.0)/(2.0**(4./3.)-2.0)
    dzeta_f = 4./3.*((1.0+zeta)**(1./3.)-(1.0-zeta)**(1./3.))/(2.0**(4./3.)-2.0)
    zeta4 = zeta*zeta*zeta*zeta
    zeta3 = zeta*zeta*zeta

    ec       = ec0 +ac*f/fzz*(1.0-zeta4)+(ec1-ec0)*f*zeta4
    drs_ec   = drs_ec0+ drs_ac*f/fzz*(1.0-zeta4)+(drs_ec1-drs_ec0)*f*zeta4
    dzeta_ec = 4.*zeta3*f*(ec0-ec1-ac/fzz)+dzeta_f*(zeta4*(ec1-ec0)+(1.0-zeta4)*ac/fzz)
    
    return ec, drs_ec
    
def pw92_g(rs, version):
    if version==0: # g(rs)=ec(rs,0)
        p      = 1.0
        A      = 0.0310907
        alpha1 = 0.21370
        beta1  = 7.5957; beta2  = 3.5876; beta3  = 1.6382; beta4  = 0.49294
        g_sign = 1.0
    elif version==1: # g(rs)=ec(rs,1)
        p      = 1.0
        A      = 0.01554535
        alpha1 = 0.20548
        beta1  = 14.1189; beta2  = 6.1977; beta3  = 3.3662; beta4  = 0.62517
        g_sign = 1.0
    elif version==2: # g(rs)=-ac(rs)
        p      = 1.0
        A      = 0.0168869
        alpha1 = 0.11125
        beta1  = 10.357; beta2  = 3.6231; beta3  = 0.88026; beta4  = 0.49671
        g_sign = 1.0
    else: # undefined value => return zeros
        g  = 0.0
        dg = 0.0

    sqrt_rs = np.sqrt(rs)
    q0  = -2.0*A*(1.0+alpha1*rs)
    q1  = 2.0*A*(beta1*sqrt_rs+beta2*rs+beta3*sqrt_rs*rs+beta4*rs**(1.0+p))
    dq1 = A*(beta1/sqrt_rs+2.0*beta2+3.0*beta3*sqrt_rs+2.0*(1.0+p)*beta4*rs**p)

    g  = g_sign*q0*np.log(1.0+1./q1)
    dg = -g_sign*2.0*A*alpha1*np.log(1.0+1./q1)-q0*dq1/(q1*q1+q1)
    
    return g, dg
    
    

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# PLOT Fx

def plot_Fx():
    FxSCAN = calc_FxSCAN()[0]
    FxTASK = calc_FxTASK()[0]
    FxmodTASK = calc_FxmodTASK()[0]
 #   graph(FxSCAN, 0, 10, 0.001)
 #   graph2(FxSCAN, FxTASK, 0, 10, 0.001)
  #  graph2(FxmodTASK, FxTASK, 0, 10, 0.001)
  #  graph3(FxmodTASK, FxTASK, FxSCAN, 0, 10, 0.001)
    ds2_FxTASK = calc_FxTASK()[1]
    graph(ds2_FxTASK, 0, 10, 0.001)
    graph2(ds2_FxTASK, FxTASK, 0, 10, 0.001)

def graph(formula, xmin, xmax, dx):
    x = np.arange(xmin, xmax+dx, dx)
    #y = formula(x) # for some f(x)
    y = formula(x,1.5) # for Fx(s,alpha=1.5)
    plt.axis([xmin-dx, xmax+dx, min(y)-dx, max(y)+dx])  # set axis to just fit given area
    plt.xlabel('s', fontsize=20)
    plt.ylabel('Fx(s,alpha)', fontsize=20)
    plt.title('FxTASK', fontsize=20)
    plt.text(xmin + 0.8*(xmax-xmin+dx), min(y) + 0.8*(max(y)-min(y)+dx), r'$\alpha=1.5$', fontsize=20)
    plot=plt.plot(x,y)  
    plt.setp(plot, color='r', linewidth=2.0, linestyle='--', label='FxTASK')
    plt.show()
    
def graph2(formula1, formula2, xmin, xmax, dx):  # plot 2 functions in 1 figure
    x = np.arange(xmin, xmax+dx, dx)
    #y = formula(x) # for some f(x)
    y1 = formula1(x,100) # for Fx(s,alpha=1.5)
    y2 = formula2(x,100)
    plt.axis([xmin-dx, xmax+dx, min(min(y1),min(y2))-dx, max(max(y1),max(y2))+dx])  # set axis to just fit given area
    plt.xlabel('s', fontsize=20)
    plt.ylabel('Fx(s,alpha)', fontsize=20)
    plt.title('FxmodTASK vs FxTASK', fontsize=20)
    #plt.text(xmin + 0.8*(xmax-xmin+dx), min(min(y1),min(y2)) + 0.8*(max(max(y1),max(y2))-min(min(y1),min(y2))+dx), r'$\alpha=1.5$', fontsize=20)
    plot1=plt.plot(x,y1)
    plt.setp(plot1, color='r', linewidth=2.0, linestyle='--', label='FxmodTASK')
    plot2=plt.plot(x,y2)
    plt.setp(plot2, color='b', linewidth=2.0, linestyle='--', label='FxTASK')
    plt.show()
    
def graph3(formula1, formula2, formula3, xmin, xmax, dx):  # plot 2 functions in 1 figure
    x = np.arange(xmin, xmax+dx, dx)
    #y = formula(x) # for some f(x)
    y1 = formula1(0.5,x) # for Fx(s,alpha=1.5)
    y2 = formula2(0.5,x)
    y3 = formula3(0.5,x)
    plt.axis([xmin-dx, xmax+dx, min(min(y1),min(y2),min(y3))-dx, max(max(y1),max(y2),max(y3))+dx])  # set axis to just fit given area
    plt.xlabel(r'$\alpha$', fontsize=20)
    plt.ylabel(r'$F_x(s,\alpha)$', fontsize=20)
    plt.title('FxmodTASK vs FxTASK vs FxSCAN', fontsize=20)
    plt.text(xmin + 0.8*(xmax-xmin+dx), min(min(y1),min(y2)) + 0.8*(max(max(y1),max(y2))-min(min(y1),min(y2))+dx), r'$s=0.5$', fontsize=20)
    plot1=plt.plot(x,y1)
    plt.setp(plot1, color='r', linewidth=2.0, linestyle='--', label='FxmodTASK')
    plot2=plt.plot(x,y2)
    plt.setp(plot2, color='b', linewidth=2.0, linestyle='--', label='FxTASK')
    plot3=plt.plot(x,y3)
    plt.setp(plot3, color='g', linewidth=2.0, linestyle='--', label='FxSCAN')
    plt.show()



# -----------------------------------------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    
    #calc_ecLSDA(1.)
    #calc_FcSCAN()
    calc_FxTASK()
    #calc_FxmodTASK()
    #plot_Fx()
    
    
    
