# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 13:36:09 2023

Functions to calculate the various Qs of a resonator

References

(1) IEEE Trans. mic. theory tech. 43, 8, 1995
    Unloaded Q measurenebts
    En-Yuan Sun

@author: Daniel S
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as mt

def CritPoints(freqs):
    f1 = freqs[0];
    f2 = freqs[1];
    f3 = freqs[2];
    f4 = freqs[3];
    
    w1 = 2*np.pi*f1;
    w2 = 2*np.pi*f2;
    w3 = 2*np.pi*f3;
    w4 = 2*np.pi*f4;

    a = 1 + ((f1**2+f2**2+2*f1*f2-4*f3*f4)/(2*(f1*f3 + f1*f4 + f2*f3 + f2*f4 + 4*f3*f4)));

    b = ((f4 - f3 - ((f1+f2)/2)**2 * (1/f4 - 1/f3))/(2*(f2-f1)))**2;

    x = np.sqrt(((b - 2*a -1) + np.sqrt((b - 2*a -1)**2 - 4*(b+a)*(a-1)))/(2*(b+a)));

    Q0 = x*(w1 + w2)/(2*abs(w1-w2))

    print("a = "+str(a))
    print("b = "+str(b))
    print("x = "+str(x))
    print("Q0 = "+str(Q0))
    
    return(Q0)
    
def threedBNew(file):
    dat = pd.read_csv(file, names = ["freq","s11","phase"]);
    
    S11r = min(dat.s11)# Return loss at resonance
    f0 = dat.freq[np.where(np.array(dat.s11)==S11r)[0][0]]
    S11rLinear = 10**(S11r/20);
    S21rLinear = (1-S11rLinear**2)**0.5;
    S213dBLinear = S21rLinear/(2**0.5);
    S113dBLinear = (1-S213dBLinear**2)**0.5;
    S113dB = 20*mt.log10(S113dBLinear);
    
    differences = np.abs(abs(dat.s11) - S113dB)
    sorted_indices = np.argsort(differences)
    closest_indices = sorted_indices[:2]
    closest_values = dat.s11[closest_indices]

    f1 = dat.freq[closest_indices[0]]
    f2 = dat.freq[closest_indices[1]]
    deltaf = abs(f2-f1);
    
    
    Q = f0/np.abs(f2-f1);
    
    text = "$TEM_{013}$ \n $f_{0}$ = "+str(round(f0,3))+" GHz\n Δf ="+str(round(deltaf,3))+" GHz\n$Q_{0}$ = "+str(round(Q,3))
    bbox_props = dict(boxstyle="square,pad=0.3", facecolor="white", edgecolor="black")

    fig, ax = plt.subplots()
    ax.plot(dat.freq,dat.s11)
    ax.scatter(f0,S11r, color = "C1")
    ax.set_title("Q by 3 dB points")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("$S_{11}$ (dB)")
    ax.axvline(f1, color = "C1", linestyle = "dashed")
    ax.axvline(f2, color = "C1", linestyle = "dashed")
    ax.text(0.7, 0.05, text, transform=ax.transAxes, horizontalalignment='left',
             verticalalignment='bottom', bbox=bbox_props)
    
    
    
    
def threedB(file):
    dat = pd.read_csv(file, names = ["freq","s11","phase"]);
    
    S11r = min(dat.s11)# Return loss at resonance
    f0 = dat.freq[np.where(np.array(dat.s11)==S11r)[0][0]]
    
    differences = np.abs(abs(dat.s11) - 3)
    sorted_indices = np.argsort(differences)
    closest_indices = sorted_indices[:2]
    closest_values = dat.s11[closest_indices]

    f1 = dat.freq[closest_indices[0]]
    f2 = dat.freq[closest_indices[1]]
    deltaf = abs(f2-f1);
    
    
    Q = f0/np.abs(f2-f1);
    
    return(Q)
    
def threedB_plot(file):
    dat = pd.read_csv(file, names = ["freq","s11","phase"]);
    
    S11r = min(dat.s11)# Return loss at resonance
    f0 = dat.freq[np.where(np.array(dat.s11)==S11r)[0][0]]
    
    differences = np.abs(abs(dat.s11) - 3)
    sorted_indices = np.argsort(differences)
    closest_indices = sorted_indices[:2]
    closest_values = dat.s11[closest_indices]

    f1 = dat.freq[closest_indices[0]]
    f2 = dat.freq[closest_indices[1]]
    deltaf = abs(f2-f1);
    
    
    Q = f0/np.abs(f2-f1);
    
    text = "$TEM_{013}$ \n $f_{0}$ = "+str(round(f0,3))+" GHz\n Δf ="+str(round(deltaf,3))+" GHz\n$Q_{0}$ = "+str(round(Q,3))
    bbox_props = dict(boxstyle="square,pad=0.3", facecolor="white", edgecolor="black")

    fig, ax = plt.subplots()
    ax.plot(dat.freq,dat.s11)
    ax.scatter(f0,S11r, color = "C1")
    ax.set_title("Q by 3 dB points")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("$S_{11}$ (dB)")
    ax.axvline(f1, color = "C1", linestyle = "dashed")
    ax.axvline(f2, color = "C1", linestyle = "dashed")
    ax.text(0.7, 0.05, text, transform=ax.transAxes, horizontalalignment='left',
             verticalalignment='bottom', bbox=bbox_props)

def VSWR(file):
    dat = pd.read_csv(file, names = ["freq","s11","phase"]);
    
    S110 = (dat.s11[0]+dat.s11[len(dat.s11)-1])/2 # Reutrun loss away from resonance
    S11r = min(dat.s11)# Return loss at resonance
    f0 = dat.freq[np.where(np.array(dat.s11)==S11r)[0][0]]

    deltaS11r = abs(S11r) - abs(S110)

    SWRr = (1+10**(-deltaS11r/20))/(1-10**(-deltaS11r/20))

    deltaSWR = 1.3209 + 0.093099*SWRr**1 + 0.248572*SWRr**2 - 0.029787*SWRr**3 +0.00129*SWRr**4

    magDeltaS111 = (deltaSWR + SWRr - 1)/(deltaSWR + SWRr + 1)

    deltaS111dB = 20*np.log10(magDeltaS111)

    S111 = deltaS111dB - abs(S110)



    differences = np.abs(abs(dat.s11) - abs(S111))
    sorted_indices = np.argsort(differences)
    closest_indices = sorted_indices[:2]
    closest_values = dat.s11[closest_indices]

    f1 = dat.freq[closest_indices[0]]
    f2 = dat.freq[closest_indices[1]]
    deltaf = abs(f1-f2)
    Q0 = f0/deltaf
    
    return(Q0)
    

def VSWR_plot(file):
    dat = pd.read_csv(file, names = ["freq","s11","phase"]);
    
    S110 = (dat.s11[0]+dat.s11[len(dat.s11)-1])/2 # Reutrun loss away from resonance
    S11r = min(dat.s11)# Return loss at resonance
    f0 = dat.freq[np.where(np.array(dat.s11)==S11r)[0][0]]

    deltaS11r = abs(S11r) - abs(S110)

    SWRr = (1+10**(-deltaS11r/20))/(1-10**(-deltaS11r/20))

    deltaSWR = 1.3209 + 0.093099*SWRr**1 + 0.248572*SWRr**2 - 0.029787*SWRr**3 +0.00129*SWRr**4

    magDeltaS111 = (deltaSWR + SWRr - 1)/(deltaSWR + SWRr + 1)

    deltaS111dB = 20*np.log10(magDeltaS111)

    S111 = deltaS111dB - abs(S110)



    differences = np.abs(abs(dat.s11) - abs(S111))
    sorted_indices = np.argsort(differences)
    closest_indices = sorted_indices[:2]
    closest_values = dat.s11[closest_indices]

    f1 = dat.freq[closest_indices[0]]
    f2 = dat.freq[closest_indices[1]]
    deltaf = abs(f1-f2)
    Q0 = f0/deltaf

    text = "$TEM_{013}$ \n $f_{0}$ = "+str(round(f0,3))+" GHz\n Δf ="+str(round(deltaf,3))+" GHz\n$Q_{0}$ = "+str(round(Q0,3))
    bbox_props = dict(boxstyle="square,pad=0.3", facecolor="white", edgecolor="black")

    fig, ax = plt.subplots()
    ax.plot(dat.freq,dat.s11)
    ax.scatter(f0,S11r, color = "C1")
    ax.set_title("Q$_{0}$ by VSWR")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("$S_{11}$ (dB)")
    ax.axvline(f1, color = "C1", linestyle = "dashed")
    ax.axvline(f2, color = "C1", linestyle = "dashed")
    ax.text(0.7, 0.05, text, transform=ax.transAxes, horizontalalignment='left',
            verticalalignment='bottom', bbox=bbox_props)
   