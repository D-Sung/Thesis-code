# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 15:33:01 2023
Analysing Q measurmenets of DNP samples in EPR resonators taken on 26/10/23.
@author: Daniel S
"""
import pandas as pd
import matplotlib.pyplot as plt
import Q_calculations as Q

path = "D:\\OneDrive - University of St Andrews\\Resonators\\EPR variant\\23-10-26\\" # Point to data

Q.VSWR_plot(path+"23102614_undercoupled_mag.csv") # Calculate Q0 by VSWR method using undercoupled data
Q.VSWR_plot(path+"23102617_undercoupled_mag.csv")
Q.VSWR_plot(path+"23102601_undercoupled_magnitude.csv")

dat10mM = pd.read_csv(path+"23102615_crit_mag.csv",names = ["freq","s11","phase"]); # Read in critically coupled data
dat25mM = pd.read_csv(path+"23102618_crit_mag.csv",names = ["freq","s11","phase"]);
dat50mM = pd.read_csv(path+"23102602_crit_magnitude.csv",names = ["freq","s11","phase"]);

f010mM = dat10mM.freq[dat10mM.s11.idxmin()]; # Calculate f0 by finding position of minimum S11
f025mM = dat25mM.freq[dat25mM.s11.idxmin()];
f050mM = dat50mM.freq[dat50mM.s11.idxmin()];

crit_freqs = [[9.3925,9.305,9.1675,9.7175],[9.3525,9.2575,9.1125,9.6325],[9.3075,9.205,9.5775,9.04]]; # Criical points read from Smith Chart on VNA for each data set
Q0crit = []; # Initialise array to store Q values

for i in range(0,len(crit_freqs)):
    Q0crit.append(Q.CritPoints(crit_freqs[i])); # Calculate Q values via crit points

# Creating and saving plots
fig1, ax = plt.subplots()
ax.set_title("EPR cavity resonance with 4-Amino-Tempo + H$_2$O")
ax.plot(dat10mM.freq,dat10mM.s11, label = "10 mM")
ax.plot(dat25mM.freq,dat25mM.s11, label = "25 mM")
ax.plot(dat50mM.freq,dat50mM.s11, label = "50 mM")
ax.text(9.4,-30,'Q$_{0}$ = '+str(round(Q0crit[0])))
ax.text(9.29,-35,'Q$_{0}$ = '+str(round(Q0crit[1])))
ax.text(9.18,-32,'Q$_{0}$ = '+str(round(Q0crit[2])))
plotrange = 0.65;
ax.set_xlim(f025mM - plotrange/2,f025mM + plotrange/2)
ax.set_xlabel("Frequency (GHz)")
ax.set_ylabel("S$_{11}$ (dB)")
ax.legend(title = "Concentration")
fig1.savefig(path+"Q0 vs Tempo concentration.svg",dpi = 1000)

fig2, ax = plt.subplots()
ax.set_title("EPR cavity resonance with 4-Amino-Tempo + H$_2$O")
ax.plot(dat10mM.freq/f010mM,dat10mM.s11, label = "10 mM")
ax.plot(dat25mM.freq/f025mM,dat25mM.s11, label = "25 mM")
ax.plot(dat50mM.freq/f050mM,dat50mM.s11, label = "50 mM")
ax.set_xlim(0.97,1.03)
ax.set_xlabel("Normalised frequency")
ax.set_ylabel("S$_{11}$ (dB)")
ax.legend(title = "Concentration")
fig2.savefig(path+"Q0 vs Tempo concentration normalised frequency.svg",dpi = 1000)
