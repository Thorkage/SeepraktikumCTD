#%%
"""
Einfaches Python Programm zur CTD-Auswertung des Seepraktikums des B.Sc. Ozeanographie im WiSe 2022/2023
Grundlage bilden SeaBird Daten. Diese wurden in unserem Fall bereits mit Loop- und Windowfilter vorprozessiert.

Zur Analyse wurden vier verschiedene Plots benutzt:
(1) Einfache Variable vs Tiefe Plots 
(2) TS Diagramme
(3) Heatmaps
(4) Schnitte über mehrere CTD Stationen

Autor: Torbjörn Kagel, torbjoern.kagel@outlook.de
"""

from cmath import nan
from pickletools import floatnl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.patheffects as pe
import pandas as pd
import os
import glob
import seawater 
import seaborn as sns
import cmocean

def clean_paths(file):
    str_path=str(file)
    str_path=str_path.replace('\\','/');str_path=str_path.replace('//','/');str_path=str_path.replace('[','');str_path=str_path.replace(']','');str_path=str_path.replace("'",'')
    return str_path

def get_start_line(file_path):
    f = open(file_path, "rb")
    data = f.read()
    k=0
    for i, line in enumerate(data.splitlines()):
        try:
            i = int(str(line)[3])
            break
        except:
            k+=1
            pass
    return k

path=os.getcwd()

#Die ctd_liste.csv enthält Metainfos über jeden Cast (bspw. Position, Uhrzeit, Tidenphase und dichtester Pegel)
ctd_list=pd.read_csv(os.path.join(path,"ctds_liste.csv"),encoding='latin1',sep=',',dayfirst=True,index_col=None)
ctd_list.set_index(ctd_list["CastNr"],inplace=True)
# ctd_list.set_index(ctd_list["CastNr"],inplace=True)

#List (cnvfiles) mit allen Pfaden der einzelnen Datendateien ersatellen
filenames = []
cnvfiles = []
cnvfiles_clean = []
d={}
fields=['julian days','pressure','temperature','conductivity','salinity','depth','oxygen saturation1','oxygen saturation2','turbidity','fluoresence','flag']

for i in ctd_list.index:
    if i < 10:
        k=0
    else:
        k=""

    #filename setzt sich aus Instrumenten Name + Datum + CastNr + letzter preprocessingschritt (wfilter) zusammen
    filenames.append("SBE19plus_01907321_????_??_??_Cast"+str(k)+str(i)+"_wfilter.cnv")
    cnvfiles.append(glob.glob((os.path.join(path,"Processed",filenames[len(filenames)-1])))) 
    str_path = clean_paths(cnvfiles[len(cnvfiles)-1])
    print(str_path)
    if str_path != "": 
        start_line = get_start_line(str_path)
        #es wird ein großes dict mit allen DataFrames der einzelnen Casts drin geschrieben
        d[i]=pd.read_csv(str_path,encoding='latin1',sep="\s+",header=None,index_col=False,names=fields,on_bad_lines='skip',skip_blank_lines=True,skiprows=start_line,comment='#')


#%% PLOT (1)

#User inputs:
desired_key = "Station"
key = "BE"
varx = "salinity"; unitx = "PSU"; limx = (25,30)
vary = "depth"; unity = "m"; limy = (0,15)
desired_entries = ctd_list.loc[ctd_list[desired_key] == key]


fig,ax = plt.subplots(dpi=100,figsize=(5,8))
n = len(desired_entries)
colors = plt.cm.jet(np.linspace(0,1,n))

j=0
for i in desired_entries.index:
    label='Cast: ' + str(ctd_list["CastNr"][i]) + ', '+ctd_list["Datum"][i] + ', ' + ctd_list["Uhrzeit [UTC]"][i][0:5] + '\n' +  ctd_list["Station"][i] + ', ' + ctd_list["Tidenzeitenphase (neu)"][i]
    ax.plot(d[i][varx],d[i][vary],
            label=label,
            color=colors[j])
    j+=1

ax.grid()
ax.set_xlabel(varx + " [" + unitx + "]")
ax.set_ylabel(vary + " [" + unity + "]")
ax.set_xlim(limx)
ax.set_ylim(limy)
ax.invert_yaxis()
fig.legend(loc='center left', bbox_to_anchor=(0.9, 0.5))
plt.show()

#%% PLOT (2)

### User inputs:
desired_key = "Tidenzeitenphase (neu)"
key = "HT"
varx = "temperature"; unitx = "Celsius"
vary = "salinity"; unity = "PSU"
varz = "depth"; unitz = "m"
desired_entries = ctd_list.loc[ctd_list[desired_key] == key]

### oder der einfache Weg:
# desired_entries = [1,2,3,4,5] #manuelle Eingabe von CastNrs, dann nur im Loop einmal die Liste einfügen.

n = len(desired_entries)
colors = plt.cm.jet(np.linspace(0,1,n))

minY=100
maxY=0
minX=100
maxX=0
minZ=100
maxZ=0

fig,ax = plt.subplots(dpi=200)
j=0
for i in desired_entries.index:
    label='Cast: ' + str(ctd_list["CastNr"][i]) + ', '+ctd_list["Datum"][i] + ', ' + ctd_list["Uhrzeit [UTC]"][i][0:5] + '\n' +  ctd_list["Station"][i] + ', ' + ctd_list["Tidenzeitenphase (neu)"][i]
    sct = ax.scatter(d[i][varx],d[i][vary],
                     c=d[i][varz],
                     cmap="Greys",
                     marker="2",
                     linewidths=0.1,
                     zorder=100)
    ax.scatter(d[i][varx],d[i][vary],
               c=colors[j],
               zorder=30,
               path_effects=[pe.Stroke(linewidth=10)],
               label=label)
    
    if np.nanmin(d[i][vary]) < minY: minY=np.nanmin(d[i][vary]) 
    if np.nanmax(d[i][vary]) > maxY: maxY=np.nanmax(d[i][vary]) 
    if np.nanmin(d[i][varx]) < minX: minX=np.nanmin(d[i][varx]) 
    if np.nanmax(d[i][varx]) > maxX: maxX=np.nanmax(d[i][varx])
    if np.nanmin(d[i][varz]) < minZ: minZ=np.nanmin(d[i][varz]) 
    if np.nanmax(d[i][varz]) > maxZ: maxZ=np.nanmax(d[i][varz])      
  
    j+=1

s = np.linspace(0.98*minX,1.02*maxX,1000)
t = np.linspace(0.98*minY,1.02*maxY,1000)
p = np.linspace(0.98*minZ,1.02*maxZ,1000)

S, T =np.meshgrid(s, t)
dens= seawater.dens(S,T,p) - 1000

#Konturplot der spezifischen Dichte
CS = plt.contour(S,T,dens, 
                 colors='gray', 
                 zorder = 5, 
                 linewidths=1)

ax.set_xlabel(varx + " [" + unitx + "]")
ax.set_ylabel(vary + " [" + unity + "]")
ax.set_xlim((0.98*minX,1.02*maxX))
ax.set_ylim((0.98*minY,1.02*maxY))
ax.clabel(CS)
# ax.set_aspect(1)
cbar = fig.colorbar(sct,
                    orientation="horizontal")

cbar.set_label(varz + " [" + unitz + "]")

fig.legend(loc='center left',
           bbox_to_anchor=(1, 0.5))
fig.suptitle(desired_key + ": " + key)
plt.tight_layout()
plt.show()

#%% PLOT (3)

'''
2D Histogram Plot, hauptsächlich um unterschiedliche Wassermassen (hier: ELBE, NORDSEE, MIXING) via CTD-Stationen zu unterscheiden.
Der ist noch bisschen ugly, damals hab ich da aus Zeitmangel viel manuell dran rumgeschraubt.
Eine Legend wäre super bspw.. Auch die einzelnen Histogramme normieren und dann ne colorbar zu plotten.
'''

varx = "salinity"; unitx = "PSU"; limx = (21,30)
vary = "temperature"; unity = "Celsius"; limy = (10,15)
watermass_key = ["NORDSEE","ELBE","MIXING"]

fig, ax = plt.subplots(dpi=150)

for key in watermass_key:
    desired_entries = pd.DataFrame()
    if key == "ELBE":
        stations = ["SW3", "NM", "NE2"]
        watermass_key_cmap = "Blues"

    elif key == "MIXING":
        stations = ["NFW"]
        watermass_key_cmap = "summer"

    elif key == "NORDSEE":
        stations = ["BE", "BM", "BW2", "NP1", "NP2", "NP3", "NP4", "NP5","TS","KFN"]
        watermass_key_cmap = "Reds"


    for j in stations:
        desired_entries = desired_entries.append(ctd_list.loc[ctd_list["Station"] == j])

    valx = []
    valy = []
    for i in desired_entries.index:
        try:
            valx.append(list(d[i][varx]))
            valy.append(list(d[i][vary]))
        except:
            pass

    valx = list(np.concatenate(valx))
    valy = list(np.concatenate(valy))

    hist = ax.hist2d(valx,valy,
            bins=30,
            cmap=watermass_key_cmap,
            cmin=1,
            range=[[limx[0], limx[1]], [limy[0], limy[1]]],
            zorder=100,
            label=key)

if varx == "salinity" and vary == "temperature":
    s = np.linspace(limx[0],limx[1],1000)
    t = np.linspace(limy[0],limy[1],1000)
    p = np.linspace(0,20,1000)

    S, T =np.meshgrid(s, t)
    dens= seawater.dens(S,T,p) - 1000

    #Konturplot der spezifischen Dichte
    CS = plt.contour(S,T,dens, 
                    colors='gray', 
                    zorder = 5, 
                    linewidths=1)
    ax.clabel(CS)

# fig.colorbar(hist[3], ax=ax, orientation="vertical")
# ax.set_aspect("equal")
ax.grid(zorder=0)
ax.set_xlabel(varx + " [" + unitx + "]")
ax.set_ylabel(vary + " [" + unity + "]")
ax.set_ylim(limy)
ax.set_xlim(limx)
#%% PLOT (4)
"""
    Hier werden manuell (ich glaube das zu automatisieren ist es nicht wert) CastNrs gewählt 
    und dann Schnitte geplottet
"""
varz= "temperature"; unitz = "PSU"; limz = (10,13); colormap=cmocean.cm.thermal
vary = "depth"; unity = "m"; limy = (0,6)
desired_entries_idx = [25,26,27,28,29]

desired_entries = ctd_list.loc[desired_entries_idx]
lats = desired_entries["Lat DecDeg"]
lons = desired_entries["LonDecDeg"]
desired_entries_df = pd.DataFrame()
zframe = pd.DataFrame()
yframe = pd.DataFrame()
d_des = {}
j=0
for i in desired_entries.index:
    zframe[i]=d[i][varz]
    yframe[i]=d[i][vary]


x = lons
y = np.mean(yframe,axis=1)
ymid = np.mean(y)
X,Y = np.meshgrid(x,y)

fig, ax = plt.subplots(dpi=150)
cont = ax.contourf(X,Y,zframe,cmap=colormap,vmin=limz[0],vmax=limz[1])
ax.set_ylim(limy)
ax.invert_yaxis()

ax.set_ylabel(vary + " [" + unity + "]")
ax.set_xlabel("Longitude")

j=0
for lon in lons:
    ax.axvline(lon,color="grey")
    plt.text(lon+0.005,ymid,desired_entries.at[desired_entries_idx[j],"Station"],rotation='vertical', ha='center', va='center')
    j+=1

cbar = fig.colorbar(cont,orientation="horizontal")
cbar.set_label(varz + " [" + unitz + "]")