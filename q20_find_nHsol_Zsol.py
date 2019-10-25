#For Q20
import os
import numpy as np
import subprocess
UVB_Q = 20
dirname='/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final'
run='/cloudy.exe'
hden= np.arange(-5,-2.999,0.01)
z = -1 #Log Scale
stcolden= 14

if not os.path.exists(dirname):
     os.makedirs(dirname)
os.chdir(dirname)

filena = dirname + '/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.in'
f=open(filena,"w+")

f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] 
        \nhden -4 vary \ngrid range from -5 to -3 with 0.01 dex steps 
        \nmetals " + "{:.2f}".format(z)+ " log 
        \nelement helium abundance 0.081632653 linear 
        \nstop column density "+"{:.1f}".format(stcolden)+" neutral H 
        \nconstant temperature, t=1e4 K [linear] 
        \nsave species column density \".spC\"no hash 
        \n\"C+\" 
        \n\"C+2\" 
        \n\"C+3\" 
        \n\"N+\" 
        \n\"N+2\" 
        \n\"N+3\" 
        \n\"N+4\" 
        \n\"O\" 
        \n\"O+\"  
        \n\"O+2\" 
        \n\"O+3\" 
        \n\"O+4\" 
        \n\"O+5\" 
        \n\"S+3\" 
        \n\"S+4\" 
        \n\"S+5\" 
        \n\"Si+\" 
        \n\"Si+2\" 
        \n\"Si+3\" 
        \nend")
f.close()

process = subprocess.Popen([run, "prog"+"{:.0f}".format(UVB_Q)+"_Oct12.in"], 
stdout =subprocess.PIPE)
process.stdout.read()

from astropy.table import Table
data_20=np.genfromtxt('/pyprog'+'{:.0f}'.format(UVB_Q) +'/Final/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.spC')
t_20=Table(data,names=['CII', 'CIII', 'CIV', 'NII', 'NIII', 
                    'NIV', 'NV', 'OI', 'OII', 'OIII', 'OIV', 
                    'OV', 'OVI', 'SIV', 'SV', 'SVI', 'SiII', 
                    'SiIII', 'SiIV'],
meta={'UVB_Q':UVB_Q,'metallicity': z,
'log(stop.colden)':stcolden})
t_20.write('/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final/colden_Oct12.fits')
read=Table.read('/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final/colden_Oct12.fits')
read.meta

#Form list of column density ratios

colratios=['CII/CIII', 'CIII/CIV', 'NII/NIII', 'NIII/NIV',
           'NIV/NV','OI/OII','OII/OIII','OIII/OIV',
           'OIV/OV','OV/OVI','SIV/SV','SV/SVI',
           'SiII/SiIII','SiIII/SiIV']

ratio_arrays=[[] for i in range(len(hden))]

row=0

for j in ratio_arrays:
    for ion in range(len(data_20[row])):
        if ion<2:
            j.append(data_20[row][ion]/data_20[row][ion+1])
        if ion>2 and ion<6:
            j.append(data_20[row][ion]/data_20[row][ion+1])
        if ion>6 and ion<12:
            j.append(data_20[row][ion]/data_20[row][ion+1])
        if ion>12 and ion<15:
            j.append(data_20[row][ion]/data_20[row][ion+1])
        if ion>15 and ion<18:
            j.append(data_20[row][ion]/data_20[row][ion+1])
    row+=1

ratio_arrays=np.array(ratio_arrays)
ratio_arrays=ratio_arrays.transpose()
ratio_arrays=ratio_arrays.tolist()
from astropy.table import Table
q20_ratios14=Table([ratio_arrays[0],ratio_arrays[1],ratio_arrays[2],ratio_arrays[3],ratio_arrays[4],
                      ratio_arrays[5],ratio_arrays[6],ratio_arrays[7],ratio_arrays[8],ratio_arrays[9],
                      ratio_arrays[10],ratio_arrays[11],ratio_arrays[12],ratio_arrays[13]],
                    names=['CII/CIII', 'CIII/CIV', 'NII/NIII', 'NIII/NIV',
                           'NIV/NV','OI/OII','OII/OIII','OIII/OIV',
                           'OIV/OV','OV/OVI','SIV/SV','SV/SVI',
                           'SiII/SiIII','SiIII/SiIV'])
q20_ratios14.write('/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final/Q20_ratios14.txt', format='ascii.tab',overwrite='True')


#Fit the log(col_ratio) vs log(nH) and find solution log(nH)

import scipy.interpolate

nHsol_allratios=[]
spline_x = np.arange(-5,-2.999999,0.00001)
for i in range(len(ratio_arrays)):
    cubic_spline = scipy.interpolate.interp1d(hden,np.log10(ratio_arrays[i]),kind='cubic')
    spline_y=cubic_spline(spline_x)
    cn=0
    for j in range(len(spline_y)):
        if cn==0:
            if round(np.log10(ratio_obs[i]),4)==round(spline_y[j],4):
                cn+=1
                nHsolspline=spline_x[j]
                nHsol_allratios.append(nHsolspline)

from astropy.table import Table
nHsol_ratios14=Table([colratios,nHsol_allratios])
nHsol_ratios14.write('/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final/nHsol_ratios14.txt', format='ascii.tab',overwrite='True')

#Store for each element the mean and std dev of nH

carbon=[nHsol_allratios[0],nHsol_allratios[1]]
c_nH=round(np.mean(carbon),2)
c_err=round(np.std(carbon),2)
nitrogen=[nHsol_allratios[2],nHsol_allratios[3],nHsol_allratios[4]]
n_nH=round(np.mean(nitrogen),2)
n_err=round(np.std(nitrogen),2)
oxygen=[nHsol_allratios[5],nHsol_allratios[6],nHsol_allratios[7],nHsol_allratios[8],nHsol_allratios[9]]
o_nH=round(np.mean(oxygen),2)
o_err=round(np.std(oxygen),2)
sulphur=[nHsol_allratios[10],nHsol_allratios[11]]
s_nH=round(np.mean(sulphur),2)
s_err=round(np.std(sulphur),2)
silicon=[nHsol_allratios[12],nHsol_allratios[13]]
si_nH=round(np.mean(silicon),2)
si_err=round(np.std(silicon),2)

mean_nH_Q20=[c_nH,n_nH,o_nH,s_nH,si_nH]

stdev_nH_Q20=[c_err,n_err,o_err,s_err,si_err]

print(mean_nH_Q20)
print(stdev_nH_Q20)

#Find solution of Z

nH=[]
for i in hden:
    nH.append(round(i,2))
colden_nHsol=[]
c_nH=mean_nH_Q20[0]
n_nH=mean_nH_Q20[1]
o_nH=mean_nH_Q20[2]
s_nH=mean_nH_Q20[3]
si_nH=mean_nH_Q20[4]

for i in range(len(nH)):
    if nH[i]==c_nH:
        print("C")
        colden_nHsol.append(data_20[:,0][i])
        colden_nHsol.append(data_20[:,1][i])
        colden_nHsol.append(data_20[:,2][i])
for i in range(len(nH)):
    if nH[i]==n_nH:
        print("N")
        colden_nHsol.append(data_20[:,3][i])
        colden_nHsol.append(data_20[:,4][i])
        colden_nHsol.append(data_20[:,5][i])
        colden_nHsol.append(data_20[:,6][i])
for i in range(len(nH)):  
    if nH[i]==o_nH:
        print("O")
        colden_nHsol.append(data_20[:,7][i])
        colden_nHsol.append(data_20[:,8][i])
        colden_nHsol.append(data_20[:,9][i])
        colden_nHsol.append(data_20[:,10][i])
        colden_nHsol.append(data_20[:,11][i])
        colden_nHsol.append(data_20[:,12][i])
for i in range(len(nH)):        
    if nH[i]==s_nH:
        print("S")
        colden_nHsol.append(data_20[:,13][i])
        colden_nHsol.append(data_20[:,14][i])
        colden_nHsol.append(data_20[:,15][i])
for i in range(len(nH)):        
    if nH[i]==si_nH:
        print("Si")
        colden_nHsol.append(data_20[:,16][i])
        colden_nHsol.append(data_20[:,17][i])
        colden_nHsol.append(data_20[:,18][i])

Zsol_20=[]
for i,j in zip(colden_nHsol,data):
    z_scaled = (j*(10**z))/i
    Zsol_20.append(z_scaled)

all_ions=['CII', 'CIII', 'CIV', 'NII', 'NIII', 
           'NIV', 'NV', 'OI', 'OII', 'OIII', 'OIV', 
           'OV', 'OVI', 'SIV', 'SV', 'SVI', 'SiII', 
           'SiIII', 'SiIV']
        
from astropy.table import Table
Zsol_allions=Table([all_ions,Zsol_20])
Zsol_allions.write('/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final/Zsol_allions.txt', format='ascii.tab',overwrite='True')

#Group metallicity solutions according to element

carbon=[Zsol_20[0],Zsol_20[1],Zsol_20[2]]
c_Z=round(np.mean(carbon),8)
c_err=round(np.std(carbon),8)
nitrogen=[Zsol_20[3],Zsol_20[4],Zsol_20[5],Zsol_20[6]]
n_Z=round(np.mean(nitrogen),8)
n_err=round(np.std(nitrogen),8)
oxygen=[Zsol_20[7],Zsol_20[8],Zsol_20[9],Zsol_20[10],Zsol_20[11],Zsol_20[12]]
o_Z=round(np.mean(oxygen),8)
o_err=round(np.std(oxygen),8)
sulphur=[Zsol_20[13],Zsol_20[14],Zsol_20[15]]
s_Z=round(np.mean(sulphur),8)
s_err=round(np.std(sulphur),8)
silicon=[Zsol_20[16],Zsol_20[17],Zsol_20[18]]
si_Z=round(np.mean(silicon),8)
si_err=round(np.std(silicon),8)

mean_Z_Q20=[c_Z,n_Z,o_Z,s_Z,si_Z]

stdev_Z_Q20=[c_err,n_err,o_err,s_err,si_err]

print("Done")
