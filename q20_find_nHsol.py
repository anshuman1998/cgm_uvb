#For Q20
import os
import numpy as np
import subprocess
UVB_Q = 20
dirname='/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final'
run='/cloudy.exe'
hden= np.arange(-5,-2.99,0.05)
z = 0.1
stcolden= 14

if not os.path.exists(dirname):
     os.makedirs(dirname)
os.chdir(dirname)

filena = dirname + '/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.in'
f=open(filena,"w+")

f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] 
        \nhden -4 vary \ngrid range from -5 to -3 with 0.05 dex steps 
        \nmetals " + "{:.4f}".format(z)+ 
        " log \nelement helium abundance 0.081632653 linear 
        \nstop column density "+"{:.1f}".format(stcolden)+
        " neutral H \nconstant temperature, t=1e4 K [linear] 
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

nHsol_allratios=[]

for i in range(len(ratio_arrays)):
    m,c=np.polyfit(hden,np.log10(ratio_arrays[i]),1)
    nHsol=(ratio_obs[i]-c)/m
    nHsol_allratios.append(nHsol)

#name_of_file= '/nH_sols_allratios.txt'
#np.savetxt(name_of_file,nHsol_allratios)

from astropy.table import Table
nHsol_ratios14=Table([colratios,nHsol_allratios])
nHsol_ratios14.write('/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final/nHsol_ratios14.txt', format='ascii.tab',overwrite='True')


print("Done")
