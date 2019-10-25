#Forming observational dataset
import os
import numpy as np
import subprocess
UVB_Q = 18
dirname='/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final/'
run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe'
hden= -4
z = -1 #in log scale

stcolden= 14
if not os.path.exists(dirname):
     os.makedirs(dirname)
os.chdir(dirname)
filena = dirname + 'prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.in'
f=open(filena,"w+")
f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q= "+"{:.0f}".format(UVB_Q)+"] 
        \nhden " + "{:.0f}".format(hden)+ 
        "\nmetals " + "{:.2f}".format(z)+ " log 
        \nelement helium abundance 0.081632653 linear 
        \nstop column density "+"{:.1f}".format(stcolden)+" neutral H 
        \nconstant temperature, t=1e4 K [linear] 
        \nsave species column density \".spC\" no hash 
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
data=np.genfromtxt('/pyprog'+'{:.0f}'.format(UVB_Q) +'/Final/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.spC')
t=Table(data,names=['CII', 'CIII', 'CIV', 'NII', 'NIII', 
                    'NIV', 'NV', 'OI', 'OII', 'OIII', 'OIV', 
                    'OV', 'OVI', 'SIV', 'SV', 'SVI', 'SiII', 
                    'SiIII', 'SiIV'],
meta={'UVB_Q':UVB_Q,'metallicity': z,'n_H':hden,
'log(stop.colden)':stcolden})
t.write('/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final/colden_Oct12.fits')
read_obs=Table.read('/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final/colden_Oct12.fits')
read_obs.meta


#Ratios to use


colratios=['CII/CIII', 'CIII/CIV', 'NII/NIII', 'NIII/NIV','NIV/NV','OI/OII','OII/OIII','OIII/OIV','OIV/OV','OV/OVI','SIV/SV','SV/SVI','SiII/SiIII','SiIII/SiIV']
ratio_obs=[]
for i in range(len(data)):
    if i<2:
        ratio_obs.append(data[i]/data[i+1])
    if i>2 and i<6:
        ratio_obs.append(data[i]/data[i+1])
    if i>6 and i<12:
        ratio_obs.append(data[i]/data[i+1])
    if i>12 and i<15:
        ratio_obs.append(data[i]/data[i+1])
    if i>15 and i<18:
        ratio_obs.append(data[i]/data[i+1])

from astropy.table import Table
observed_ratios=Table([colratios,ratio_obs])
observed_ratios.write('/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final/obs_ratios.txt', format='ascii.tab',overwrite='True')

print("Done")
