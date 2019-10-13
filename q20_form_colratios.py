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
        \nhden -4 vary \ngrid range from -5 to -3 with 0.05 dex steps \nmetals " + "{:.4f}".format(z)+ 
        " log \nelement helium abundance 0.081632653 linear \nstop column density "+"{:.1f}".format(stcolden)+
        " neutral H \nconstant temperature, t=1e4 K [linear] \nsave species column density \".spC\" no hash 
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

#ratio_arrays=np.array(ratio_arrays)
#ratio_arrays=ratio_arrays.transpose()
from astropy.table import Table
q20_ratios=Table([colratios,
                      ratio_arrays[0],ratio_arrays[1],ratio_arrays[2],ratio_arrays[3],ratio_arrays[4],
                      ratio_arrays[5],ratio_arrays[6],ratio_arrays[7],ratio_arrays[8],ratio_arrays[9],
                      ratio_arrays[10],ratio_arrays[11],ratio_arrays[12],ratio_arrays[13],ratio_arrays[14],
                      ratio_arrays[15],ratio_arrays[16],ratio_arrays[17],ratio_arrays[18],ratio_arrays[19],
                      ratio_arrays[20],ratio_arrays[21],ratio_arrays[22],ratio_arrays[23],ratio_arrays[24],
                      ratio_arrays[25],ratio_arrays[26],ratio_arrays[27],ratio_arrays[28],ratio_arrays[29],
                      ratio_arrays[30],ratio_arrays[31],ratio_arrays[32],ratio_arrays[33],ratio_arrays[34],
                      ratio_arrays[35],ratio_arrays[36],ratio_arrays[37],ratio_arrays[38],ratio_arrays[39],
                      ratio_arrays[40] 
                      ])
q20_ratios.write('/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final/Q20_ratios.txt', format='ascii.tab',overwrite='True')

print("Done")
