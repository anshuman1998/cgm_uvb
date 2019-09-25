#To call and run the Cloudy program

import os
import numpy as np
from scipy import interpolate
import subprocess

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from astropy.table import Table

#Iterative program for Q20
UVB_Q = 20
dirname='/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q) 
run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe' 

stcolden= 14 
#minstddev=1.0 
hden=np.arange(-4.6,-3.61,0.05)   #Decide the grid of nH

if not os.path.exists(dirname):      
    os.makedirs(dirname) 
os.chdir(dirname)
   
z=np.arange(0.01,1.02,0.05)
metal=[]



for i in z:
        metal.append(10 ** i)  
#Converting to linear metallicity values

name=['CII','CIII','CIV','CV','MgII','NeVIII',
'NII','NIII','NIV','OII','OIII','OIV','OV',
'OVI','SiII','SiIII','SiIV','SIV','SV','SVI'] 
#Ions for which calculation is performed.

print("Begin Iterative Program over n_H and Z")

for h in hden:
    print(h)
    filena = dirname + '/prog_metal'+'{:.0f}'.format(UVB_Q)+'_nH_'+'{:.2f}'.format(h)+'.in'     
    f=open(filena,"w+")     
    f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] 
            \nhden "+"{:.2f}".format(h)+" 
            \nmetals 0.11 log vary 
            \ngrid range 0.01 to 1.01 with 0.05 dex steps 
            \nelement helium abundance 0.081632653 linear 
            \nstop column density "+"{:.1f}".format(stcolden)+"neutral H 
            \nconstant temperature, t=1e4 K [linear] 
            \nsave species column density \".spC\" no hash 
            \n\"C+\" 
            \n\"C+2\" 
            \n\"C+3\" 
            \n\"C+4\" 
            \n\"Mg+\" 
            \n\"Ne+7\" 
            \n\"N+\" 
            \n\"N+2\" 
            \n\"N+3\" 
            \n\"O+\" 
            \n\"O+2\" 
            \n\"O+3\" 
            \n\"O+4\" 
            \n\"O+5\" 
            \n\"Si+\" 
            \n\"Si+2\" 
            \n\"Si+3\" 
            \n\"S+3\" 
            \n\"S+4\" 
            \n\"S+5\" 
            \nend")     
    f.close()     
    process = subprocess.Popen([run, "prog_metal"+"{:.0f}".format(UVB_Q)+'_nH_'
    +'{:.2f}'.format(h)+'.in'],stdout =subprocess.PIPE)     
    process.stdout.read()

print("Simulation run complete!")

ioncol=['red','blue','green','yellow','salmon',
'orchid','lime','dodgerblue','darkgoldenrod',
'fuchsia','sienna','dimgrey','navy','tomato','olive',
'black','cyan','darkorange','thistle','purple']
#Graphing colours used

for ion in range(len(name)):
    hsol=[]
    metsol=[]
    for h in hden:
        trial = np.genfromtxt("/Users/anshumanacharya/Downloads/
        c17.01/source/pyprog20/prog_metal20_nH_"+'{:.2f}'.format(h)+".spC")
        obs_vals=[0.452850,22.434700,19.903900,17.462100,0.001982,0.002153,
        0.076940,7.407480,9.796510,0.297943,23.139700,73.066600,20.292600,
        4.076360,0.010686,0.255324,1.030160,1.687530,0.517160,0.889996] 

        #Observed column densities of the 20 ions divided by 10^12 
        #for log(Z)=0.1,log(nH)=-4

        m,c= np.polyfit(metal,(trial[:,ion]/(10 ** 12)),1)
        #Slope and y intercept of linear fitting equation

        zsol=(obs_vals[ion]-c)/m
        #Solution value of metallicity for a metal

        if zsol>0:
            hsol.append(h)
            metsol.append(zsol)
       
       #Only save solutions with realistic metallicity values
       #That is, non-negative metallicity solutions

    print("Finished for:",name[ion])

    
    #n_H and Z solutions for all metals
    #20 files generated, one for each metal

    gridsol_Q20=Table([hsol,metsol],names=['n_H','Z'])
    gridsol_Q20.write("/Users/anshumanacharya/Downloads/
    c17.01/source/pyprog20/20Ions/gridsol_Q20_"+name[ion]+
    "_.txt", format='ascii.tab',overwrite='True')

    #Graphical representation of all solutions
    
    plt.scatter(metsol,hsol,label=name[ion],color=ioncol[ion])
    

print("Done!")

plt.ylabel('log(nH)')
plt.xlabel('Metallicity')

ax=plt.gca()
ax.set_xlim(0,3.0)
plt.legend(title='2D Grid For 20 ions',loc='upper center',
bbox_to_anchor=(0.5,1.05), ncol = 4,fancybox=True, shadow=True)
plt.show(block=False)
