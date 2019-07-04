#To call and run the Cloudy program
import os
import numpy as np
import subprocess
#for Q14
UVB_Q = 14 
dirname='/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q) 
run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe' 
stcolden= 14 
if not os.path.exists(dirname):      
    os.makedirs(dirname) 
os.chdir(dirname)     
filena = dirname + '/prog_metal'+'{:.0f}'.format(UVB_Q)+'.in'     
f=open(filena,"w+")     
f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] \nhden -4 \nmetals 0.11 log vary \ngrid range 0.01 to 1.01 with 0.05 dex steps \nelement helium abundance 0.081632653 linear \nstop column density "+"{:.1f}".format(stcolden)+" neutral H \nconstant temperature, t=1e4 K [linear] \nsave species column density \".spC\" no hash \n\"C+\" \n\"C+2\" \n\"C+3\" \n\"C+4\" \n\"Mg+\" \n\"Ne+7\" \n\"N+\" \n\"N+2\" \n\"N+3\" \n\"O+\" \n\"O+2\" \n\"O+3\" \n\"O+4\" \n\"O+5\" \n\"Si+\" \n\"Si+2\" \n\"Si+3\" \n\"S+3\" \n\"S+4\" \n\"S+5\" \nend")     
f.close()     
process = subprocess.Popen([run, "prog_metal"+"{:.0f}".format(UVB_Q)+".in"],stdout =subprocess.PIPE)     
process.stdout.read()

#To graph lin. colden/ 10^12 versus lin. Z
#Find slope and intercept for all ions
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
trial = np.genfromtxt("/Users/anshumanacharya/Downloads/c17.01/source/pyprog14/prog_metal14.spC")
name=[
'CII',
'CIII',
'CIV',
'CV',
'MgII',
'NeVIII',
'NII',
'NIII',
'NIV',
'OII',
'OIII',
'OIV',
'OV',
'OVI',
'SiII',
'SiIII',
'SiIV',
'SIV',
'SV',
'SVI']

z=np.arange(0.01,1.02,0.05)
metal=[]
slope=[]
intercept=[]
for i in z:
    metal.append(10 ** i)
for i in range(len(name)):
    plt.scatter(metal,(trial[:,i]/(10 ** 12)),label=name[i], marker='o')
    m,c= np.polyfit(metal,(trial[:,i]/(10 ** 12)),1)
    slope.append(m)
    intercept.append(c)
for i in range(len(slope)):
    y=slope[i]*np.asarray(metal)+ intercept[i]
    plt.plot(metal,y)
#plt.yscale('log')
plt.xlabel('Z/Zo')
plt.ylabel('Column Density/10^12')
plt.legend(title='Q14, n_H=-4, const. temperature 10000 K')
from astropy.table import Table
mc_colden_vs_Z=Table([name, slope, intercept],names=('Ion', 'Slope','Intercept'))
mc_colden_vs_Z.write("/Users/anshumanacharya/Desktop/Figs/mc_coldenvsZ_Q14_nH4.txt", format='ascii.tab',overwrite='True')
plt.show(block=False)