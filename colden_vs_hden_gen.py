#For getting column density ratios CIII/CIV and OIII/OIV versus hydrogen density
#First all models are run and the column density ratio data is stored.
#For Q20
import os
import numpy as np
import subprocess
UVB_Q = 20
dirname='/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final'
run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe'
hden= np.arange(-5,-2.999,0.01)
z = -1 #in log scale
stcolden= 14
if not os.path.exists(dirname):
     os.makedirs(dirname)
os.chdir(dirname)
filena = dirname + '/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.in'
f=open(filena,"w+")
f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] \nhden -4 vary \ngrid range from -5 to -3 with 0.01 dex steps \nmetals " + "{:.2f}".format(z)+ " log \nelement helium abundance 0.081632653 linear \nstop column density "+"{:.1f}".format(stcolden)+" neutral H \nconstant temperature, t=1e4 K [linear] \nsave species column density \".spC\" no hash \n\"C+\" \n\"C+2\" \n\"C+3\" \n\"N+\" \n\"N+2\" \n\"N+3\" \n\"N+4\" \n\"O\" \n\"O+\"  \n\"O+2\" \n\"O+3\" \n\"O+4\" \n\"O+5\" \n\"S+3\" \n\"S+4\" \n\"S+5\" \n\"Si+\" \n\"Si+2\" \n\"Si+3\" \nend")
f.close()
process = subprocess.Popen([run, "prog"+"{:.0f}".format(UVB_Q)+"_Oct12.in"], 
stdout =subprocess.PIPE)
process.stdout.read()

from astropy.table import Table
data_20=np.genfromtxt('/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q) +'/Final/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.spC')

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

#Store the ratios of the ions obtained in the simulation run

ratio_arrays=np.array(ratio_arrays)
ratio_arrays=ratio_arrays.transpose()
ratio_arrays=ratio_arrays.tolist()

c3_c4_q20=ratio_arrays[1]
o3_o4_q20=ratio_arrays[7]

print("Done for 20!!")

#For Q18 check
import os
import numpy as np
import subprocess
UVB_Q = 18
dirname='/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final'
run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe'
hden= np.arange(-5,-2.999,0.01)
z = -1
stcolden= 14
if not os.path.exists(dirname):
     os.makedirs(dirname)
os.chdir(dirname)
filena = dirname + '/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12check.in'
f=open(filena,"w+")
f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] \nhden -4 vary \ngrid range from -5 to -3 with 0.01 dex steps \nmetals " + "{:.2f}".format(z)+ " log \nelement helium abundance 0.081632653 linear \nstop column density "+"{:.1f}".format(stcolden)+" neutral H \nconstant temperature, t=1e4 K [linear] \nsave species column density \".spC\" no hash \n\"C+\" \n\"C+2\" \n\"C+3\" \n\"N+\" \n\"N+2\" \n\"N+3\" \n\"N+4\" \n\"O\" \n\"O+\"  \n\"O+2\" \n\"O+3\" \n\"O+4\" \n\"O+5\" \n\"S+3\" \n\"S+4\" \n\"S+5\" \n\"Si+\" \n\"Si+2\" \n\"Si+3\" \nend")
f.close()
process = subprocess.Popen([run, "prog"+"{:.0f}".format(UVB_Q)+"_Oct12check.in"], 
stdout =subprocess.PIPE)
process.stdout.read()

from astropy.table import Table
data_18=np.genfromtxt('/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q) +'/Final/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12check.spC')

colratios=['CII/CIII', 'CIII/CIV', 'NII/NIII', 'NIII/NIV',
           'NIV/NV','OI/OII','OII/OIII','OIII/OIV',
           'OIV/OV','OV/OVI','SIV/SV','SV/SVI',
           'SiII/SiIII','SiIII/SiIV']

ratio_arrays=[[] for i in range(len(hden))]

row=0

for j in ratio_arrays:
    for ion in range(len(data_18[row])):
        if ion<2:
            j.append(data_18[row][ion]/data_18[row][ion+1])
        if ion>2 and ion<6:
            j.append(data_18[row][ion]/data_18[row][ion+1])
        if ion>6 and ion<12:
            j.append(data_18[row][ion]/data_18[row][ion+1])
        if ion>12 and ion<15:
            j.append(data_18[row][ion]/data_18[row][ion+1])
        if ion>15 and ion<18:
            j.append(data_18[row][ion]/data_18[row][ion+1])
    row+=1

ratio_arrays=np.array(ratio_arrays)
ratio_arrays=ratio_arrays.transpose()
ratio_arrays=ratio_arrays.tolist()

c3_c4_q18=ratio_arrays[1]
o3_o4_q18=ratio_arrays[7]

print("Done for 18!!")

#For Q16
import os
import numpy as np
import subprocess
UVB_Q = 16
dirname='/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final'
run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe'
hden= np.arange(-5,-2.999,0.01)
z = -1
stcolden= 14
if not os.path.exists(dirname):
     os.makedirs(dirname)
os.chdir(dirname)
filena = dirname + '/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.in'
f=open(filena,"w+")
f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] \nhden -4 vary \ngrid range from -5 to -3 with 0.01 dex steps \nmetals " + "{:.2f}".format(z)+ " log \nelement helium abundance 0.081632653 linear \nstop column density "+"{:.1f}".format(stcolden)+" neutral H \nconstant temperature, t=1e4 K [linear] \nsave species column density \".spC\" no hash \n\"C+\" \n\"C+2\" \n\"C+3\" \n\"N+\" \n\"N+2\" \n\"N+3\" \n\"N+4\" \n\"O\" \n\"O+\"  \n\"O+2\" \n\"O+3\" \n\"O+4\" \n\"O+5\" \n\"S+3\" \n\"S+4\" \n\"S+5\" \n\"Si+\" \n\"Si+2\" \n\"Si+3\" \nend")
f.close()
process = subprocess.Popen([run, "prog"+"{:.0f}".format(UVB_Q)+"_Oct12.in"], 
stdout =subprocess.PIPE)
process.stdout.read()

from astropy.table import Table
data_16=np.genfromtxt('/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q) +'/Final/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.spC')


#Calculate the ratios of ions from simulation run

colratios=['CII/CIII', 'CIII/CIV', 'NII/NIII', 'NIII/NIV',
           'NIV/NV','OI/OII','OII/OIII','OIII/OIV',
           'OIV/OV','OV/OVI','SIV/SV','SV/SVI',
           'SiII/SiIII','SiIII/SiIV']

ratio_arrays=[[] for i in range(len(hden))]

row=0

for j in ratio_arrays:
    for ion in range(len(data_16[row])):
        if ion<2:
            j.append(data_16[row][ion]/data_16[row][ion+1])
        if ion>2 and ion<6:
            j.append(data_16[row][ion]/data_16[row][ion+1])
        if ion>6 and ion<12:
            j.append(data_16[row][ion]/data_16[row][ion+1])
        if ion>12 and ion<15:
            j.append(data_16[row][ion]/data_16[row][ion+1])
        if ion>15 and ion<18:
            j.append(data_16[row][ion]/data_16[row][ion+1])
    row+=1

#Store the ratios of the ions obtained in the simulation run

ratio_arrays=np.array(ratio_arrays)
ratio_arrays=ratio_arrays.transpose()
ratio_arrays=ratio_arrays.tolist()

c3_c4_q16=ratio_arrays[1]
o3_o4_q16=ratio_arrays[7]

print("Done for 16!!")

#For Q15
import os
import numpy as np
import subprocess
UVB_Q = 15
dirname='/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final'
run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe'
hden= np.arange(-5,-2.999,0.01)
z = -1
stcolden= 14
if not os.path.exists(dirname):
     os.makedirs(dirname)
os.chdir(dirname)
filena = dirname + '/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.in'
f=open(filena,"w+")
f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] \nhden -4 vary \ngrid range from -5 to -3 with 0.01 dex steps \nmetals " + "{:.2f}".format(z)+ " log \nelement helium abundance 0.081632653 linear \nstop column density "+"{:.1f}".format(stcolden)+" neutral H \nconstant temperature, t=1e4 K [linear] \nsave species column density \".spC\" no hash \n\"C+\" \n\"C+2\" \n\"C+3\" \n\"N+\" \n\"N+2\" \n\"N+3\" \n\"N+4\" \n\"O\" \n\"O+\"  \n\"O+2\" \n\"O+3\" \n\"O+4\" \n\"O+5\" \n\"S+3\" \n\"S+4\" \n\"S+5\" \n\"Si+\" \n\"Si+2\" \n\"Si+3\" \nend")
f.close()
process = subprocess.Popen([run, "prog"+"{:.0f}".format(UVB_Q)+"_Oct12.in"], 
stdout =subprocess.PIPE)
process.stdout.read()

from astropy.table import Table
data_15=np.genfromtxt('/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q) +'/Final/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.spC')

colratios=['CII/CIII', 'CIII/CIV', 'NII/NIII', 'NIII/NIV',
           'NIV/NV','OI/OII','OII/OIII','OIII/OIV',
           'OIV/OV','OV/OVI','SIV/SV','SV/SVI',
           'SiII/SiIII','SiIII/SiIV']

ratio_arrays=[[] for i in range(len(hden))]

row=0

for j in ratio_arrays:
    for ion in range(len(data_15[row])):
        if ion<2:
            j.append(data_15[row][ion]/data_15[row][ion+1])
        if ion>2 and ion<6:
            j.append(data_15[row][ion]/data_15[row][ion+1])
        if ion>6 and ion<12:
            j.append(data_15[row][ion]/data_15[row][ion+1])
        if ion>12 and ion<15:
            j.append(data_15[row][ion]/data_15[row][ion+1])
        if ion>15 and ion<18:
            j.append(data_15[row][ion]/data_15[row][ion+1])
    row+=1

#Store the ratios of the ions obtained in the simulation run

ratio_arrays=np.array(ratio_arrays)
ratio_arrays=ratio_arrays.transpose()
ratio_arrays=ratio_arrays.tolist()

c3_c4_q15=ratio_arrays[1]
o3_o4_q15=ratio_arrays[7]

print("Done for 15!!")


#For Q14
import os
import numpy as np
import subprocess
UVB_Q = 14
dirname='/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final'
run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe'
hden= np.arange(-5,-2.999,0.01)
z = -1
stcolden= 14
if not os.path.exists(dirname):
     os.makedirs(dirname)
os.chdir(dirname)
filena = dirname + '/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.in'
f=open(filena,"w+")
f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] \nhden -4 vary \ngrid range from -5 to -3 with 0.01 dex steps \nmetals " + "{:.2f}".format(z)+ " log \nelement helium abundance 0.081632653 linear \nstop column density "+"{:.1f}".format(stcolden)+" neutral H \nconstant temperature, t=1e4 K [linear] \nsave species column density \".spC\" no hash \n\"C+\" \n\"C+2\" \n\"C+3\" \n\"N+\" \n\"N+2\" \n\"N+3\" \n\"N+4\" \n\"O\" \n\"O+\"  \n\"O+2\" \n\"O+3\" \n\"O+4\" \n\"O+5\" \n\"S+3\" \n\"S+4\" \n\"S+5\" \n\"Si+\" \n\"Si+2\" \n\"Si+3\" \nend")
f.close()
process = subprocess.Popen([run, "prog"+"{:.0f}".format(UVB_Q)+"_Oct12.in"], 
stdout =subprocess.PIPE)
process.stdout.read()

from astropy.table import Table
data_14=np.genfromtxt('/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q) +'/Final/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.spC')

#Calculate the ratios of ions from simulation run

colratios=['CII/CIII', 'CIII/CIV', 'NII/NIII', 'NIII/NIV',
           'NIV/NV','OI/OII','OII/OIII','OIII/OIV',
           'OIV/OV','OV/OVI','SIV/SV','SV/SVI',
           'SiII/SiIII','SiIII/SiIV']

ratio_arrays=[[] for i in range(len(hden))]

row=0

for j in ratio_arrays:
    for ion in range(len(data_14[row])):
        if ion<2:
            j.append(data_14[row][ion]/data_14[row][ion+1])
        if ion>2 and ion<6:
            j.append(data_14[row][ion]/data_14[row][ion+1])
        if ion>6 and ion<12:
            j.append(data_14[row][ion]/data_14[row][ion+1])
        if ion>12 and ion<15:
            j.append(data_14[row][ion]/data_14[row][ion+1])
        if ion>15 and ion<18:
            j.append(data_14[row][ion]/data_14[row][ion+1])
    row+=1

#Store the ratios of the ions obtained in the simulation run

ratio_arrays=np.array(ratio_arrays)
ratio_arrays=ratio_arrays.transpose()
ratio_arrays=ratio_arrays.tolist()

c3_c4_q14=ratio_arrays[1]
o3_o4_q14=ratio_arrays[7]

print("Done for 14!!")


#For Q17
import os
import numpy as np
import subprocess
UVB_Q = 17
dirname='/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final'
run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe'
hden= np.arange(-5,-2.999,0.01)
z = -1
stcolden= 14
if not os.path.exists(dirname):
     os.makedirs(dirname)
os.chdir(dirname)
filena = dirname + '/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.in'
f=open(filena,"w+")
f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] \nhden -4 vary \ngrid range from -5 to -3 with 0.01 dex steps \nmetals " + "{:.2f}".format(z)+ " log \nelement helium abundance 0.081632653 linear \nstop column density "+"{:.1f}".format(stcolden)+" neutral H \nconstant temperature, t=1e4 K [linear] \nsave species column density \".spC\" no hash \n\"C+\" \n\"C+2\" \n\"C+3\" \n\"N+\" \n\"N+2\" \n\"N+3\" \n\"N+4\" \n\"O\" \n\"O+\"  \n\"O+2\" \n\"O+3\" \n\"O+4\" \n\"O+5\" \n\"S+3\" \n\"S+4\" \n\"S+5\" \n\"Si+\" \n\"Si+2\" \n\"Si+3\" \nend")
f.close()
process = subprocess.Popen([run, "prog"+"{:.0f}".format(UVB_Q)+"_Oct12.in"], 
stdout =subprocess.PIPE)
process.stdout.read()

from astropy.table import Table
data_17=np.genfromtxt('/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q) +'/Final/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.spC')

#Calculate the ratios of ions from simulation run

colratios=['CII/CIII', 'CIII/CIV', 'NII/NIII', 'NIII/NIV',
           'NIV/NV','OI/OII','OII/OIII','OIII/OIV',
           'OIV/OV','OV/OVI','SIV/SV','SV/SVI',
           'SiII/SiIII','SiIII/SiIV']

ratio_arrays=[[] for i in range(len(hden))]

row=0

for j in ratio_arrays:
    for ion in range(len(data_17[row])):
        if ion<2:
            j.append(data_17[row][ion]/data_17[row][ion+1])
        if ion>2 and ion<6:
            j.append(data_17[row][ion]/data_17[row][ion+1])
        if ion>6 and ion<12:
            j.append(data_17[row][ion]/data_17[row][ion+1])
        if ion>12 and ion<15:
            j.append(data_17[row][ion]/data_17[row][ion+1])
        if ion>15 and ion<18:
            j.append(data_17[row][ion]/data_17[row][ion+1])
    row+=1

#Store the ratios of the ions obtained in the simulation run

ratio_arrays=np.array(ratio_arrays)
ratio_arrays=ratio_arrays.transpose()
ratio_arrays=ratio_arrays.tolist()

c3_c4_q17=ratio_arrays[1]
o3_o4_q17=ratio_arrays[7]

print("Done for 17!!")

#For Q19
import os
import numpy as np
import subprocess
UVB_Q = 19
dirname='/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q)+'/Final'
run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe'
hden= np.arange(-5,-2.999,0.01)
z = -1
stcolden= 14
if not os.path.exists(dirname):
     os.makedirs(dirname)
os.chdir(dirname)
filena = dirname + '/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.in'
f=open(filena,"w+")
f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] \nhden -4 vary \ngrid range from -5 to -3 with 0.01 dex steps \nmetals " + "{:.2f}".format(z)+ " log \nelement helium abundance 0.081632653 linear \nstop column density "+"{:.1f}".format(stcolden)+" neutral H \nconstant temperature, t=1e4 K [linear] \nsave species column density \".spC\" no hash \n\"C+\" \n\"C+2\" \n\"C+3\" \n\"N+\" \n\"N+2\" \n\"N+3\" \n\"N+4\" \n\"O\" \n\"O+\"  \n\"O+2\" \n\"O+3\" \n\"O+4\" \n\"O+5\" \n\"S+3\" \n\"S+4\" \n\"S+5\" \n\"Si+\" \n\"Si+2\" \n\"Si+3\" \nend")
f.close()
process = subprocess.Popen([run, "prog"+"{:.0f}".format(UVB_Q)+"_Oct12.in"], 
stdout =subprocess.PIPE)
process.stdout.read()

from astropy.table import Table
data_19=np.genfromtxt('/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q) +'/Final/prog'+'{:.0f}'.format(UVB_Q)+'_Oct12.spC')

#Calculate the ratios of ions from simulation run

colratios=['CII/CIII', 'CIII/CIV', 'NII/NIII', 'NIII/NIV',
           'NIV/NV','OI/OII','OII/OIII','OIII/OIV',
           'OIV/OV','OV/OVI','SIV/SV','SV/SVI',
           'SiII/SiIII','SiIII/SiIV']

ratio_arrays=[[] for i in range(len(hden))]

row=0

for j in ratio_arrays:
    for ion in range(len(data_19[row])):
        if ion<2:
            j.append(data_19[row][ion]/data_19[row][ion+1])
        if ion>2 and ion<6:
            j.append(data_19[row][ion]/data_19[row][ion+1])
        if ion>6 and ion<12:
            j.append(data_19[row][ion]/data_19[row][ion+1])
        if ion>12 and ion<15:
            j.append(data_19[row][ion]/data_19[row][ion+1])
        if ion>15 and ion<18:
            j.append(data_19[row][ion]/data_19[row][ion+1])
    row+=1

#Store the ratios of the ions obtained in the simulation run

ratio_arrays=np.array(ratio_arrays)
ratio_arrays=ratio_arrays.transpose()
ratio_arrays=ratio_arrays.tolist()

c3_c4_q19=ratio_arrays[1]
o3_o4_q19=ratio_arrays[7]

print("Done for 19!!")
