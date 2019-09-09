import os
import numpy as np
import subprocess

from astropy.table import Table

#Run Cloudy to generate test observations

UVB_Q = 18 #Q18 Quasar Source
dirname='/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q)

run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe'
hden= -4
z = 0.1
stcolden= 14
if not os.path.exists(dirname):
     os.makedirs(dirname)
os.chdir(dirname)
filena = dirname + '/prog'+'{:.0f}'.format(UVB_Q)+'.in'
f=open(filena,"w+")
f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] 
        \nhden" + "{:.0f}".format(hden)+ 
        "\nmetals" + "{:.2f}".format(z)+ " log 
        \nelement helium abundance 0.081632653 linear 
        \nstop column density"+"{:.1f}".format(stcolden)+" neutral H 
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

#Saving column density values for 20 metals

f.close()
process = subprocess.Popen([run, "prog"+"{:.0f}".format(UVB_Q)+".in"], 
                           stdout =subprocess.PIPE)
process.stdout.read()

data=np.genfromtxt('/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+
'{:.0f}'.format(UVB_Q)+'/prog'+'{:.0f}'.format(UVB_Q)+'.spC')

#Save Column density data on astropy table

t=Table(data,names=['CII','CIII','CIV','CV','MgII','NeVIII',
'NII','NIII','NIV','OII','OIII','OIV','OV','OVI','SiII',
'SiIII','SiIV','SIV','SV','SVI'],
meta={'UVB_Q':UVB_Q, 'metallicity': z,'n_H':hden,'log(stop.colden)':stcolden})

t.write('/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q)+
        '/colden.fits')
read_obs=Table.read('/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+
'{:.0f}'.format(UVB_Q)+'/colden.fits')
read_obs.meta
