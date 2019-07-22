#Observations using Q18, with log(nH)=-4 and log(Z)=0.1

import osimport numpy as npimport subprocessUVB_Q = 18dirname='/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q)run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe'hden= -4z = 0.1stcolden= 14if not os.path.exists(dirname):     os.makedirs(dirname)os.chdir(dirname)filena = dirname + '/prog'+'{:.0f}'.format(UVB_Q)+'.in'f=open(filena,"w+")f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] \nhden" + "{:.0f}".format(hden)+ "\nmetals" + "{:.2f}".format(z)+ " log \nelement helium abundance 0.081632653 linear \nstop column density"+"{:.1f}".format(stcolden)+"neutral H \nsave species column density \".spC\" no hash \n\"C+\" \n\"C+2\" \n\"C+3\" \n\"C+4\" \n\"Mg+\" \n\"Ne+7\" \n\"N+\" \n\"N+2\" \n\"N+3\" \n\"O+\" \n\"O+2\" \n\"O+3\" \n\"O+4\" \n\"O+5\" \n\"Si+\" \n\"Si+2\" \n\"Si+3\" \n\"S+3\" \n\"S+4\" \n\"S+5\" \nend")f.close()process = subprocess.Popen([run, "prog"+"{:.0f}".format(UVB_Q)+".in"], stdout =subprocess.PIPE)process.stdout.read()


#for Q20
UVB_Q = 20
dirname='/Users/anshumanacharya/Downloads/c17.01/source/pyprog'+'{:.0f}'.format(UVB_Q) 
run='/Users/anshumanacharya/Downloads/c17.01/source/cloudy.exe' 
stcolden= 14 
minstddev=1.0
hden=np.arange(-4.2,-4.1,0.01)   #Decide the grid of nH
if not os.path.exists(dirname):      
    os.makedirs(dirname) 
os.chdir(dirname)   
z=np.arange(0.01,1.02,0.05)
metal=[]
print("Begin Iterative Program")
for i in z:
        metal.append(10 ** i)  #Converting to linear metallicity values
name=['CII','CIII','CIV','CV','MgII','NeVIII','NII','NIII','NIV','OII','OIII','OIV','OV','OVI','SiII','SiIII','SiIV',
          'SIV','SV','SVI'] #Ions for which calculation is performed.
for h in hden:
    filena = dirname + '/prog_metal'+'{:.0f}'.format(UVB_Q)+'_nH_'+'{:.2f}'.format(h)+'.in'     
    f=open(filena,"w+")     
    f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q="+"{:.0f}".format(UVB_Q)+"] \nhden "+"{:.2f}".format(h)+" \nmetals 0.11 log vary \ngrid range 0.01 to 1.01 with 0.05 dex steps \nelement helium abundance 0.081632653 linear \nstop column density "+"{:.1f}".format(stcolden)+" neutral H \nconstant temperature, t=1e4 K [linear] \nsave species column density \".spC\" no hash \n\"C+\" \n\"C+2\" \n\"C+3\" \n\"C+4\" \n\"Mg+\" \n\"Ne+7\" \n\"N+\" \n\"N+2\" \n\"N+3\" \n\"O+\" \n\"O+2\" \n\"O+3\" \n\"O+4\" \n\"O+5\" \n\"Si+\" \n\"Si+2\" \n\"Si+3\" \n\"S+3\" \n\"S+4\" \n\"S+5\" \nend")     
    f.close()     
    process = subprocess.Popen([run, "prog_metal"+"{:.0f}".format(UVB_Q)+'_nH_'+'{:.2f}'.format(h)+'.in'],stdout =subprocess.PIPE)     
    process.stdout.read()

    trial = np.genfromtxt("/Users/anshumanacharya/Downloads/c17.01/source/pyprog20/prog_metal20_nH_"+'{:.2f}'.format(h)+".spC")
    obs_vals=[0.452850,22.434700,19.903900,17.462100,0.001982,0.002153,0.076940,7.407480,9.796510,0.297943,23.139700,73.066600,20.292600,4.076360,0.010686,
              0.255324,1.030160,1.687530,0.517160,0.889996] #Observed column densities of the 20 ions divided by 10^12 for log(Z)=0.1,log(nH)=-4
    
    #Finish running the Cloudy program for log(metallicity) varying from 0.01 to 1.01 for each nH value
    print("Ran for n_H= ", h)
    lin_met=[]
    for i in range(len(name)):
        m,c= np.polyfit(metal,(trial[:,i]/(10 ** 12)),1)
        lin_met.append((obs_vals[i]-c)/m)  #Find the required metallicity value to match with observed data
        
    lin_met=np.array(lin_met)
    if lin_met.all()>0:               #Consider only those cases which are physical, i.e., linear metallicity is positive
        stdev=np.std(lin_met)
        print(stdev)
        if stdev<minstddev:         #Finding the nH value where the standard deviation of req. metallicity is lowest.
            minstddev=stdev
            metallicity=lin_met
            mean=np.mean(lin_met)
            n_H=h         
name.append("Mean Z")
name.append("Stddev Z")
metallicity =metallicity.tolist()
metallicity.append(mean)
metallicity.append(minstddev)
#Save the result
from astropy.table import Table
soln_Q20=Table([name,metallicity],names=('Ion', 'Metallicity'))
soln_Q20.write("/Users/anshumanacharya/Downloads/c17.01/source/pyprog20/soln_Q20_nH_"+'{:.2f}'.format(n_H)+".txt", format='ascii.tab',overwrite='True')
print("Done!")