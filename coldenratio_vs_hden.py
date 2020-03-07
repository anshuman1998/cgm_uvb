#Plotting column density ratio CIII/CIV and OIII/OIV versus hydrogen density

import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':18})

alpha = [-2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4]
plt.plot(hden,np.log10(c3_c4_q14),color='palegreen')
plt.plot(hden,np.log10(o3_o4_q14),color='lightskyblue',linestyle='--')

plt.plot(hden,np.log10(c3_c4_q15),color='lightgreen')
plt.plot(hden,np.log10(o3_o4_q15),color='deepskyblue',linestyle='--')

plt.plot(hden,np.log10(c3_c4_q16),color='limegreen')
plt.plot(hden,np.log10(o3_o4_q16),color='dodgerblue',linestyle='--')

plt.plot(hden,np.log10(c3_c4_q17),color='forestgreen')
plt.plot(hden,np.log10(o3_o4_q17),color='royalblue',linestyle='--')

plt.plot(hden,np.log10(c3_c4_q18),color='green')
plt.plot(hden,np.log10(o3_o4_q18),color='blue',linestyle='--')

plt.plot(hden,np.log10(c3_c4_q19),color='darkgreen')
plt.plot(hden,np.log10(o3_o4_q19),color='mediumblue',linestyle='--')

plt.plot(hden,np.log10(c3_c4_q20),color='darkolivegreen')
plt.plot(hden,np.log10(o3_o4_q20),color='navy',linestyle='--')

plt.axvline(x=-4,color='crimson',label='Test Observation Values')

a= [None]
b=[None]
plt.plot(b,a,color='lightgrey',label='UVB=Q14')
plt.plot(b,a,color='black',label='UVB=Q20')
plt.plot(b,a,color='green',label='CIII/CIV',marker='*')
plt.plot(b,a,color='blue',label='OIII/OIV',marker='*',linestyle='--')
c3_c4_obs = np.log10(1.1233751021002576)
o3_o4_obs = np.log10(0.3158421652224946)
plt.scatter(-4,c3_c4_obs,color='green',marker='*',s=80)
plt.scatter(-4,o3_o4_obs,color='blue',marker='*',s=80)

plt.legend(title='CIII/CIV and OIII/OIV for different UVB',loc='upper center',bbox_to_anchor=(0.5,1.17), ncol = 2,fancybox=True, shadow=True)
plt.ylabel('log(Column Density ratio)')
plt.xlabel('log(Hydrogen Density)')

plt.show(block=False)

