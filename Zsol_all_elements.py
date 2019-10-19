%matplotlib qt
import matplotlib
import matplotlib.pyplot as plt

elements=['Carbon','Nitrogen','Oxygen','Sulphur','Silicon']

plt.plot(elements,mean_Z_Q20,color='crimson', label='UVB=Q20',marker='o')
plt.errorbar(elements,mean_Z_Q20,yerr=stdev_Z_Q20,fmt='o',color='crimson', alpha=0.5,capsize=10)

plt.plot(elements,mean_Z_Q19,color='goldenrod', label='UVB=Q19',marker='o')
plt.errorbar(elements,mean_Z_Q19,yerr=stdev_Z_Q19,fmt='o',color='goldenrod', alpha=0.5,capsize=10)

plt.plot(elements,mean_Z_Q18,color='dodgerblue', label='UVB=Q18',marker='o', alpha=0.5)
plt.errorbar(elements,mean_Z_Q18,yerr=0.0,fmt='o',color='dodgerblue')

plt.plot(elements,mean_Z_Q17,color='deeppink', label='UVB=Q17',marker='o')
plt.errorbar(elements,mean_Z_Q17,yerr=stdev_Z_Q17,fmt='o',color='deeppink', alpha=0.5,capsize=10)

plt.plot(elements,mean_Z_Q16,color='forestgreen', label='UVB=Q16',marker='o')
plt.errorbar(elements,mean_Z_Q16,yerr=stdev_Z_Q16,fmt='o',color='forestgreen', alpha=0.5,capsize=10)

plt.plot(elements,mean_Z_Q15,color='darkorange', label='UVB=Q15',marker='o')
plt.errorbar(elements,mean_Z_Q15,yerr=stdev_Z_Q15,fmt='o',color='darkorange', alpha=0.5,capsize=10)

plt.plot(elements,mean_Z_Q14,color='mediumorchid', label='UVB=Q14',marker='o')
plt.errorbar(elements,mean_Z_Q14,yerr=stdev_Z_Q14,fmt='o',color='mediumorchid',alpha=0.5,capsize=10)

plt.axhline(y=(10**0.1),color='black')
plt.legend(title='Metallicity Solutions for various Elements',loc='upper center',bbox_to_anchor=(0.5,1.1), 
           ncol = 4,fancybox=True, shadow=True)

plt.ylabel('Metallicity')
plt.xlabel('Elements')
plt.show(block=False)
