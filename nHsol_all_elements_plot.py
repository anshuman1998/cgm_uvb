%matplotlib qt
import matplotlib
import matplotlib.pyplot as plt

elements=['Carbon','Nitrogen','Oxygen','Sulphur','Silicon']
plt.scatter(elements,mean_nH_Q20,color='crimson', label='UVB=Q20')
plt.errorbar(elements,mean_nH_Q20,yerr=stdev_nH_Q20,fmt='o',color='crimson')

plt.scatter(elements,mean_nH_Q19,color='goldenrod', label='UVB=Q19')
plt.errorbar(elements,mean_nH_Q19,yerr=stdev_nH_Q19,fmt='o',color='goldenrod')

plt.scatter(elements,mean_nH_Q18,color='dodgerblue', label='UVB=Q18')
plt.errorbar(elements,mean_nH_Q18,yerr=stdev_nH_Q18,fmt='o',color='dodgerblue')

plt.scatter(elements,mean_nH_Q17,color='deeppink', label='UVB=Q17')
plt.errorbar(elements,mean_nH_Q17,yerr=stdev_nH_Q17,fmt='o',color='deeppink')

plt.scatter(elements,mean_nH_Q16,color='forestgreen', label='UVB=Q16')
plt.errorbar(elements,mean_nH_Q16,yerr=stdev_nH_Q16,fmt='o',color='forestgreen')

plt.scatter(elements,mean_nH_Q15,color='darkorange', label='UVB=Q15')
plt.errorbar(elements,mean_nH_Q15,yerr=stdev_nH_Q15,fmt='o',color='darkorange')

plt.scatter(elements,mean_nH_Q14,color='mediumorchid', label='UVB=Q14')
plt.errorbar(elements,mean_nH_Q14,yerr=stdev_nH_Q14,fmt='o',color='mediumorchid')

plt.axhline(y=-4,color='black')
plt.legend(title='log(Hydrogen Density) Solutions for various Elements',loc='upper center',bbox_to_anchor=(0.5,1.1), ncol = 4,fancybox=True, shadow=True)
plt.ylabel('log(nH Solution)')
plt.xlabel('Elements')
plt.show(block=False)

