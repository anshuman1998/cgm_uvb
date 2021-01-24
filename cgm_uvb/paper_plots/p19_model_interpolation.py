import numpy as np
from scipy.interpolate import interp1d
import astropy.table as tab

z  =  np.genfromtxt("bkgthick.out", skip_header= 10, max_rows=1)
data  =  np.loadtxt("bkgthick.out", skiprows = 11)
wave = data[:, 0]
data_new =  data[:, 1:]

# remove repeated wavelengths
wave_sorted, inds = np.unique(wave, return_index =  True)
sorted_data  = np.zeros((len(inds), np.shape(data_new)[1]))

for j in range(len(inds)):
    sorted_data[j, :] = data_new[ inds[j], :]

# now do the interpolation
f = interp1d(z, sorted_data)

my_z = np.arange(0.1, 1.1, 0.2)

folder_to_store = '/home/vikram/cgm_uvb/cgm_uvb/paper_plots/p19_ebl'
for z_val in my_z:
    redshift = z_val
    new_file_name = folder_to_store + '/P19_EBL_z_{:.2f}.fits'.format(redshift)
    uvb_table =  tab.Table(meta= {'z': redshift,
                                  'Wave' : 'Angstrom',
                                  'Jnu' : 'erg/s/cm^2/Hz/Sr',
                                  'info' : 'FG20 rescaled model'})
    uvb_table.add_column(wave_sorted, name = 'Wave')
    uvb_table.add_column(f(z_val), name = 'Jnu')

    uvb_table.write(new_file_name, overwrite = True)
    print(new_file_name)






"""
ipython terminal checks
import numpy as np
z  =  np.genfromtxt("bkgthick.out", skip_header= 10, max_rows=1)
z
data  =  np.loadtxt("bkgthick.out", skiprows = 11)
wave = data[:, 0]
wave
wave.shape
data_new =  data[:, 1:]
z.shape
data_new.shape
from scipy.interpolate import interp1d
f = interp1d(z, data_new, axis = 0)
f = interp1d(z, data_new)
ebl0 = data_new[:, 0]
ebl0.shape
z[0]
z[1]
f(0.04912)
data_new[:, 1] -  f(0.04912)
z[2] +z[3]
(z[2] +z[3])/2
ebl  = (data_new[:, 2] +  data_new[:, 3])/2
f(0.12765) - ebl
f(0.1276499999) - ebl

"""