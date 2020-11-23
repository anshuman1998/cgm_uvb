import astropy.table as tab
import numpy as np

fg20  = tab.Table.read('fg20_spec_nu.dat', format = 'ascii')
fg20.remove_row(0) # removing z row (first value of which is garbage)

fg20_orig  = tab.Table.read('fg20_spec_nu.dat', format = 'ascii')

folder_to_store = '/home/vikram/cloudy_run/scratch/fg_20_fits_files'

col = fg20_orig.colnames
col.pop(0)
z_array = []
for name in col:
    z_array.append(fg20_orig[name][0])


for iz, colname in zip(z_array, col):
    redshift = iz
    new_file_name = folder_to_store + '/FG20_EBL_z_{:.2f}.fits'.format(redshift)
    uvb_table =  tab.Table(meta= {'z': redshift,
                                  'ene' : 'Ryd',
                                  'EBL' : 'erg/s/cm^2/Hz/Sr',
                                  'info' : 'FG20 rescaled model'})
    uvb_table.add_column(fg20['col1'], name = 'ene')
    uvb_table.add_column(1e-21 *fg20[colname], name = 'Jnu')

    uvb_table.write(new_file_name, overwrite = True)
    print(new_file_name)

