import astropy.table as tab
import os
import glob

basepath =  os.getcwd()

uvb_array= [14, 15, 16, 17, 18, 19, 20]

file_names = tab.Table()
for Q in uvb_array:
    q_name = 'Q{}'.format(Q)
    dir_path =  basepath + '/' + q_name
    all_files = sorted(glob.glob(dir_path + '/*.txt'))
    file_names.add_column(all_files, name = q_name)

trial_path =  basepath + '/Q18'
fname = sorted(glob.glob(trial_path + '/*.txt'))
z_array = []
for i in fname :
    t =  (i.split('z_')[-1]).split('.txt')[0]
    z_array.append(float(t))

file_names.add_column(z_array, name = 'z')

for iz in range(len(z_array)):
    redshift = file_names['z'][iz]
    new_file_name = basepath + '/KS19_EBL_z_{:.1f}.fits'.format(redshift)
    uvb_table =  tab.Table(meta= {'z': redshift,
                                  'Wave' : 'Angstrom',
                                  'EBL' : 'erg/s/cm^2/Hz/Sr',
                                  'info' : 'KS19 EBL models Q18 is fiducial'})
    for Q in uvb_array:
        q_name = 'Q{}'.format(Q)
        file_to_copy = file_names[q_name][iz]
        uvb_to_copy  = tab.Table.read(file_to_copy, format = 'ascii')
        if Q == uvb_array[0]:
            uvb_table.add_column(uvb_to_copy['col1'], name = 'Wave' )

        uvb_table.add_column(uvb_to_copy['col2'], name =q_name)

    uvb_table.write(new_file_name, overwrite = True)

