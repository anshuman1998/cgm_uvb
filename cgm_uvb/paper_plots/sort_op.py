import numpy as np
import astropy.table as tab

#path  = '/home/vikram/cloudy_run/diff_op/photo_NH15'
path  = '/home/vikram/cloudy_run/diff_op/hybrid_NH15'


uvb = ['KS18', 'HM12',  'P19', 'FG20']
uvb_Q = [14, 15, 16, 17, 18, 19, 20]

uvb_models =[]
the_Q_values = []
for background in uvb:
    if background == 'KS18':
        for q in uvb_Q:
            uvb_models.append(background)
            the_Q_values.append(q)
    else:
        q = 18
        uvb_models.append(background)
        the_Q_values.append(q)



# combine files
#-- imp: the structure and columns are identical

#---- read the identical columns and store them
q=18
uvb = 'KS18'
#file_op = path + '/photo_model_{}_Q{}_lsf_inference.fits'.format(uvb, q)
logT= 5.5
file_op = path + '/hybrid_model_{}_Q{}_logT{:.0f}_lsf_inference.fits'.format(uvb, q, logT*100)


data_common = tab.Table.read(file_op)

new_tab = tab.Table()

new_tab.add_columns([data_common['true_uvb'], data_common['true_Q'], data_common['true_nH'], data_common['true_Z']],
                    names = ('true_uvb', 'true_Q', 'true_nH', 'true_logZ'))

#----add last two columns later
ncolname_array0 =[]
zcolname_array0 = []
for uvb, q in zip(uvb_models, the_Q_values):
#    file_op = path + '/photo_model_{}_Q{}_lsf_inference.fits'.format(uvb, q)
    file_op = path + '/hybrid_model_{}_Q{}_logT{:.0f}_lsf_inference.fits'.format(uvb, q, logT*100)

    data = tab.Table.read(file_op)
    print(len(data), uvb, q)
    col_name_nH = 'nH_'+ uvb + '_Q{}'.format(q)
    ncolname_array0.append(col_name_nH)
    col_name_Z = 'Z_'+ uvb + '_Q{}'.format(q)
    zcolname_array0.append(col_name_Z)
    # note both are in log
    new_tab.add_columns([data['nH'], data['Z']], names = (col_name_nH, col_name_Z))


#----two name array

print('now adding dmax')
# delta max calculations
n_dmax = [] #----for all 9 => 7*KS + FG20 +P19
z_dmax = []
n_dmax_KS = []
z_dmax_KS = []
n_dmax_HM = []
z_dmax_HM = []
print(len(data_common), '===> this no')
for i in range(len(data_common)):
    ncolname_array = ncolname_array0 [:]
    zcolname_array = zcolname_array0 [:]

    n_dmax_HM.append(max([new_tab[k][i] for k in ncolname_array]) - np.min([new_tab[k][i] for k in ncolname_array]))
    z_dmax_HM.append(max([new_tab[k][i] for k in zcolname_array]) - np.min([new_tab[k][i] for k in zcolname_array]))

    ncolname_array.remove('nH_HM12_Q18')
    zcolname_array.remove('Z_HM12_Q18')
    n_dmax.append(max([new_tab[k][i] for k in ncolname_array]) - np.min([new_tab[k][i] for k in ncolname_array]))
    z_dmax.append(max([new_tab[k][i] for k in zcolname_array]) - np.min([new_tab[k][i] for k in zcolname_array]))

    ncolname_array.remove('nH_P19_Q18')
    zcolname_array.remove('Z_P19_Q18')
    ncolname_array.remove('nH_FG20_Q18')
    zcolname_array.remove('Z_FG20_Q18')
    n_dmax_KS.append(max([new_tab[k][i] for k in ncolname_array]) - np.min([new_tab[k][i] for k in ncolname_array]))
    z_dmax_KS.append(max([new_tab[k][i] for k in zcolname_array]) - np.min([new_tab[k][i] for k in zcolname_array]))

new_tab.add_columns([n_dmax_HM, z_dmax_HM, n_dmax, z_dmax, n_dmax_KS, z_dmax_KS],
                    names=('n_dmax_HM', 'z_dmax_HM', 'n_dmax', 'z_dmax', 'n_dmax_KS', 'z_dmax_KS'))

new_tab.add_columns([data_common['n_ions'], data_common['ions']], names = ('n_ions', 'ions'))

final_filename = path+ '/all_combined_logT550.fits'
new_tab.write(final_filename, overwrite = True)
