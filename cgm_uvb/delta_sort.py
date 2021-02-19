import numpy as np
import astropy.table as tab

def find_med(true_uvb, true_Q = 18, true_nH = 1e-4, true_logZ = -1, nion = 8,
             table_file ='/home/vikram/cloudy_run/diff_op/photo_NH15/all_combined.fits' ):

    ks_nH = []
    ks_Z = []
    all_nH = []
    all_Z = []

    data = tab.Table.read(table_file)

    if data != None:
        data = data [data ['n_ions'] == nion]

    data_sort = data [data['true_uvb'] == true_uvb]
    data_sort = data_sort [data_sort['true_Q'] == true_Q]
    data_sort = data_sort [data_sort['true_logZ'] == true_logZ]
    data_sort = data_sort [data_sort['true_nH'] == true_nH]
    print(len(data_sort), ': number of points')


    uvb_array = ['KS18',  'P19', 'FG20']
    uvb_Q = [14, 15, 16, 17, 18, 19, 20]
    uvb_models = []
    the_Q_values = []
    for background in uvb_array:
        if background == 'KS18':
            for q in uvb_Q:
                uvb_models.append(background)
                the_Q_values.append(q)
        else:
            q = 18
            uvb_models.append(background)
            the_Q_values.append(q)



    for uvb, q_val in zip(uvb_models, the_Q_values):
        col_name_nH = 'nH_' + uvb + '_Q{}'.format(q_val)
        col_name_Z = 'Z_' + uvb + '_Q{}'.format(q_val)

        if uvb == 'KS18':
            ks_nH.append(np.median(data_sort[col_name_nH]))
            all_nH.append(np.median(data_sort[col_name_nH]))
            ks_Z.append(np.median(data_sort[col_name_Z]))
            all_Z.append(np.median(data_sort[col_name_Z]))
        else:
            all_nH.append(np.median(data_sort[col_name_nH]))
            all_Z.append(np.median(data_sort[col_name_Z]))

    del_n_ks = np.max(ks_nH) - np.min(ks_nH)
    del_z_ks = np.max(ks_Z) - np.min(ks_Z)
    del_n_all = np.max(all_nH) - np.min(all_nH)
    del_z_all = np.max(all_Z) - np.min((all_Z))

    print('KS', del_n_ks, del_z_ks)
    print('ALL', del_n_all, del_z_all)


    return all_nH, all_Z


def find_med_all( true_nH = 1e-4, true_logZ = -1,
             table_file ='/home/vikram/cloudy_run/diff_op/photo_NH15/all_combined.fits' ):



    uvb_array = ['KS18',  'P19', 'FG20']
    uvb_Q = [14, 15, 16, 17, 18, 19, 20]
    uvb_models = []
    the_Q_values = []
    for background in uvb_array:
        if background == 'KS18':
            for q in uvb_Q:
                uvb_models.append(background)
                the_Q_values.append(q)
        else:
            q = 18
            uvb_models.append(background)
            the_Q_values.append(q)

    nion_array = [ 8]

    data = tab.Table.read(table_file)
    data = data [data['true_logZ'] == true_logZ]
    data = data [data ['true_nH'] == true_nH]

    full_z = []
    full_n = []

    del_n_array = []
    del_z_array = []

    for true_uvb, true_Q in zip(uvb_models, the_Q_values):
        #data_sort = tab.Table()
        data_sort_in = data[data['true_uvb'] == true_uvb]
        data_sort_in = data_sort_in[data_sort_in['true_Q'] == true_Q]

        for nion in nion_array:
            data_sort = data_sort_in[data_sort_in['n_ions'] == nion]

            ks_nH = []
            ks_Z = []
            all_nH = []
            all_Z = []


            #print(len(data_sort), ': number of points')

            for uvb, q_val in zip(uvb_models, the_Q_values):
                col_name_nH = 'nH_' + uvb + '_Q{}'.format(q_val)
                col_name_Z = 'Z_' + uvb + '_Q{}'.format(q_val)

                full_n.append(np.median(data_sort[col_name_nH]))
                full_z.append(np.median(data_sort[col_name_Z]))

                all_nH.append(np.median(data_sort[col_name_nH]))
                all_Z.append(np.median(data_sort[col_name_Z]))

                if uvb == 'KS18':
                    ks_nH.append(np.median(data_sort[col_name_nH]))
                    ks_Z.append(np.median(data_sort[col_name_Z]))


            del_n_ks = np.max(ks_nH) - np.min(ks_nH)
            del_z_ks = np.max(ks_Z) - np.min(ks_Z)
            del_n_all = np.max(all_nH) - np.min(all_nH)
            del_z_all = np.max(all_Z) - np.min((all_Z))
            del_n_array.append(del_n_all)
            del_z_array.append(del_z_all)

            #print('for true', true_uvb, true_Q, 'KS', del_n_ks, del_z_ks, 'nion:', nion)
            print('for true', true_uvb, true_Q, 'AL', del_n_all, del_z_all, 'nion:', nion)

    return full_n, full_z, del_n_array, del_z_array



def make_plot():

    return 


#a , b = find_med('KS18')
#for n in list([1e-3, 1e-4, 1e-5]):
#print(n)import matplotlib.pyplot as plt

#nn = 1e-4
#n, z= find_med_all(true_nH= nn, true_logZ=0)

#------------------ use following code for the work
nH = [1e-5, 1e-4, 1e-3]
lgZ = [-2, -1, 0]

for nn in nH:
    for zz in lgZ:
        n, z, dn, dz = find_med_all(true_nH=nn, true_logZ=zz)

        x = np.median((np.array(dn) / 2))
        y = np.median((np.array(dz) / 2))
        plt.errorbar([np.log10(nn)], [zz], [x], [y], alpha=0.5)
plt.show()



