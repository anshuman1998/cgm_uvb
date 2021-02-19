from cgm_uvb.find_uvb_scaling import get_HI_photoionization_rate
from cgm_uvb.find_uvb_scaling import get_uvb
import numpy as np
import astropy.table as tab
store_path = '/home/vikram/cgm_uvb/cgm_uvb/paper_plots/gamma_HI_files'

uvb = ['KS18', 'HM12',  'P19', 'FG20']
uvb_Q = [14, 15, 16, 17, 18, 19, 20]

my_z = np.arange(0.0, 1.01, 0.1)

for model in uvb:
    if model == 'KS18':
        for q in uvb_Q:
            fname  = 'Gamma_HI_KS18_Q{}.fits'.format(q)
            filename  = store_path + '/' + fname
            gamma_HI_array = []
            for z in my_z:
                w, j, ih_assumed = get_uvb(uvb=model, Q=q, z=z)
                gamma_HI = get_HI_photoionization_rate(w, j, IH_assumed=ih_assumed)
                gamma_HI_array.append(gamma_HI)

            # writing table
            g_table = tab.Table(meta={'uvb': model + 'Q{}'.format(q),
                                      'z': 'redshift',
                                      'g': '1/s'})
            g_table.add_column(my_z, name='z')
            g_table.add_column(gamma_HI_array, name='g')

            g_table.write(filename, overwrite=True)
            print('done', model, 'Q', q)

    else:
        fname = 'Gamma_HI_{}.fits'.format(model)
        filename = store_path + '/' + fname
        gamma_HI_array = []
        for z in my_z:
            w, j, ih_assumed = get_uvb(uvb=model, z=z)
            gamma_HI = get_HI_photoionization_rate(w, j, IH_assumed=ih_assumed)
            gamma_HI_array.append(gamma_HI)

        # writing table
        g_table = tab.Table(meta={'uvb': model,
                                  'z': 'redshift',
                                  'g': '1/s'})
        g_table.add_column(my_z, name='z')
        g_table.add_column(gamma_HI_array, name='g')

        g_table.write(filename, overwrite=True)
        print('done', model)







