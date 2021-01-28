import astropy.table as tab
import numpy as np
from cgm_uvb.find_uvb_scaling import find_uvb_scaling

def find_max_numbers(filename):
    print('for file: ', filename)

    data =  tab.Table.read(filename)
    # check if files are in log
    if data ['nH'][0] < -1:
        log = True
    else:
        log = False

    max_n_mod =  data['uvb'][np.argmax(data['nH'])]
    min_n_mod =  data['uvb'][np.argmin(data['nH'])]
    if log :
        diff_n_all  = np.max(data ['nH']) - np.min(data ['nH'])
        diff_n_ks  = data ['nH'][0] - data ['nH'][6]
    else:
        print('-- lin')
        data['nH'] = np.log10(data['nH'])
        diff_n_all  = np.max(data ['nH']) - np.min(data ['nH'])
        diff_n_ks  = data ['nH'][0] - data ['nH'][6]

    print('delta nH = ', diff_n_all, 'for models:', max_n_mod, '-', min_n_mod)
    print('delta nH = ', diff_n_ks, 'for KS models only:',  data['uvb'][0], '-', data['uvb'][6])



    max_n_mod =  data['uvb'][np.argmax(data['Z'])]
    min_n_mod =  data['uvb'][np.argmin(data['Z'])]
    if log :
        diff_n_all  = np.max(data ['Z']) - np.min(data ['Z'])
        diff_n_ks  = data ['Z'][0] - data ['Z'][6]

    else:
        data['Z'] = np.log10(data['Z'])
        diff_n_all  = np.max(data ['Z']) - np.min(data ['Z'])
        diff_n_ks  = data ['Z'][0] - data ['Z'][6]

    print('delta Z = ', diff_n_all, 'for models:', max_n_mod, '-', min_n_mod)
    print('delta Z = ', diff_n_ks, 'for KS models only:', data['uvb'][0], '-', data['uvb'][6])

    return



print('-------------------for the photoinized model all uvb')

lsf_file_name = '/home/vikram/cloudy_run/figures/2D/NH15_log_lsf_out.fits'
find_max_numbers(lsf_file_name)

print('-------------------for the hybrid model all uvb')

lsf_file_name = '/home/vikram/cloudy_run/figures/hybrid/NH15_log_lsf_hybrid_T550.fits'
find_max_numbers(lsf_file_name)

print('-------------------for the photoinized model with rescaled uvb')
lsf_file_name = '/home/vikram/cloudy_run/figures/rescaled/rescaled_NH15_log_lsf_out.fits'
find_max_numbers(lsf_file_name)

print('-------------------for the hybrid model with rescaled uvb')
lsf_file_name = '/home/vikram/cloudy_run/figures/rescaled_hybrid/NH15_log_lsf_hybrid_T550.fits'
find_max_numbers(lsf_file_name)



#-------------uvb scaling
print('--------------- for UVB scaling')
uvb_Q = [14, 15, 16, 17, 18, 19, 20]

for q_model in uvb_Q:
    scaling_factor  =  find_uvb_scaling(uvb = 'KS18', uvb_Q = q_model)
    print(scaling_factor, np.log10(scaling_factor), q_model, 'KS19')

for i in ['HM12',  'P19', 'FG20']:
    scaling_factor_uvb = find_uvb_scaling(uvb=i)
    print(scaling_factor_uvb, np.log10(scaling_factor_uvb), i)
