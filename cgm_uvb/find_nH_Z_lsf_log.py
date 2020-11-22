import numpy as np
import astropy.table as tab


def get_true_model(model_q, Q= 18, nH = 1e-4):
    """
    :param model: The data where Q18 model is stored
    :return: a row of ion column densities at n_H = 1e-4 cm^-2
    """
    model = model_q.split('_Q')[0] + '_Q{}.fits'.format(Q)
    data = tab.Table.read(model)
    true_ion_col = data [data['hden'] == nH]
   # print(true_ion_col)

    return true_ion_col

def find_density_LSF(model_Q, *ions_to_use):
    model = tab.Table.read(model_Q)
    obs_ion_col = get_true_model(model_Q)

    # filter to use specific ions
    least_square_array = np.zeros(len(model))
    hden_array = np.array(model['hden'])
    for ion in ions_to_use:
        least_square_array += (model[ion] - obs_ion_col[ion])**2

    print(np.min(least_square_array))
    print(hden_array[np.argmin(least_square_array)])

    return  hden_array, least_square_array


def find_nH_and_Z_LSF_log(model_Q,  ions_to_use, reference_log_metal = -1.0, true_Q = 18, true_nH = 1e-4):
    #print(model_Q)
    model = tab.Table.read(model_Q)
    obs_ion_col = get_true_model(model_Q, Q= true_Q, nH= true_nH)
    #for i in ions_to_use:
    #    print(np.log10(obs_ion_col[i][0]), i)

    hden_array = np.array(model['hden'])

    metal_array = 10**(np.arange(-2, 1.01, 0.01))
    len_metal = len(metal_array)

    number_of_ions = len(ions_to_use)
    # filter to use specific ions
    least_square_2D = np.zeros((len(model), len_metal))
    for i in range(len_metal):
        metal_scaling_linear = metal_array[i] / 10**reference_log_metal
        least_square_array = np.zeros(len(model))
        for ion in ions_to_use:
            least_square_array += (np.log10(model[ion] * metal_scaling_linear) - np.log10(obs_ion_col[ion])) ** 2

        least_square_2D[:, i] = least_square_array/number_of_ions

    ind = np.unravel_index(np.argmin(least_square_2D, axis=None), least_square_2D.shape)

    #print('nH =', hden_array[ind[0]], 'Z = ', metal_array[ind[1]])
    #print('LS value= ', np.min(least_square_2D))


    return  hden_array[ind[0]], metal_array[ind[1]], np.min(least_square_2D)

ions_to_use= ['C+3', 'N+3', 'Si+3', 'O+5', 'C+2']
true_Q =18
outpath = '/home/vikram/cloudy_run/figures'

outfile = outpath + '/NH14_log_lsf_out.fits'

q=[14, 15, 16, 17, 18, 19, 20]
narray = []
zarray = []
value_ls=[]
for q_num in q:
    model = '/home/vikram/cloudy_run/anshuman/try_Q{}.fits'.format(q_num)
    nH, Z, min_LS= find_nH_and_Z_LSF_log(model, ions_to_use)
    narray.append(np.log10(nH))
    zarray.append(np.log10(Z))
    value_ls.append(min_LS)

outdata = tab.Table([q, narray, zarray, value_ls], names =('Q', 'nH', 'Z', 'LS'))
outdata.write(outfile, overwrite =  True)


