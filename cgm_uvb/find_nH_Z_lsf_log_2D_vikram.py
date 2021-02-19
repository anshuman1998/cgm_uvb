import numpy as np
import astropy.table as tab
from scipy.interpolate import interp2d


def get_true_model(model_path, Q= 18, true_nH = 1e-4,  logZ = -1, uvb = 'KS18'):
    """
    :param model: The data where Q18 model is stored
    :return: a row of ion column densities at n_H = 1e-4 cm^-2
    """
    model = model_path + '/try_{}_Q{}_Z{:.0f}.fits'.format(uvb, Q, (logZ+4)*100)
    data = tab.Table.read(model)
    true_ion_col = data [data['hden'] == true_nH]
   # print(true_ion_col)

    return true_ion_col

#----model interpolation
def get_interp_func(model_path, ions_to_use, Q_uvb, uvb = 'KS18'):
    logZ = np.around(np.arange(-3, 1, 0.05), decimals = 2) # hardcoded
    #get nH array
    logZ_try = -1
    model_try = model_path + '/try_{}_Q{}_Z{:.0f}.fits'.format(uvb, Q_uvb, (logZ_try+4)*100)
    model = tab.Table.read(model_try)
    lognH = np.log10(np.array(model['hden']))

    interpolation_function_list = []
    for ion in ions_to_use:
        z = np.zeros((len(lognH), len(logZ)))
        for i in range(len(logZ)):
            model = model_path + '/try_{}_Q{}_Z{:.0f}.fits'.format(uvb, Q_uvb, (logZ[i]+4)*100)
            d = tab.Table.read(model)
            d[ion][d[ion] == 0 ] = 1e-15 # for avoiding log10 (0) error
            z [:, i] = np.log10(d[ion]) #--- for log - log interpolation
        f = interp2d(lognH, logZ, z.T)
        interpolation_function_list.append(f)

    return interpolation_function_list




def get_LSF(model_path, Q_uvb, model_uvb, ions_to_use, true_Q, true_uvb):
    # run_mcmc(model_Q= model, ions_to_use= ions)
    # ------------------ here is a way to run code
    number_of_ions = len(ions_to_use)

    # get interpolation functions for model
    interp_logf = get_interp_func(model_path = model_path, ions_to_use = ions_to_use, Q_uvb = Q_uvb, uvb = model_uvb)

    print('interpolation done')


    # ---- the nine combinations
    true_nH_array = [1e-4]
    true_logZ_array  = [-1]
    for true_nH in true_nH_array:
        for true_logZ in true_logZ_array:
            # getting lsf array
            lognH_true = np.log10(true_nH)
            lognH_array = np.arange(lognH_true - 1, lognH_true + 1, 0.01)
            logZ_array = np.arange(true_logZ - 1, true_logZ + 1, 0.01)
            number_nH = len(lognH_array)
            number_Z = len(logZ_array)

            least_square_2D = np.zeros((number_nH, number_Z))
            # get true data
            data_col_all = get_true_model(model_path, Q=true_Q, true_nH=true_nH, logZ =true_logZ, uvb= true_uvb)
            # converting astropy table row to a list
            data_col = []
            for name in ions_to_use:
                data_col.append(data_col_all[name][0])

            obs_ion_col = np.log10(np.array(data_col))
            #print(obs_ion_col)

            for i in range(number_nH):
                val_lognH = lognH_array[i]
                for j in range(number_Z):
                    val_logZ = logZ_array[j]
                    col = []
                    for k in range(number_of_ions):
                        # print('==>', i, lognH, logT)
                        # print(interp_logf[i](lognH, logT), i, lognH, logT)
                        col_mod = interp_logf[k](val_lognH, val_logZ)[0]
                        col.append(col_mod)

                    model_col = np.array(col)
                    least_square_2D[i, j] = np.sum((model_col - obs_ion_col) ** 2)

            ind = np.unravel_index(np.argmin(least_square_2D, axis=None), least_square_2D.shape)
            solution_n = lognH_array[ind[0]]
            solution_Z = logZ_array[ind[1]]
            min_val_lsf = np.min(least_square_2D)


            print(lognH_array[ind[0]], logZ_array[ind[1]], np.min(least_square_2D), 'for (nH, Z)',
                np.log10(true_nH), true_logZ, 'true uvb: Q', true_Q, true_uvb, 'model uvb: Q', Q_uvb, model_uvb )

    return solution_n, solution_Z, min_val_lsf



#ions_to_use= ['C+3', 'N+3', 'Si+3', 'O+5', 'C+2']

#ions_to_use= ['C+3', 'N+3', 'N+4', 'O+2', 'O+5']

ions_to_use= ['O+5', 'N+4', 'C+3', 'Si+3', 'Ne+7']

model_path  = '/home/vikram/cloudy_run/rescaled_metal_NH15'
outpath = '/home/vikram/cloudy_run/figures/rescaled'
outfile = outpath + '/rescaled_NH15_log_lsf_out.fits'

#model_path  = '/home/vikram/cloudy_run/metal_NH15_new'
#outpath = '/home/vikram/cloudy_run/figures/2D'
#outfile = outpath + '/NH15_log_lsf_out.fits'


uvb_array = ['KS18', 'KS18', 'KS18', 'KS18', 'KS18', 'KS18', 'KS18', 'P19', 'FG20', 'HM12']
Q_array= [14, 15, 16, 17, 18, 19, 20, 18, 18, 18]

narray = []
zarray = []
value_ls=[]
for uvb, q in zip(uvb_array, Q_array):
    nH, Z, min_LS= get_LSF(model_path = model_path, Q_uvb = q, model_uvb= uvb,
        ions_to_use = ions_to_use, true_Q =18, true_uvb = 'KS18')
    narray.append(nH)
    zarray.append(Z)
    value_ls.append(min_LS)

uvb_column = ['Q14', 'Q15', 'Q16', 'Q17', 'Q18', 'Q19', 'Q20', 'P19', 'FG20', 'HM12']
outdata = tab.Table([uvb_column, Q_array, narray, zarray, value_ls], names =('uvb', 'Q', 'nH', 'Z', 'LS'))

outdata.write(outfile, overwrite =  True)

