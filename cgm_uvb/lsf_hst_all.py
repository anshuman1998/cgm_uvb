#------------------
# for reducing  numpy threads
import os
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
#-----------------
import itertools
import numpy as np
import astropy.table as tab
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import emcee
import corner


#----data
def get_true_model(model_Q, Q= 17):
    """
    :param model: The data where Q18 model is stored
    :return: a row of ion column densities at n_H = 1e-4 cm^-2
    """
    #model = model_Q.split('_Q')[0] + '_Q{}.fits'.format(Q)
    data = tab.Table.read(model_Q)
    true_ion_col = data [data['hden'] == 1e-4]
   # print(true_ion_col)

    return true_ion_col

#----model interpolation
def get_interp_func(model_Qpath, ions_to_use):
    number_of_ions = len(ions_to_use)

    model = tab.Table.read(model_Qpath)
    sorted_model = model[ions_to_use]
    hden_array = np.array(model['hden'])

    model_touple = ()
    for j in range(number_of_ions):
        model_touple += (sorted_model[ions_to_use[j]],)

    # interpolating in log log scale
    logf = interp1d(np.log10(hden_array), np.log10(model_touple), fill_value='extrapolate')

    return logf


def get_list_of_qualified_ions(true_Q = 17, threshold_col = 11.5,
                               logZ = None, logT= None, true_nH = 1e-4 ):

    ions = ["C+", "C+2", "C+3",
            "N+2", "N+3", "N+4",
            "O+", "O+2", "O+3", "O+4", "O+5",
            "S+",
            "S+2", "S+3", "S+4", "S+5",
            "Si+",
            "Si+2", "Si+3",
            "Mg+",
            "Ne+7",
            "Fe+"
            ]

    #column_densities = get_true_model(uvb=true_uvb, Q=true_Q, model_path=model_path,true_nH=true_nH, logT=logT, logZ=logZ)

    # get the Nion for true UVB for the ions to use
    model_true = '/home/vikram/cloudy_run/hst/NH15/Q17/Q{}.fits'.format(true_Q)
    column_densities = get_true_model(model_Q = model_true, Q=true_Q)

    list_of_qualified_ions = []
    for i in ions:
        if np.log10(column_densities[i]) > threshold_col:
            list_of_qualified_ions.append(i)

    print('Number of qualified ions', len(list_of_qualified_ions), ' for log N >', threshold_col )

    return list_of_qualified_ions



def get_a_set_of_ions_for_inference(list_of_qualified_ions, num_ions =8, total_unique_comb = 165, seed = 123):

    # total combinations
    all_unique_combinations = list(itertools.combinations(list_of_qualified_ions, num_ions))

    # choose randomly if
    if len(all_unique_combinations) >= total_unique_comb:
        random.seed(seed)
        final_list = random.sample(all_unique_combinations, total_unique_comb)
    else:
        print('Warning: number of unique combination {} are less '
              'that we want {}'.format(len(all_unique_combinations), total_unique_comb))
        final_list = all_unique_combinations

    print('total number of combinations are', len(all_unique_combinations), 'for num of ions  =', num_ions)

    return final_list




def get_LSF_phot(ions_to_use, model_path, Q_uvb, model_uvb, true_Q, true_uvb, interpolated_grid_3D, true_nH_array, true_logZ_array):
    # run_mcmc(model_Q= model, ions_to_use= ions)
    # ------------------ here is a way to run code
    number_of_ions = len(ions_to_use)

    #------------default array: KEEP AS IT IS IN THE  inference_for_photoionized_cloud sub
    # interpolating grid
    lognH_array = np.arange(-6, -1.999, 0.01)
    logZ_array = np.arange(-3, 1.0001, 0.01)
    number_nH = len(lognH_array)
    number_Z = len(logZ_array)
    #-----------------------------

    # define arrays
    make_nH_true_array = []
    make_Z_true_array = []
    inferred_nH_array = []
    inferred_Z_array= []


    true_Q = 17
    model_true = '/home/vikram/cloudy_run/hst/NH15/Q17/Q{}.fits'.format(true_Q)
    data_col_all = get_true_model(model_Q = model_true, Q=true_Q)

    data_col = []
    for name in ions_to_use:
        try:
            data_col.append(data_col_all[name][0])
        except:
            print(name, data_col_all)

    obs_ion_col = np.log10(np.array(data_col))
    # print(obs_ion_col)

    # set least square zero array
    least_square_2D = np.zeros((number_nH, number_Z))
    # calculate least squares
    for k in range(number_of_ions):
        least_square_2D += (interpolated_grid_3D[:, :, k] - obs_ion_col[k]) ** 2

    # getting index for th min of least square
    ind = np.unravel_index(np.argmin(least_square_2D, axis=None), least_square_2D.shape)

    inferred_nH = lognH_array[ind[0]]
    inferred_Z = logZ_array[ind[1]]

    # now write the files

    inferred_nH_array.append(inferred_nH)
    inferred_Z_array.append(inferred_Z)
    make_nH_true_array.append(true_nH)
    make_Z_true_array.append(true_logZ)
    print(':->', model_uvb, Q_uvb, number_of_ions, true_nH, true_logZ, 'num_ions:', number_of_ions)



    """

    for true_nH in true_nH_array:
        for true_logZ in true_logZ_array:
            # getting lsf array
            lognH_true = np.log10(true_nH)


            # get true data
            data_col_all = get_true_model(uvb= true_uvb, model_path=model_path, Q=true_Q, true_nH=true_nH, logZ =true_logZ)
            # converting astropy table row to a list
            data_col = []
            for name in ions_to_use:
                try:
                    data_col.append(data_col_all[name][0])
                except:
                    print(name, data_col_all)

            obs_ion_col = np.log10(np.array(data_col))
            #print(obs_ion_col)

            # set least square zero array
            least_square_2D = np.zeros((number_nH, number_Z))
            # calculate least squares
            for k in range(number_of_ions):
                least_square_2D += (interpolated_grid_3D[:, :, k]  - obs_ion_col[k])**2

            # getting index for th min of least square
            ind = np.unravel_index(np.argmin(least_square_2D, axis=None), least_square_2D.shape)

            inferred_nH = lognH_array[ind[0]]
            inferred_Z = logZ_array[ind[1]]

            # now write the files

            inferred_nH_array.append(inferred_nH)
            inferred_Z_array.append(inferred_Z)
            make_nH_true_array.append(true_nH)
            make_Z_true_array.append(true_logZ)
            print(':->', model_uvb, Q_uvb, number_of_ions, true_nH, true_logZ, 'num_ions:', number_of_ions)


    # create remaining array for table
    size_of_array = len(inferred_nH_array)
    model_uvb_array = [model_uvb] *size_of_array
    moldel_Q_array = [Q_uvb] *size_of_array
    make_true_uvb_array  = [true_uvb] *size_of_array
    make_true_Q_array = [true_Q] * size_of_array
    number_of_ions_array = [number_of_ions] *size_of_array

    # for ion names
    ion_name_string = ions_to_use[0]
    for i in range(number_of_ions - 1):
        ion_name_string = ion_name_string + '_' + ions_to_use[i + 1]

    ion_name_array = [ion_name_string]* size_of_array

    #print(ion_name_array)
    #xx = (make_true_uvb_array, make_true_Q_array, make_nH_true_array, make_Z_true_array, model_uvb_array,
    #moldel_Q_array, inferred_nH_array, inferred_Z_array, number_of_ions_array, ion_name_array)
    #for j in xx:
    #    print(j)


    res_tab = tab.Table([make_true_uvb_array, make_true_Q_array, make_nH_true_array, make_Z_true_array, model_uvb_array,
                         moldel_Q_array, inferred_nH_array, inferred_Z_array, number_of_ions_array, ion_name_array],
                        names = ('true_uvb', 'true_Q', 'true_nH', 'true_Z', 'model_uvb', 'model_Q', 'nH', 'Z', 'n_ions', 'ions'))
    
    """

    return #res_tab



def inference_for_photoionized_cloud():

    """
    # get the ion list from a default true model
    default_true_uvb = 'KS18'
    default_true_Q = 17
    default_true_nH = 1e-4
    default_true_logZ  = -1
    qualified_ion_list= get_list_of_qualified_ions(true_Q= default_true_Q,
                                                   true_nH= default_true_nH, logZ= default_true_logZ)

    """

    # eight selected ions
    qualified_ion_list = ['C+', 'C+3', 'N+2', 'N+4', 'O+5', 'Si+', 'Si+2', 'Si+3']
    reference_log_metal = -1
    n = 3
    true_Q = 17
    true_nH = -4
    true_logZ = -4


    print(qualified_ion_list)


    # interpolating grid
    lognH_array = np.arange(-5, -1.999, 0.02) # grid in log nH
    logZ_array = np.arange(-3, 1.0001, 0.02)
    number_nH = len(lognH_array)
    number_Z = len(logZ_array)
    number_of_qualified_ions = len(qualified_ion_list)
    # setting grid to a zero value everywhere
    grid_3D = np.zeros((number_nH, number_Z, number_of_qualified_ions))

    mod_Q = [14, 15, 16, 17, 18, 19, 20]

    for model_Q in mod_Q:
        print('model Q  ----> ', model_Q)

        # get the full interpolation function list => time consuming part first
        print('performing interpolation')
        model_path = '/home/vikram/cloudy_run/hst/NH15/Q{}/Q{}.fits'.format(model_Q, model_Q)
        full_interp_function_list = get_interp_func(model_Qpath=model_path, ions_to_use=qualified_ion_list)
        print('interpolation done')

        # get interpolated 2D array

        for i in range(number_nH):
            val_lognH = lognH_array[i]
            col_ref_Z = full_interp_function_list(val_lognH)
            for j in range(number_Z):
                val_logZ = logZ_array[j]
                metal_scaling_linear = 10 ** val_logZ / 10 ** reference_log_metal
                grid_3D[i, j, :] = np.log10(10 ** col_ref_Z * metal_scaling_linear)

        print('3D grid is done')

        filename_string = 'photo_model_KS18_Q{}_lsf_inference_hst'.format(model_Q)
        #output_file_name = outpath + '/' + filename_string + '.fits'
        result_table = tab.Table()

        ion_set = get_a_set_of_ions_for_inference(list_of_qualified_ions=qualified_ion_list, num_ions=n)

        sigma_noise = 0.001
        noise = np.random.normal(0, sigma_noise, len(ion_set[0]))

        for ions_to_use in ion_set:
            print('==== runnig for a set of ion', ions_to_use)

            # get the 3D grids functions for the ions to use
            trimed_3D_grid = np.zeros((number_nH, number_Z, n))
            for m in range(n):
                ion_index_in_qualified_array = qualified_ion_list.index(ions_to_use[m])
                trimed_3D_grid[:, :, m] = grid_3D[:, :, ion_index_in_qualified_array]

            # lsf
            model_true = '/home/vikram/cloudy_run/hst/NH15/Q17/Q17.fits'
            data_col_all = get_true_model(model_Q=model_true)

            data_col = []
            for name in ions_to_use:
                try:
                    data_col.append(data_col_all[name][0])
                except:
                    print(name, data_col_all)

            data_col = 10 ** (np.log10(data_col) + noise)

            obs_ion_col = np.log10(np.array(data_col))
            # print(obs_ion_col)

            # set least square zero array
            least_square_2D = np.zeros((number_nH, number_Z))
            # calculate least squares
            for k in range(len(ions_to_use)):
                least_square_2D += (trimed_3D_grid[:, :, k] - obs_ion_col[k]) ** 2

            # getting index for th min of least square
            ind = np.unravel_index(np.argmin(least_square_2D, axis=None), least_square_2D.shape)

            inferred_nH = lognH_array[ind[0]]
            inferred_Z = logZ_array[ind[1]]

            # now write the files

            #inferred_nH_array.append(inferred_nH)
            #inferred_Z_array.append(inferred_Z)
            #make_nH_true_array.append(true_nH)
            #make_Z_true_array.append(true_logZ)
            print(':->', model_Q, 'nH:', true_nH, inferred_nH, 'Z:', true_logZ, inferred_Z, 'num_ions:', len(ions_to_use))















#    return

inference_for_photoionized_cloud()


"""
    for n in number_of_ions_array:
        #fits_filename = 'photo_n_ions_{}.fits'.format(n)
       print('----------->for n number of ions', n, 'uvb_model', model_uvb, model_Q)

        ion_set = get_a_set_of_ions_for_inference(list_of_qualified_ions= qualified_ion_list, num_ions= n, total_unique_comb= total_ion_comb)

        for ions_to_use in ion_set:

            print('==== runnig for a set of ion')

            # get the 3D grids functions for the ions to use
            trimed_3D_grid = np.zeros((number_nH, number_Z, n))
            for m in range(n):
                ion_index_in_qualified_array  = qualified_ion_list.index(ions_to_use[m])
                trimed_3D_grid[:, :, m]  = grid_3D[:, :, ion_index_in_qualified_array]


            if true_uvb_model == 'all':
                true_Q = [14, 15, 16, 17, 18, 19, 20]
                uvb_models = ['KS18', 'FG20', 'P19', 'HM12']

                for true_uvb in uvb_models:
                    if true_uvb == 'KS18':
                        for Q in true_Q:
                            res = get_LSF_phot(ions_to_use= ions_to_use, model_path= model_path, Q_uvb=model_Q, model_uvb= model_uvb,
                                               true_uvb= true_uvb, true_Q=Q, interpolated_grid_3D= trimed_3D_grid,
                                               true_nH_array= true_nH_array, true_logZ_array= true_logZ_array)

                            result_table = tab.vstack([result_table, res])


                    else:
                        Q  =18
                        res = get_LSF_phot(ions_to_use=ions_to_use, model_path=model_path, Q_uvb=model_Q, model_uvb= model_uvb,
                                           true_uvb=true_uvb, true_Q=Q, interpolated_grid_3D=trimed_3D_grid,
                                           true_nH_array= true_nH_array, true_logZ_array= true_logZ_array)

                        result_table = tab.vstack([result_table, res])


            #get_LSF_phot(ions_to_use, model_path, Q_uvb, model_uvb, true_Q, true_uvb, interp_logf):

            else:
                'code runs for all uvb models by default'


    result_table.write(output_file_name, overwrite = True)







def get_LSF_phot(ions_to_use, model_path, Q_uvb, true_Q = 17): #, interpolated_grid_3D, true_nH_array, true_logZ_array):

    #
    number_of_ions = len(ions_to_use)
    
    # get the Nion for true UVB for the ions to use
    model_true = '/home/vikram/cloudy_run/hst/NH15/Q17/Q{}.fits'.format(true_Q)
    data_col_all = get_true_model(model_Q = model_true, Q=true_Q)
    # converting astropy table row to a list
    data_col = []
    for name in ions_to_use:
        data_col.append(data_col_all[name][0])

    obs_ion_col = np.log10(np.array(data_col))


    # interpolating grid
    lognH_array = np.arange(-6, -1.999, 0.01)
    logZ_array = np.arange(-3, 1.0001, 0.01)
    number_nH = len(lognH_array)
    number_Z = len(logZ_array)
    
    # setting grid to a zero value everywhere
    grid_3D = np.zeros((number_nH, number_Z, number_of_qualified_ions))

    # get interpolated 2D array

    for i in range(number_nH):
        val_lognH = lognH_array[i]
        # print('at n_H array ind: ', i)
        for j in range(number_Z):
            val_logZ = logZ_array[j]
            for k in range(number_of_qualified_ions):
                # print('==>', i, lognH, logT)
                # print(interp_logf[i](lognH, logT), i, lognH, logT)
                col_mod = full_interp_function_list[k](val_lognH, val_logZ)[0]
                grid_3D[i, j, k] = col_mod

    print('3D grid is done')

    
    
    # set least square zero array
    least_square = np.zeros((number_nH, number_Z))
    # calculate least squares
    for k in range(number_of_ions):
        least_square += (interpolated_grid_3D[:, :, k] - obs_ion_col[k]) ** 2
        
        
    

    
    model_Q = give_a_file_location
    

    #------------default array: KEEP AS IT IS IN THE  inference_for_photoionized_cloud sub
    # interpolating grid
    lognH_array = np.arange(-6, -1.999, 0.01)
    logZ_array = np.arange(-3, 1.0001, 0.01)
    number_nH = len(lognH_array)
    number_Z = len(logZ_array)
    #-----------------------------

    # define arrays
    make_nH_true_array = []
    make_Z_true_array = []
    inferred_nH_array = []
    inferred_Z_array= []

    for true_nH in true_nH_array:
        for true_logZ in true_logZ_array:
            # getting lsf array
            lognH_true = np.log10(true_nH)


            # get true data
            data_col_all = get_true_model(uvb= true_uvb, model_path=model_path, Q=true_Q, true_nH=true_nH, logZ =true_logZ)
            # converting astropy table row to a list
            data_col = []
            for name in ions_to_use:
                try:
                    data_col.append(data_col_all[name][0])
                except:
                    print(name, data_col_all)

            obs_ion_col = np.log10(np.array(data_col))
            #print(obs_ion_col)

            # set least square zero array
            least_square_2D = np.zeros((number_nH, number_Z))
            # calculate least squares
            for k in range(number_of_ions):
                least_square_2D += (interpolated_grid_3D[:, :, k]  - obs_ion_col[k])**2

            # getting index for th min of least square
            ind = np.unravel_index(np.argmin(least_square_2D, axis=None), least_square_2D.shape)

            inferred_nH = lognH_array[ind[0]]
            inferred_Z = logZ_array[ind[1]]

            # now write the files

            inferred_nH_array.append(inferred_nH)
            inferred_Z_array.append(inferred_Z)
            make_nH_true_array.append(true_nH)
            make_Z_true_array.append(true_logZ)
            print(':->', model_uvb, Q_uvb, number_of_ions, true_nH, true_logZ, 'num_ions:', number_of_ions)


    # create remaining array for table
    size_of_array = len(inferred_nH_array)
    model_uvb_array = [model_uvb] *size_of_array
    moldel_Q_array = [Q_uvb] *size_of_array
    make_true_uvb_array  = [true_uvb] *size_of_array
    make_true_Q_array = [true_Q] * size_of_array
    number_of_ions_array = [number_of_ions] *size_of_array

    # for ion names
    ion_name_string = ions_to_use[0]
    for i in range(number_of_ions - 1):
        ion_name_string = ion_name_string + '_' + ions_to_use[i + 1]

    ion_name_array = [ion_name_string]* size_of_array

    #print(ion_name_array)
    #xx = (make_true_uvb_array, make_true_Q_array, make_nH_true_array, make_Z_true_array, model_uvb_array,
    #moldel_Q_array, inferred_nH_array, inferred_Z_array, number_of_ions_array, ion_name_array)
    #for j in xx:
    #    print(j)


    res_tab = tab.Table([make_true_uvb_array, make_true_Q_array, make_nH_true_array, make_Z_true_array, model_uvb_array,
                         moldel_Q_array, inferred_nH_array, inferred_Z_array, number_of_ions_array, ion_name_array],
                        names = ('true_uvb', 'true_Q', 'true_nH', 'true_Z', 'model_uvb', 'model_Q', 'nH', 'Z', 'n_ions', 'ions'))

    return res_tab




#----for mcmc
def log_likelihood(theta, interp_logf, obs_ion_col, col_err, reference_log_metal = -1.0):
    
    #For a gaussian distributed errors
    #:param theta: parameters [nH, Z]
    #:param x: data x
    #:param y: data y
    #:param yerr: data err
    #:return:
    lognH, logZ =  theta
    # get metal ion column density for n_H and Z = 0.1
    col = 10 ** interp_logf(lognH)
    # scale the column densities by the metallicity Z
    metal_scaling_linear = 10 ** logZ / 10 ** reference_log_metal
    model_col = np.log10(col * metal_scaling_linear)

    lnL = -0.5 * np.sum(np.log(2 * np.pi * col_err ** 2) + (obs_ion_col - model_col) ** 2 / col_err ** 2)

    return lnL

def log_prior(theta):
    lognH, logZ =  theta
    # flat prior
    if -6 < lognH < -2 and -2 < logZ < 1 :
        return 0.0
    return -np.inf

def log_posterior(theta, interp_func, data_col, sigma_col):
    log_p = log_prior(theta) + \
            log_likelihood(theta, interp_logf = interp_func, obs_ion_col = data_col, col_err = sigma_col)

    return log_p


def run_mcmc(model_Q, ions_to_use, true_Q =18, figname = 'test.pdf', same_error = True, sigma_noise = 0.1):
    # run_mcmc(model_Q= model, ions_to_use= ions)
    # ------------------ here is a way to run code
    truths = [-4, -1]  # (lognH, logZ) true values
    number_of_ions = len(ions_to_use)


    #np.random.seed(0)
    if same_error:
        sigma_col = sigma_noise * np.ones(number_of_ions)
        noise  = np.random.normal(0, sigma_noise, number_of_ions)
        print(noise)
    else:
        sigma_col = np.random.uniform(0.1, 0.3, number_of_ions)


    data_col_all = get_true_model(model_Q, Q=true_Q)

    # converting astropy table row to a list
    data_col = []
    for name in ions_to_use:
        data_col.append(data_col_all[name][0])

    print(data_col)
    # shifting the column values
    data_col =  10**(np.log10(data_col) + noise)
    print(data_col)

    print(np.log10(data_col), sigma_col)

    interp_logf = get_interp_func(model_Q, ions_to_use)

    # Here we'll set up the computation. emcee combines multiple "walkers",
    # each of which is its own MCMC chain. The number of trace results will
    # be nwalkers * nsteps

    ndim = 2  # number of parameters in the model
    nwalkers = 50  # number of MCMC walkers
    nsteps = 5000  # number of MCMC steps to take

    # set theta near the maximum likelihood, with
    n_guess = np.random.uniform(-5, -2, nwalkers)
    z_guess = np.random.uniform(-2, 1, nwalkers)
    starting_guesses = np.vstack((n_guess, z_guess)).T  # initialise at a tiny sphere

    # Here's the function call where all the work happens:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(interp_logf, np.log10(data_col), sigma_col))
    sampler.run_mcmc(starting_guesses, nsteps, progress=True)

    # find out number of steps
    tau = sampler.get_autocorr_time()  # number of steps needed to forget the starting position
    #print(tau)
    thin = int(np.mean(tau) / 2)  # use this number for flattning the sample as done below
    #thin = 100
    flat_samples = sampler.get_chain(discard=thin * 20, thin= 5, flat=True)
    # we are discarding some initial steps roughly 5 times the autocorr_time steps
    # then we thin by about half the autocorrelation time steps for plotting => one does not have to do this step

    labels = ['log nH', 'log Z']
    uvb_q= 17 #int((model_Q.split('try_Q')[-1]).split('.fits')[0])

    if uvb_q == true_Q:
        fig = corner.corner(flat_samples, labels=labels, truths=truths, quantiles=[0.16, 0.5, 0.84],
            show_titles=True, title_kwargs={"fontsize": 12})
    else:
        fig = corner.corner(flat_samples, labels=labels, quantiles=[0.16, 0.5, 0.84],
            show_titles=True, title_kwargs={"fontsize": 12})

    fig.savefig(figname)

    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print(labels[i], '=', mcmc[1], q[0], q[1])


    return flat_samples, ndim


ions_to_use= ['C+3', 'N+3', 'Si+3', 'O+5', 'C+2']
true_Q =17
outpath = '/home/vikram/cloudy_run/hst'

outfile = outpath + '/NH15.fits'



uvb_array= [17]
out_tab =  tab.Table()
for uvb_q in uvb_array:
    model_Q = '/home/vikram/cloudy_run/hst/NH15/Q17/Q{}.fits'.format(uvb_q)
    name = 'Q17'
    figname = outpath + '/' + name + '.pdf'

    flat_samples, ndim = run_mcmc(model_Q=model_Q, ions_to_use=ions_to_use, true_Q=true_Q, figname=figname)
    # to efficiently save numpy array
    save_file_name = outpath + '/' + name
    np.save(save_file_name, flat_samples)

    out =[[uvb_q]]
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        out.append([mcmc[1]])
        out.append([q[0]])
        out.append([q[1]])

    print(out)
    t = tab.Table(out, names = ('Q', 'nH', 'n16', 'n84', 'Z', 'Z16', 'Z84'))
    out_tab = tab.vstack((out_tab, t))


out_tab.write(outfile, overwrite = True)

"""