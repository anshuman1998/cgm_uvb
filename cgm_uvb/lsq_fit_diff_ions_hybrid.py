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
import random
from scipy.interpolate import interp2d
import multiprocessing as mp


seed = 123


def get_true_model(uvb, model_path, Q = 18, true_nH=1e-4, logT = None, logZ = None):
    """
    :param uvb:
    :param output_path:
    :param Q:
    :param nH:
    :param logZ:
    :param logT:
    :return:
    """
    # for photoionized case:
    if logZ != None:
        fname = (logZ + 4) * 100
        model_file = model_path + '/try_{}_Q{}_Z{:.0f}.fits'.format(uvb, Q, fname)

    if logT != None:
        tname = logT * 100
        model_file = model_path + '/try_{}_Q{}_logT{:.0f}.fits'.format(uvb, Q, tname)

    data = tab.Table.read(model_file)
    true_ion_col = data[data['hden'] == true_nH]

    return true_ion_col



#----model interpolation
def get_interp_func_for_all_hybrid(model_path, ions_to_use, Q_uvb = 18, uvb = 'KS18', logT = 5.5):
    logZ = np.arange(-3, 1.1, 0.02) # hardcoded
    #get nH array
    model_try = model_path + '/try_{}_Q{}_logT{:.0f}.fits'.format(uvb, Q_uvb, logT*100)
    model = tab.Table.read(model_try)
    lognH = np.log10(np.array(model['hden']))

    reference_log_metal = -1

    interpolation_function_list = []
    for ion in ions_to_use:
        #model[ion][model[ion] == 0] = 1e-15  # for avoiding log10 (0) error
        z = np.zeros((len(lognH), len(logZ)))
        for i in range(len(logZ)):
            metal_scaling_linear = 10**logZ[i]/ 10 ** reference_log_metal
            z [:, i] = np.log10(model[ion]*metal_scaling_linear) #--- for log - log interpolation
        f = interp2d(lognH, logZ, z.T)
        interpolation_function_list.append(f)

    return interpolation_function_list



def get_LSF_hybrid(ions_to_use, model_path, Q_uvb, model_uvb, true_Q, true_uvb, interpolated_grid_3D, true_nH_array,
                   true_logZ_array, logT= 5.5, reference_log_metal = -1):


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

    for true_nH in true_nH_array:
        # getting lsf array
        # get true data---- >------- be careful for hybrid model i.e keep logZ = None
        data_col_all = get_true_model(uvb=true_uvb, model_path=model_path, Q=true_Q, true_nH=true_nH, logT=logT)
        # converting astropy table row to a list
        data_col = []
        for name in ions_to_use:
            try:
                data_col.append(data_col_all[name][0])
            except:
                print(name, data_col_all)

        obs_ion_col_Zref = np.array(data_col)
        # print(obs_ion_col)


        for true_logZ in true_logZ_array:

            metal_scaling_linear = 10**true_logZ/ 10 ** reference_log_metal
            obs_ion_col = np.log10(obs_ion_col_Zref*metal_scaling_linear)

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

def get_list_of_qualified_ions(true_uvb, true_Q = 18, model_path = '', threshold_col = 11,
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

    column_densities = get_true_model(uvb=true_uvb, Q=true_Q, model_path=model_path,
                                      true_nH=true_nH, logT=logT, logZ=logZ)

    list_of_qualified_ions = []
    for i in ions:
        if np.log10(column_densities[i]) > threshold_col:
            list_of_qualified_ions.append(i)

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



def inference_for_hybrid_cloud(model_uvb = 'KS18', model_Q = 18, true_uvb_model = 'all',
                                     true_nH_array= [1e-5, 1e-4, 1e-3],
                                     true_logZ_array = [-2,-1, 0],
                                     number_of_ions_array = [2, 3, 4, 5, 6, 7, 8],
                                     model_path='/home/vikram/cloudy_run/hybrid_NH15',
                                     outpath  = '/home/vikram/tmp/new', logT = 5.5,
                                     total_ion_comb = 165):

    # get the ion list from a default true model
    default_true_uvb = 'KS18'
    default_true_Q = 18
    default_true_nH = 1e-4
    default_true_logZ  = None
    qualified_ion_list= get_list_of_qualified_ions(true_uvb = default_true_uvb, true_Q= default_true_Q,
                                                   true_nH= default_true_nH, logZ= default_true_logZ,
                                                   model_path= model_path, logT = logT)

    print(qualified_ion_list)

    # get the full interpolation function list => time consuming part first
    print('performing interpolation')
    full_interp_function_list = get_interp_func_for_all_hybrid(model_path=model_path, ions_to_use=qualified_ion_list,
                                                              Q_uvb=model_Q, uvb= model_uvb, logT= logT)
    print('interpolation done')

    # interpolating grid
    lognH_array = np.arange(-6, -1.999, 0.01)
    logZ_array = np.arange(-3, 1.0001, 0.01)
    number_nH = len(lognH_array)
    number_Z = len(logZ_array)
    number_of_qualified_ions = len(qualified_ion_list)
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


    filename_string = 'hybrid_model_{}_Q{}_logT{:.0f}_lsf_inference'.format(model_uvb, model_Q, logT*100)
    output_file_name = outpath + '/' + filename_string +'.fits'
    result_table =  tab.Table()

    for n in number_of_ions_array:
        #fits_filename = 'photo_n_ions_{}.fits'.format(n)
        print('----------->for n number of ions', n, 'uvb_model', model_uvb, model_Q)

        ion_set = get_a_set_of_ions_for_inference(list_of_qualified_ions= qualified_ion_list, num_ions= n, total_unique_comb= total_ion_comb)

        for ions_to_use in ion_set:

            print('==== running for a set of ion')

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
                            res = get_LSF_hybrid(ions_to_use= ions_to_use, model_path= model_path, Q_uvb=model_Q, model_uvb= model_uvb,
                                               true_uvb= true_uvb, true_Q=Q, interpolated_grid_3D= trimed_3D_grid,
                                               true_nH_array= true_nH_array, true_logZ_array= true_logZ_array)

                            result_table = tab.vstack([result_table, res])


                    else:
                        Q  =18
                        res = get_LSF_hybrid(ions_to_use=ions_to_use, model_path=model_path, Q_uvb=model_Q, model_uvb= model_uvb,
                                           true_uvb=true_uvb, true_Q=Q, interpolated_grid_3D=trimed_3D_grid,
                                           true_nH_array= true_nH_array, true_logZ_array= true_logZ_array)

                        result_table = tab.vstack([result_table, res])


            #get_LSF_phot(ions_to_use, model_path, Q_uvb, model_uvb, true_Q, true_uvb, interp_logf):

            else:
                'code runs for all uvb models by default'


    result_table.write(output_file_name, overwrite = True)


#----------test
#inference_for_photoionized_cloud(model_path='/home/vikram/cloudy_run/metal_NH15_new', total_ion_comb=10)

def run_parallel(model_uvb, model_Q):
    outpath  = '/home/vikram/cloudy_run/diff_op/hybrid_NH15'
    inference_for_hybrid_cloud(model_uvb= model_uvb, model_Q= model_Q, outpath= outpath)
    return

# runnning in parallel
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


print(uvb_models, '==== models all')
print(the_Q_values, '===Q val')

pool = mp.Pool(processes=5)
results = [pool.apply_async(run_parallel, args=(for_uvb_model, for_Q,)) for for_uvb_model, for_Q
           in zip(uvb_models, the_Q_values)]
output = [p.get() for p in results]

