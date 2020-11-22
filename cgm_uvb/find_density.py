import numpy as np
import astropy.table as tab


def get_true_model(model_q, Q= 18, nH = 1-4):
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


def find_nH_and_Z_LSF(model_Q,  *ions_to_use, reference_log_metal = -1.0):
    #print(model_Q)
    model = tab.Table.read(model_Q)
    obs_ion_col = get_true_model(model_Q)
    #for i in ions_to_use:
    #    print(np.log10(obs_ion_col[i][0]), i)

    hden_array = np.array(model['hden'])

    metal_array = 10**(np.arange(-2, 1.01, 0.02))
    len_metal = len(metal_array)

    number_of_ions = len(ions_to_use)
    # filter to use specific ions
    least_square_2D = np.zeros((len(model), len_metal))
    for i in range(len_metal):
        metal_scaling_linear = metal_array[i] / 10**reference_log_metal
        least_square_array = np.zeros(len(model))
        for ion in ions_to_use:
            least_square_array += (model[ion] * metal_scaling_linear - obs_ion_col[ion]) ** 2

        least_square_2D[:, i] = least_square_array/number_of_ions

    ind = np.unravel_index(np.argmin(least_square_2D, axis=None), least_square_2D.shape)

    #print('nH =', hden_array[ind[0]], 'Z = ', metal_array[ind[1]])
    #print('LS value= ', np.min(least_square_2D))


    return  hden_array[ind[0]], metal_array[ind[1]], np.min(least_square_2D)


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

def bootstrap_error(model_Q,  *ions_to_use, reference_log_metal = -1.0, sample_size= 50, seed= 0):
    print(model_Q)
    f_nH, f_Z, array = find_nH_and_Z_LSF(model_Q, *ions_to_use, reference_log_metal = reference_log_metal)

    np.random.seed(seed)
    ions =  np.array(ions_to_use)
    number_of_ions = len(ions_to_use)

    sigma_n = 0
    sigma_z = 0
    for i in range(sample_size):

        c = np.random.choice(number_of_ions, number_of_ions)
        while len(np.unique(c)) == 1:
            c =  np.random.choice(number_of_ions, number_of_ions)
        rand_ions = ions[c]
        nH, Z, array = find_nH_and_Z_LSF(model_Q, *rand_ions, reference_log_metal =reference_log_metal)
        #print(i, c, nH, Z)
        sigma_n += (f_nH - nH)**2
        sigma_z += (f_Z -Z)**2
        #if (i%50)==0:
        #    print(f_nH, np.sqrt(sigma_n / i), ' :density', i)
        #    print(f_Z, np.sqrt(sigma_z / i), ': metallicity', i)

    e_n = np.sqrt(sigma_n / sample_size)
    e_z = np.sqrt(sigma_z / sample_size)

    print(f_nH, e_n, e_n/f_nH, ':density', f_Z, e_z, e_z/f_Z, ':metal')
    #print( ': metallicity')


    return f_nH, sigma_n/sample_size, f_Z, sigma_z/sample_size



ions_to_use= ['C+3', 'N+3', 'Si+3', 'O+5', 'C+2']
ions= ['C+3', 'C+2', 'Si+3', 'Si+2', 'O+5']

q=[14, 15, 16, 17, 18, 19, 20]
"""

#res = tab.Table(meta = {'stop_NHI' : 15})
for q_num in q:
    model = '/home/vikram/cloudy_run/anshuman_NH15/try_Q{}.fits'.format(q_num)
    nH, Z, min_LS= find_nH_and_Z_LSF_log(model, *ions_to_use)
    print(q_num, np.log10(nH), np.log10(Z), min_LS, np.log10(np.sqrt(min_LS)))
    t= tab.Table([[q_num], [nH], [Z], [min_LS]], names = ('Q', 'nH', 'Z', 'min_LS'))
    #res = tab.vstack((res, t))

#print(res)
#res.write('log_LS_out_16.fits', overwrite= True)
"""
ions_to_use=['C+3', 'N+2', 'Si+3', 'S+2', 'O+5']
ions_to_use=['C+3', 'N+2', 'Si+3',  'O+5']
ions_to_use=['C+3', 'C+2']
np.random.seed(123)

ions_to_use=['C+3', 'N+2']
print(ions_to_use)

true_nH = [1e-5, 1e-4, 1e-3]
for n in true_nH:
    for Q in q:
        #print(Q, n)
        narray = []
        zarray = []
        for q_num in q:
            model = '/home/vikram/cloudy_run/anshuman_NH15/try_Q{}.fits'.format(q_num)
            nH, Z, min_LS = find_nH_and_Z_LSF_log(model, ions_to_use, true_Q=Q, true_nH=n)
            narray.append(nH)
            zarray.append(Z)
            # print(q_num, np.log10(nH), np.log10(Z), min_LS, np.log10(np.sqrt(min_LS)))
            # t = tab.Table([[q_num], [nH], [Z], [min_LS]], names=('Q', 'nH', 'Z', 'min_LS'))i
        logdiffn =np.log10(max(narray)) - np.log10(min(narray))
        diffn =(max(narray)) - (min(narray))

        logdiffz= np.log10(max(zarray)) - np.log10(min(zarray))
        diffz= (max(zarray)) - (min(zarray))
        print('diff n', logdiffn , logdiffz, Q, n)
        #print('diff n', logdiffn / (-np.log10(narray)), logdiffz / (-np.log10(zarray)), Q, n)
        #print('diff nH', logdiffn, diffn,  logdiffn /(-np.log10(n)), diffn/ n)
        #print('diff Z', logdiffz, diffz, logdiffz/ 1, diffz/0.1)




"""

model ='/home/vikram/cloudy_run/anshuman/try_Q14.fits'
ions= ['C+3', 'C+2', 'Si+3', 'Si+2', 'O+5']

#nH, Z, array= find_nH_and_Z_LSF(model, *ions)

model ='/home/vikram/cloudy_run/anshuman/try_Q20.fits'
ions= ['C+3', 'N+3', 'Si+3', 'O+5', 'C+2']

#nH, Z, array= find_nH_and_Z_LSF(model, *ions)

q=[14, 15, 16, 17, 18, 19, 20]

for q_num in q:
    model = '/home/vikram/cloudy_run/anshuman/try_Q{}.fits'.format(q_num)
    nH, Z, min_LS= find_nH_and_Z_LSF(model, *ions)
    print(q_num, nH, Z, min_LS/len(ions), np.log10(np.sqrt(min_LS/len(ions))))

print('for NH = 15')
for q_num in q:
    model = '/home/vikram/cloudy_run/anshuman_NH15/try_Q{}.fits'.format(q_num)
    nH, Z, min_LS= find_nH_and_Z_LSF(model, *ions)
    print(q_num, nH, Z, min_LS/len(ions), np.log10(np.sqrt(min_LS/len(ions))))

print('for NH = 16')
for q_num in q:
    model = '/home/vikram/cloudy_run/anshuman_NH16/try_Q{}.fits'.format(q_num)
    nH, Z, min_LS= find_nH_and_Z_LSF(model, *ions)
    print(q_num, nH, Z, min_LS/len(ions), np.log10(np.sqrt(min_LS/len(ions))))
    
s = 200
model ='/home/vikram/cloudy_run/anshuman/try_Q20.fits'
solution  = bootstrap_error(model, *ions, sample_size= s)
model ='/home/vikram/cloudy_run/anshuman/try_Q19.fits'
solution  = bootstrap_error(model, *ions, sample_size= s)
model ='/home/vikram/cloudy_run/anshuman/try_Q18.fits'
solution  = bootstrap_error(model, *ions, sample_size= s)
model ='/home/vikram/cloudy_run/anshuman/try_Q17.fits'
solution  = bootstrap_error(model, *ions, sample_size= s)
model ='/home/vikram/cloudy_run/anshuman/try_Q16.fits'
solution  = bootstrap_error(model, *ions, sample_size= s)
model ='/home/vikram/cloudy_run/anshuman/try_Q15.fits'
solution  = bootstrap_error(model, *ions, sample_size= s)
model ='/home/vikram/cloudy_run/anshuman/try_Q14.fits'
solution  = bootstrap_error(model, *ions, sample_size= s)

"""


#solution  = bootstrap_error(model, *ions, sample_size=80)



"""
    ions = ["H", "H+",
            "He", "He+", "He+2",
            "C", "C+", "C+2", "C+3", "C+4", "C+5",
            "N", "N+", "N+2", "N+3", "N+4",
            "O", "O+", "O+2", "O+3", "O+4", "O+5", "O+6", "O+7",
            "S", "S+", "S+2", "S+3", "S+4", "S+5",
            "Si", "Si+", "Si+2", "Si+3", "Si+4"]
"""




