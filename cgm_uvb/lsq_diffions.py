import itertools
import numpy as np
import astropy.table as tab

def get_true_model(Q= 18, nH = 1e-4):
    """
    :param model: The data where Q18 model is stored
    :return: a row of ion column densities at n_H = 1e-4 cm^-2
    """
    ks = [14,15,16,17,18,19,20]
    if Q=='FG20':
        model = '/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/output/try_FG20.fits'
        #model = '/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/hybrid/try_FG20_Q18_logT550.fits'
        data = tab.Table.read(model)
        true_ion_col = data [data['hden'] == nH]
    
    elif Q=='P19':
        model = '/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/output/try_P19.fits'
        #model = '/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/hybrid/try_P19_Q18_logT550.fits'
        data = tab.Table.read(model)
        true_ion_col = data [data['hden'] == nH]

    elif Q in ks:
        model = '/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/output/try_Q{}.fits'.format(Q)
        #model = '/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/hybrid/try_KS18_Q{}_logT550.fits'.format(Q)
        data = tab.Table.read(model)
        true_ion_col = data [data['hden'] == nH]
    return true_ion_col

def find_nH_and_Z_LSF_log(model_Q,  ions_to_use, reference_log_metal = -1.0, true_Q = 18, true_nH = 1e-4):
    #print(model_Q)
    model = tab.Table.read(model_Q)
    obs_ion_col = get_true_model(Q= true_Q, nH= true_nH)
    hden_array = np.array(model['hden'])

    metal_array = 10**(np.arange(-2, 1.01, 0.01))
    len_metal = len(metal_array)
    number_of_ions = len(ions_to_use)
    least_square_2D = np.zeros((len(model), len_metal))
    
    for i in range(len_metal):
        metal_scaling_linear = metal_array[i] / 10**reference_log_metal
        least_square_array = np.zeros(len(model))
        for ion in ions_to_use:
            least_square_array += (np.log10(model[ion] * metal_scaling_linear) - np.log10(obs_ion_col[ion])) ** 2

        least_square_2D[:, i] = least_square_array/number_of_ions

    ind = np.unravel_index(np.argmin(least_square_2D, axis=None), least_square_2D.shape)

    return  hden_array[ind[0]], metal_array[ind[1]], np.min(least_square_2D)

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

def combination(true_Q=18,true_nH=1e-4,no=8,setno=0,ions_wecanuse=ions):
    ion_combis=list(set(itertools.combinations(ions_wecanuse, no)))
    combno = len(ion_combis)
    if len(ion_combis)<setno: #To ensure that the number of combinations requested is always less than possible no. of combinations
        print("Warning, invalid call!",len(ion_combis))
        return 0
    return ion_combis[setno],combno

#Getting the ions that we can use and number of iterations
obs_ion_col = get_true_model(Q=18,nH=1e-4)
ions_wecanuse =[]                      #This is the number of ions that we can use, assuming true obs as Q18, nH=1e-4
for i in ions:
    if np.log10(obs_ion_col[i][0])>11.0:
        ions_wecanuse.append(i)

one_combo,combno = combination(true_Q=18,true_nH=1e-4,no=8,setno=0,ions_wecanuse=ions_wecanuse)
print(combno)
if combno>165:
    combno=165 #Ensuring number of draws is 165 at most
q=[14, 15, 16, 17, 18, 19, 20,'FG20','P19']
#oneq = [18]
ks = [14,15,16,17,18,19,20]

true_nH = [1e-5, 1e-4, 1e-3]
#true_nH=[1e-5]
mdn = []
mdz = []
uvb = []
tnH = []
ion = []

nq14=[]
zq14=[]
nq15=[]
zq15=[]
nq16=[]
zq16=[]
nq17=[]
zq17=[]
nq18=[]
zq18=[]
nq19=[]
zq19=[]
nq20=[]
zq20=[]
nqfg=[]
zqfg=[]
nqp=[]
zqp=[]
h = [nq14,nq15,nq16,nq17,nq18,nq19,nq20,nqfg,nqp]
z = [zq14,zq15,zq16,zq17,zq18,zq19,zq20,zqfg,zqp]

print("Start")
for n in true_nH:
    for Q in q:
        print(Q, n)
        logdiffn = []
        logdiffz = []
        
        for i in range(combno):
            narray = []
            zarray = []
            #Choosing ions for the model from the list made above
            no = 8
            ions_to_use,ionno=combination(true_Q=Q,true_nH=n,no=no,setno=0,ions_wecanuse=ions_wecanuse)


            for q_num,nfile,zfile in zip(q,h,z):
                if q_num=='FG20':
                    model = '/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/output/try_FG20.fits'
                    #model = '/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/hybrid/try_FG20_Q18_logT550.fits'
                    if Q in ks:
                        nH, Z, min_LS = find_nH_and_Z_LSF_log(model, ions_to_use, true_Q=Q, true_nH=n)
                    if Q=='FG20':
                        nH, Z, min_LS = find_nH_and_Z_LSF_log(model, ions_to_use, true_Q=Q, true_nH=n)
                    if Q=='P19':
                        nH, Z, min_LS = find_nH_and_Z_LSF_log(model, ions_to_use, true_Q=Q, true_nH=n)
                    narray.append(nH)
                    zarray.append(Z)
                    nfile.append(str(round(np.log10(narray[7]),4)))
                    zfile.append(str(round(np.log10(zarray[7]),4)))
                    
                elif q_num=='P19':
                    model = '/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/output/try_P19.fits'
                    #model = '/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/hybrid/try_P19_Q18_logT550.fits'
                    if Q in ks:
                        nH, Z, min_LS = find_nH_and_Z_LSF_log(model, ions_to_use, true_Q=Q, true_nH=n)
                    if Q=='FG20':
                        nH, Z, min_LS = find_nH_and_Z_LSF_log(model, ions_to_use, true_Q=Q, true_nH=n)
                    if Q=='P19':
                        nH, Z, min_LS = find_nH_and_Z_LSF_log(model, ions_to_use, true_Q=Q, true_nH=n)
                    narray.append(nH)
                    zarray.append(Z)
                    nfile.append(str(round(np.log10(narray[8]),4)))
                    zfile.append(str(round(np.log10(zarray[8]),4)))
                    
                elif q_num in ks:
                    model = '/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/output/try_Q{}.fits'.format(q_num)
                    #model = '/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/hybrid/try_KS18_Q{}_logT550.fits'.format(q_num)
                    if Q in ks:
                        nH, Z, min_LS = find_nH_and_Z_LSF_log(model, ions_to_use, true_Q=Q, true_nH=n)
                    if Q=='FG20':
                        nH, Z, min_LS = find_nH_and_Z_LSF_log(model, ions_to_use, true_Q=Q, true_nH=n)
                    if Q=='P19':
                        nH, Z, min_LS = find_nH_and_Z_LSF_log(model, ions_to_use, true_Q=Q, true_nH=n)
                    narray.append(nH)
                    zarray.append(Z)
                    
                    nfile.append(str(round(np.log10(nH),4)))
                    zfile.append(str(round(np.log10(Z),4)))
                    
            logdiffn.append(np.log10(max(narray)) - np.log10(min(narray)))
            logdiffz.append(np.log10(max(zarray)) - np.log10(min(zarray)))
            mdn.append(str(round(np.log10(max(narray)),4) - round(np.log10(min(narray)),4)))
            mdz.append(str(round(np.log10(max(zarray)),4) - round(np.log10(min(zarray)),4)))
            
            uvb.append(Q)
            tnH.append(n)
            ion.append(str(ions_to_use[0]+'_'+ions_to_use[1]+'_'+ions_to_use[2]+'_'+ions_to_use[3]+'_'+ions_to_use[4]+'_'+
                           ions_to_use[5]+'_'+ions_to_use[6]+'_'+ions_to_use[7]))#+'_'+ions_to_use[8]))
        
from astropy.table import Table

mdn = np.array(mdn)
mdz = np.array(mdz)
nana = []
zaza = []
for i,j in zip(mdn,mdz):
    nana.append(str(round(float(i),4)))
    zaza.append(str(round(float(j),4)))
    
diff_res=Table([uvb,tnH, nq14,nq15,nq16,nq17,nq18,nq19,nq20,nqfg,nqp,nana,
                zq14,zq15,zq16,zq17,zq18,zq19,zq20,zqfg,zqp, zaza,ion],
               names=('Q','true_nH','N_Q14','N_Q15','N_Q16','N_Q17','N_Q18','N_Q19','N_Q20','N_FG20','N_P19','Max_diff_nH',
                     'Z_Q14','Z_Q15','Z_Q16','Z_Q17','Z_Q18','Z_Q19','Z_Q20','Z_FG20','Z_P19','Max_diff_Z','Ions_Used'))
diff_res.write('/Users/anshumanacharya/Downloads/cgm_uvb_master/mytry/full_{}ions.txt'.format(no), format='ascii.tab',overwrite='True')

print("Done")
