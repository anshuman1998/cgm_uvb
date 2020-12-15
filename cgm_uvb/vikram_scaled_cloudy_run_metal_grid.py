#------------------
# for reducing  numpy threads
import os
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
#-----------------
from cgm_uvb.cloudy_run import write_input
from cgm_uvb.cloudy_run import run
from cgm_uvb.cloudy_run import store_table
from cgm_uvb.cloudy_run import cloudy_params_defaults
import multiprocessing as mp
import numpy as np
from cgm_uvb.write_uvb_in_cloudy_format import write_uvb_in_cloudy_format
import astropy.table as tab
from cgm_uvb.find_uvb_scaling import find_uvb_scaling


def uvb_files(file_path, **kwargs):
    if kwargs['uvb'] == 'FG20':
        path_to_store = file_path
        ebl_file_name =  'ebl_fg20_z{:.2f}.txt'.format(kwargs['z'])
        ebl_file_name_with_path = path_to_store + '/' + ebl_file_name
        fg20_file_path_and_name = os.getcwd() + '/paper_plots/fg20_fits_files' + '/FG20_EBL_z_{:.2f}.fits'.format(kwargs['z'])
        uvb_statement = 'TABLE SED \"{}\" \n'.format(ebl_file_name)

        if not os.path.exists(ebl_file_name_with_path):
            norm_statement = write_uvb_in_cloudy_format(fg20_file_path_and_name, FG20 = True, outfilename = ebl_file_name_with_path)
        else:
            print('file exists', ebl_file_name_with_path)

        print(uvb_statement)

    if kwargs['uvb'] == 'P19':
        path_to_store = file_path
        ebl_file_name =  'ebl_p19_z{:.2f}.txt'.format(kwargs['z'])
        ebl_file_name_with_path = path_to_store + '/' + ebl_file_name
        p19_file_path_and_name = os.getcwd() + '/paper_plots/p19_ebl' + '/P19_EBL_z_{:.2f}.fits'.format(kwargs['z'])
        if not os.path.exists(p19_file_path_and_name):
            print('file {} does not exist, generate one'.format(p19_file_path_and_name))
        #    generate file ===> add this part later for now see if all files are there
        uvb_statement = 'TABLE SED \"{}\" \n'.format(ebl_file_name)

        if not os.path.exists(ebl_file_name_with_path):
            norm_statement = write_uvb_in_cloudy_format(p19_file_path_and_name, P19 = True, outfilename = ebl_file_name_with_path)
        else:
            print('file exists', ebl_file_name_with_path)

        print(uvb_statement)

    return

def run_parallel(logZ, uvb_Q, uvb, scale_value):
    # for vikram
    cloudy_path = '/home/vikram/c17.02'
    fname = (logZ+4)*100
    input_File = '/mnt/quasar2/vikram/cloudy_run/rescaled_metal_NH15/try_{}_Q{}_Z{:.0f}.in'.format(uvb, uvb_Q, fname)
    print(uvb, 'Q=', uvb_Q, 'Z=', logZ)

    # write input file and run cloudy
    ions, params = cloudy_params_defaults(uvb = uvb, uvb_Q=uvb_Q, log_hden=[-6, -2, 0.02], stop_NHI = 15, T = None,
        metal = logZ, uvb_scale = scale_value, sequential = True)
    write_input(input_File, *ions, **params)
    run(cloudy_path=cloudy_path, input_file=input_File)

    # write output tables
    output_filename = input_File.split('.in')[0] +  '.spC'
    fits_filename = input_File.split('.in')[0] + '.fits'
    store_table(ions=ions, output_file=output_filename, fits_filename=fits_filename)

    return

# runnning in parallel
#uvb_array= [14, 15, 16, 17, 18, 19, 20]
logZ_array = np.around(np.arange(-3, 1, 0.05), decimals = 2)
uvb = ['KS18', 'HM12',  'P19', 'FG20']
uvb_Q = [14, 15, 16, 17, 18, 19, 20]

scale_ks  =  []
for q_model in uvb_Q:
    scaling_factor  =  find_uvb_scaling(uvb = 'KS18', uvb_Q = q_model)
    scale_ks.append(scaling_factor)
    print(scale_ks, q_model)

logZ = []
uvb_models =[]
the_Q_values = []
scaling_factor_array = []
for background in uvb:
    if background == 'KS18':
        for q, the_scale in zip(uvb_Q, scale_ks):
            for metal in logZ_array:
                uvb_models.append(background)
                the_Q_values.append(q)
                logZ.append(metal)
                scaling_factor_array.append(the_scale)

    else:
        q = 18
        scaling_factor_uvb  =  find_uvb_scaling(uvb = background)
        print(scaling_factor_uvb, background)
        for metal in logZ_array:
            uvb_models.append(background)
            the_Q_values.append(q)
            logZ.append(metal)
            scaling_factor_array.append(scaling_factor_uvb)


#-----write uvb fg and hm in cloudy format first
path = '/mnt/quasar2/vikram/cloudy_run/rescaled_metal_NH15'

kwagrs = {'uvb' : 'P19', 'z' : 0.2}
uvb_files(path, **kwagrs)

kwagrs = {'uvb' : 'FG20', 'z' : 0.2}
uvb_files(path, **kwagrs)


pool = mp.Pool(processes=6)
results = [pool.apply_async(run_parallel, args=(Z, Q, mod, scale_val,)) for  Z, Q, mod, scale_val
           in zip(logZ, the_Q_values, uvb_models, scaling_factor_array)]
output = [p.get() for p in results]

