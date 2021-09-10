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

def run_parallel(logZ, stop_NHI, uvb_Q, uvb):
    # for vikram
    cloudy_path = '/home/vikram/c17.02'
    fname = (logZ+4)*100
    input_File = '/scratch/vikram/cloudy_run/metal_NH{:0.0f}/try_{}_Q{}_Z{:.0f}.in'.format(stop_NHI, uvb, uvb_Q, fname)
    print(uvb, 'Q=', uvb_Q, 'Z=', logZ)

    # write input file and run cloudy
    ions, params = cloudy_params_defaults(uvb = uvb, uvb_Q=uvb_Q, log_hden=[-6, -2, 0.02], stop_NHI = stop_NHI, T = None, metal = logZ,
                                          sequential = True)
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
NHI = [14, 18]

logZ = []
log_NHI = []
uvb_models =[]
the_Q_values = []
for column in NHI:
    for background in uvb:
        if background=='KS18':
            for q in uvb_Q:
                for metal in logZ_array:
                    uvb_models.append(background)
                    the_Q_values.append(q)
                    logZ.append(metal)
                    log_NHI.append(column)
        else:
            q = 18
            for metal in logZ_array:
                uvb_models.append(background)
                the_Q_values.append(q)
                logZ.append(metal)
                log_NHI.append(column)

#-----write uvb fg and hm in cloudy format first
#path = '/mnt/quasar2/vikram/cloudy_run/metal_NH17_new'
path  = '/scratch/vikram/cloudy_run/metal_NH14'

kwagrs = {'uvb' : 'P19', 'z' : 0.2}
uvb_files(path, **kwagrs)

kwagrs = {'uvb' : 'FG20', 'z' : 0.2}
uvb_files(path, **kwagrs)

#print(len(logZ))

pool = mp.Pool(processes=168)
results = [pool.apply_async(run_parallel, args=(Z, stop_column, Q, mod,)) for  Z, stop_column, Q, mod in zip(logZ, log_NHI, the_Q_values, uvb_models)]
output = [p.get() for p in results]


