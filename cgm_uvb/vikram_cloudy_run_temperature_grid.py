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

def run_parallel(logTemp, uvb_Q = 18):
    # for vikram
    cloudy_path = '/home/vikram/c17.02'
    input_File = '/home/vikram/cloudy_run/temperature_NH15/try_Q{}_T{:.0f}.in'.format(uvb_Q, logTemp*100)

    # write input file and run cloudy
    ions, params = cloudy_params_defaults(uvb_Q=uvb_Q, log_hden=[-6, -2, 0.01], stop_NHI = 15, T = 10**logTemp,
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
logTemp_array = np.arange(3.64, 6.02, 0.04)

pool = mp.Pool(processes=4)
results = [pool.apply_async(run_parallel, args=(T, )) for  T in logTemp_array]
output = [p.get() for p in results]

"""
uvb_Q=19
cloudy_path = '/home/vikram/c17.02'
input_File = '/home/vikram/cloudy_run/try.in'

# write input file and run cloudy
ions, params = cloudy_params_defaults(uvb_Q=uvb_Q, log_hden= [-5, -3, 1])
write_input(input_File, *ions, **params)
run(cloudy_path= cloudy_path, input_file= input_File)

# write output tables
output_filename =  input_File.split('.in')[0] + '.spC'
fits_filename = input_File.split('.in')[0] + '_Q{}'.format(uvb_Q) + '.fits'
store_table(ions= ions, output_file= output_filename, fits_filename= fits_filename)
"""