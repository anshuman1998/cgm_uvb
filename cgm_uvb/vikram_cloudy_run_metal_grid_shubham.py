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
import astropy.table as tab


def run_parallel(logZ):
    uvb = 'KS18'
    uvb_Q = 18
    # for vikram
    cloudy_path = '/home/vikram/c17.02'
    fname = (logZ+4)*100
    input_File = '/mnt/quasar2/vikram/cloudy_run/shubham/try_{}_Q{}_Z{:.0f}.in'.format(uvb, uvb_Q, fname)
    print(uvb, 'Q=', uvb_Q, 'Z=', logZ)

    # write input file and run cloudy
    ions, params = cloudy_params_defaults(uvb = uvb, uvb_Q=uvb_Q, z=0.39047, log_hden=[-6, -2, 0.02],
        stop_NHI = 18.85, T = None, metal = logZ, sequential = True)
    write_input(input_File, *ions, **params)
    run(cloudy_path=cloudy_path, input_file=input_File)

    # write output tables
    output_filename = input_File.split('.in')[0] +  '.spC'
    fits_filename = input_File.split('.in')[0] + '.fits'
    store_table(ions=ions, output_file=output_filename, fits_filename=fits_filename)

    return

# runnning in parallel
logZ_array = np.around(np.arange(-3, 1, 0.05), decimals = 2)

logZ = []

for metal in logZ_array:
    logZ.append(metal)

#-----write uvb fg and hm in cloudy format first
path = '/mnt/quasar2/vikram/cloudy_run/shubham'


pool = mp.Pool(processes=20)
results = [pool.apply_asyn(run_parallel, args=(Z,)) for  Z in zip(logZ)]
output = [p.get() for p in results]

