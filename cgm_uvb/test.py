#------------------
# for reducing  numpy threads
import os
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
#-----------------

import time
import numpy as np
import multiprocessing as mp

def run_parallel(x):
    time.sleep(2)
    f = open("vikram{}.txt".format(x+1), "a")
    f.write("Now the file has more content!")
    f.close()
    print('process = ', x+ 1 )

    return




num_proc= 56
pool = mp.Pool(processes=num_proc)
results = [pool.apply_async(run_parallel, args=(i,))for  i in np.arange(num_proc)]
output = [p.get() for p in results]