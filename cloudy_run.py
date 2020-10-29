# A simple code to run cloudy and store its output

import os
import numpy as np
import subprocess



def run(cloudy_path, input_file):
    """
    :param cloudy_path: the path where your cloudy files are stored
    :param input_file: the input file to run
    :return:
    """

    run_command =  cloudy_path + '/source/cloudy.exe'
    print(run_command, input_file)
    process = subprocess.Popen([run_command, input_file], stdout=subprocess.PIPE)
    process.stdout.read()

    return


def write_input():
    xy =3

    return

"""
UVB_Q = 20
dirname = '/home/vikram/cloudy_run/final'
run = '/home/vikram/c17.02/source/cloudy.exe'
hden = np.arange(-5, -2.999, 1)  # Hydrogen density grid
z = -1  # in log scale
stcolden = 14
if not os.path.exists(dirname):
    os.makedirs(dirname)

os.chdir(dirname)

filena = dirname + '/prog' + '{:.0f}'.format(UVB_Q) + '_Oct5.in'
f = open(filena, "w+")
f.write("TABLE KS18 redshift = 0.2 [scale= 1][Q=" + "{:.0f}".format(UVB_Q)
        + "] \nhden -4 vary \ngrid range from -5 to -3 with 1 dex steps \nmetals "
        + "{:.2f}".format(z) + " log \nelement helium abundance 0.081632653 linear \nstop column density "
        + "{:.1f}".format(
    stcolden) + " neutral H \nconstant temperature, t=1e4 K [linear] \nsave species column density \".spC\" no hash \n\"C+\" \n\"C+2\" \n\"C+3\" \n\"N+\" \n\"N+2\" \n\"N+3\" \n\"N+4\" \n\"O\" \n\"O+\"  \n\"O+2\" \n\"O+3\" \n\"O+4\" \n\"O+5\" \n\"S+3\" \n\"S+4\" \n\"S+5\" \n\"Si+\" \n\"Si+2\" \n\"Si+3\" \nend")
f.close()
process = subprocess.Popen([run, "prog" + "{:.0f}".format(UVB_Q) + "_Oct5.in"],
                           stdout=subprocess.PIPE)
process.stdout.read()

outfilename = filena.split('.in')[0] + '.spC'
data_20 = np.genfromtxt(outfilename)


print("Done!!")
"""
