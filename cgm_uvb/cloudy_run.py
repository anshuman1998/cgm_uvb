# A simple code to run cloudy and store its output

import os
import numpy as np
import astropy.table as tab
import subprocess
from cgm_uvb.write_uvb_in_cloudy_format import write_uvb_in_cloudy_format

def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def run(cloudy_path, input_file):
    """
    :param cloudy_path: the path where your cloudy files are stored
    :param input_file: the input file (with full path) to run
    :return:
    """
    # find the original path
    basepath = os.getcwd()

    # go to the directory of input file
    os.chdir(os.path.dirname(input_file))

    # input file name
    file_name = os.path.basename(input_file)

    run_command =  cloudy_path + '/source/cloudy.exe'
    print('now running:', run_command, file_name)
    process = subprocess.Popen([run_command, file_name], stdout=subprocess.PIPE)
    process.stdout.read()

    # come back to the original path
    os.chdir(basepath)

    return


def write_input(file_name, *args, **kwargs):

    """
    :param file_name: the input filename where cloudy commands will be written
    :param args: these will be ions to store
    :param kwargs: the dict with keywords
        :keyword uvb (if KS18 then): uvb_scale, uvb_Q, z,
            hden_vary (True/False): if False give log_hden value
                else give log_hden1, log_hden2, log_hden_step
            log_metal,
            scale_He,
            stop_logNHI,
            constant_T,
            out_file_ext
    :return:
    """

    f = open(file_name, "w+")

    if kwargs['uvb'] == 'KS18':
        uvb_statement =  'TABLE {} redshift = {} [scale = {}] [Q = {}] \n'.format(
            kwargs['uvb'], kwargs['z'], kwargs['uvb_scale'], kwargs['uvb_Q'])
        f.write(uvb_statement)

    if kwargs['uvb'] == 'HM12':
        uvb_statement =  'TABLE {} redshift = {} [scale = {}] \n'.format(
            kwargs['uvb'], kwargs['z'], kwargs['uvb_scale'])
        f.write(uvb_statement)

    if kwargs['uvb'] == 'FG20':
        path_to_store = os.path.dirname(file_name)
        ebl_file_name =  'ebl_fg20_z{:.2f}.txt'.format(kwargs['z'])
        ebl_file_name_with_path = path_to_store + '/' + ebl_file_name
        fg20_file_path_and_name = os.getcwd() + '/paper_plots/fg20_fits_files' + '/FG20_EBL_z_{:.2f}.fits'.format(kwargs['z'])
        uvb_statement = 'TABLE SED \"{}\" \n'.format(ebl_file_name)

        if not os.path.exists(ebl_file_name_with_path):
            print('creating file :{}'.format(fg20_file_path_and_name))
            write_uvb_in_cloudy_format(fg20_file_path_and_name, FG20 = True, outfilename = ebl_file_name_with_path)

        x = tab.Table.read(fg20_file_path_and_name)
        data = tab.Table([12398.0 / x['Wave'] / 13.6057, x['Jnu'] * 2.99792e18 * np.pi * 4 / x['Wave']],
            names=('Ryd', 'nufnu'))
        ind = find_nearest_index(data['Ryd'], 1.000)
        norm_statement = 'nuf(nu) = {} [at {} Ryd]\n'.format(np.log10(data['nufnu'][ind] * kwargs['uvb_scale']),
            data['Ryd'][ind])

        f.write(uvb_statement)
        f.write(norm_statement)

    if kwargs['uvb'] == 'P19':
        path_to_store = os.path.dirname(file_name)
        ebl_file_name =  'ebl_p19_z{:.2f}.txt'.format(kwargs['z'])
        ebl_file_name_with_path = path_to_store + '/' + ebl_file_name
        p19_file_path_and_name = os.getcwd() + '/paper_plots/p19_ebl' + '/P19_EBL_z_{:.2f}.fits'.format(kwargs['z'])
        if not os.path.exists(p19_file_path_and_name):
            print('file {} does not exist, generate one'.format(p19_file_path_and_name))
        #    generate file ===> add this part later for now see if all files are there
        uvb_statement = 'TABLE SED \"{}\" \n'.format(ebl_file_name)

        if not os.path.exists(ebl_file_name_with_path):
            print('creating file: {}'.format(p19_file_path_and_name))
            write_uvb_in_cloudy_format(p19_file_path_and_name, P19 = True, outfilename = ebl_file_name_with_path)

        x = tab.Table.read(p19_file_path_and_name)
        data = tab.Table([12398.0 / x['Wave'] / 13.6057, x['Jnu'] * 2.99792e18 * np.pi * 4 / x['Wave']],
            names=('Ryd', 'nufnu'))
        ind = find_nearest_index(data['Ryd'], 1.000)
        norm_statement = 'nuf(nu) = {} [at {} Ryd]\n'.format(np.log10(data['nufnu'][ind] * kwargs['uvb_scale']),
            data['Ryd'][ind])


        f.write(uvb_statement)
        f.write(norm_statement)

    if kwargs['hden_vary'] == True:
        density_statement = 'hden {} vary \n'.format(kwargs['log_hden'])
        f.write(density_statement)
        if kwargs['sequential'] == True:
            variation_statement = 'grid sequential range from {} to {} with {} dex step \n'.format(
                kwargs['log_hden1'], kwargs['log_hden2'], kwargs['log_hden_step'])
            f.write(variation_statement)
        else:
            variation_statement = 'grid range from {} to {} with {} dex step \n'.format(
                kwargs['log_hden1'], kwargs['log_hden2'], kwargs['log_hden_step'])
            f.write(variation_statement)
    else:
        density_statement = 'hden {} \n'.format(kwargs['log_hden'])
        f.write(density_statement)

    if kwargs['abundances'] is not None:
        abd_statement = 'abundances \"{}\" \n'.format(kwargs['abundances'])
        f.write(abd_statement)

    metal_statement =  'metals {:.2f} log \n'.format(kwargs['log_metal'])
    f.write(metal_statement)

    if 'scale_He' in kwargs.keys():
        scale_He_statement ='element helium abundance {} linear \n'.format(kwargs['scale_He'])
        f.write(scale_He_statement)

    stop_statement = 'stop column density {} neutral H \n'.format(kwargs['stop_logNHI'])
    f.write(stop_statement)

    if kwargs['constant_T'] is not None:
        temp_statement =  'constant temperature, t={} K [linear] \n'.format(kwargs['constant_T'])
        f.write(temp_statement)

    # new line
    save_hydrogen = 'save hydrogen conditions \".hydro\" last no clobber \n'
    f.write(save_hydrogen)

    if 'out_file_ext' in kwargs.keys():
        out_file_extension = kwargs['out_file_ext']
    else:
        out_file_extension = '.spC'

    save_statement = 'save species column density \"{}\" no hash \n'.format(out_file_extension)
    f.write(save_statement)

    for ion in args:
        write_ion = "\"{}\" \n".format(ion)
        f.write(write_ion)

    f.write('end')

    f.close()

    return


# this is the part one needs to change if one wants to change the cloudy program
def cloudy_params_defaults(uvb_Q = 18, uvb_scale = 1, log_hden = [-4, -4], hden_vary=True, uvb = 'KS18', z=0.2,
                           T = None, metal = -1, stop_NHI = 15, abundances = 'solar_GASS10.abn', sequential = False):

    cloudy_params = {'uvb': uvb, 'z' : z, 'uvb_scale': uvb_scale, 'uvb_Q' : uvb_Q,
                     'hden_vary' : hden_vary,
                     'log_metal': metal,
                     'constant_T': T,
                     'stop_logNHI': stop_NHI,
                     'scale_He': 0.081632653,
                     'abundances' : abundances,
                     'sequential': sequential}
    print(cloudy_params)

    if hden_vary :
        cloudy_params['log_hden'] = log_hden[0]
        cloudy_params['log_hden1'] = log_hden[0]
        cloudy_params['log_hden2'] = log_hden[1]
        cloudy_params['log_hden_step'] = log_hden[2]
    else:
        cloudy_params['log_hden'] = log_hden[0]


    ions = ["H", "H+",
            "He", "He+", "He+2",
            "C", "C+", "C+2", "C+3", "C+4", "C+5",
            "N", "N+", "N+2", "N+3", "N+4",
            "O", "O+", "O+2", "O+3", "O+4", "O+5", "O+6", "O+7",
            "S", "S+", "S+2", "S+3", "S+4", "S+5",
            "Si", "Si+", "Si+2", "Si+3", "Si+4",
            "Mg", "Mg+", "Mg+2",
            "Ne", "Ne+", "Ne+2", "Ne+3", "Ne+4", "Ne+5", "Ne+6", "Ne+7", "Ne+8",
            "Fe", "Fe+", "Fe+2",
            "Na", "Na+", "Na+2"]


    return ions, cloudy_params


def store_table(ions, output_file, fits_filename = None):
    # get hydrogen density array
    hydro_file = output_file.split('.')[0] + '.hydro'
    hydro =  tab.Table.read(hydro_file, format = 'ascii')
    hden_array = np.unique(hydro['HDEN'])

    # read the cloudy output files
    cloudy_output = tab.Table.read(output_file, format = 'ascii')

    for old_name, new_name in zip(cloudy_output.colnames, ions):
        cloudy_output.rename_column (old_name, new_name)

    cloudy_output.add_column(hden_array, name = 'hden' )

    if fits_filename == None:
        fits_filename = output_filename.split('.')[0] + 'fits'

    cloudy_output.write (fits_filename, overwrite = True)

    return


"""
Example run : 
#----give this
uvb_Q=20
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

