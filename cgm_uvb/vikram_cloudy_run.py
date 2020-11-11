from cgm_uvb.cloudy_run import write_input
from cgm_uvb.cloudy_run import run
from cgm_uvb.cloudy_run import store_table
from cgm_uvb.cloudy_run import cloudy_params_defaults


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