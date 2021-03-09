import astropy.table as tab
import numpy as np

def find_max_o6_frac(filename, Z = 0.1, log_o_abd = -3.3098):

    data  = tab.Table.read(filename)
    tot_o = (data['H']+ data['H+'])*Z*10**log_o_abd
    frac_array = data['O+5']/tot_o
    max_arg =np.argmax(frac_array)
    print(frac_array[max_arg], data['hden'][max_arg], filename.split('/')[-1])

    return


path = '/home/vikram/cloudy_run/metal_NH15_new'

uvb = ['KS18', 'HM12',  'P19', 'FG20']
uvb_Q = [14, 15, 16, 17, 18, 19, 20]

uvb_models =[]
the_Q_values = []
for background in uvb:
    if background == 'KS18':
        for q in uvb_Q:
            uvb_models.append(background)
            the_Q_values.append(q)
    else:
        q = 18
        uvb_models.append(background)
        the_Q_values.append(q)


logZ = -1
fname = (logZ+4)*100

for uvb_name, q_name in zip(uvb_models, the_Q_values):
    filename = 'try_{}_Q{}_Z{:.0f}.fits'.format(uvb_name, q_name, fname)
    file_with_path = path + '/' + filename
    find_max_o6_frac(filename= file_with_path)

logT = 5.5
print('for hot gas')
path = '/home/vikram/cloudy_run/hybrid_NH15'
fname = logT * 100
for uvb_name, q_name in zip(uvb_models, the_Q_values):

    filename = 'try_{}_Q{}_logT{:.0f}.fits'.format(uvb_name, q_name, fname)
    file_with_path = path + '/' + filename
    find_max_o6_frac(filename= file_with_path)
