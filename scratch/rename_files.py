import os
import glob

path  = '/home/vikram/cloudy_run/rescaled_hybrid_NH15/'
fnames = glob.glob(path+ '/t*.*')


for f in fnames:
    x= f.split('.')
    y = x[0]+ '0.' + x[1]

    num = x[0].split('logT')[-1]
    if num == '52':
        y = x[0] + '5.' + x[1]
    if num == '58':
        y =  x[0].split('58')[0] + '575.' + x[1]

    os.rename(f, y)
    print('renamed', f, 'into', y)

