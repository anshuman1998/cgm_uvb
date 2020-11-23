import astropy.table as tab
from astropy.io import ascii
import numpy as np


def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def write_uvb_in_cloudy_format (uvb_table, KS19_model = None, FG20 =False, P19 = False, outfilename = 'ebl.test'):
    x = tab.Table.read(uvb_table)
    if KS19_model is not None:
        data = tab.Table([12398.0 / x['Wave'] / 13.6057, x['Q14'] * 2.99792e18 * np.pi * 4 / x['Wave']],
                         names=('Ryd', 'nufnu'))

    if FG20:
        data = tab.Table([12398.0 / x['Wave'] / 13.6057, x['Jnu'] * 2.99792e18 * np.pi * 4 / x['Wave']],
                         names=('Ryd', 'nufnu'))

    if P19:
        data = tab.Table([12398.0 / x['Wave'] / 13.6057, x['Jnu'] * 2.99792e18 * np.pi * 4 / x['Wave']],
                         names=('Ryd', 'nufnu'))


    ascii.write(data, outfilename, overwrite=True, format='no_header')

    ind = find_nearest_index(data['Ryd'], 1.000)
    statement = 'nuf(nu) = {} [at {} Ryd]\n'.format(np.log10(data['nufnu'][ind]), data['Ryd'][ind])
    print(statement)

    with open(outfilename) as f:
        lines = f.readlines()

    # lines # ['This is the first line.\n', 'This is the second line.\n']
    lines[0] = lines[0].split('\n')[0] + " nufnu\n"
    # lines # ["This is the line that's replaced.\n", 'This is the second line.\n']

    with open(outfilename, "w") as f:
        f.writelines(lines)

    f.close()

    return statement