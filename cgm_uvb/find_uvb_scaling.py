import numpy as np
import astropy.table as tab
from scipy.interpolate import  interp1d
import astropy.constants as const
import astropy.units as u


#todo make this code work for an array of frequencies
def photoionization_cross_section(nu, element = 'HI'):
    """
    :param nu: works on only one number
    :param element: {'HI', 'HeII'}
    :return:
    """
    h = (const.h).to(u.eV* u.s).value # h in eV s
    IH = 13.602 # eV ionization potential for H

    if element == 'HI':
        Z = 1
    # check the helium part later
    if element == 'HeII':
        Z = 2

    sigma_0 = 6.304e-18/Z**2
    y = h*nu/IH/Z**2
    #print(y)
    if y > 1:
        x = (y - 1)**0.5
        cross_section = sigma_0 * y**(-4) * np.exp(4 - 4*np.arctan(x)/x)/(1-np.exp(-2*np.pi/x))
    else:
        cross_section = 0

    return cross_section

def get_HI_photoionization_rate(wavelength, intensity):

    # just to simple trapazoid integration
    c = const.c.value # in m/s by default
    h = (const.h).to(u.erg * u.s).value # h in erg s

    nu_array = c*1e10/wavelength # wavelength in Angstrom

    sigma_array = []
    for nu in nu_array:
        sigma = photoionization_cross_section(nu)
        sigma_array.append(sigma)

    sigma_array= np.array(sigma_array)

    integral_array  =  np.pi * 4 * intensity * sigma_array /(nu_array*h)

    integration  = np.trapz (integral_array, nu_array)

    return -1.0*integration

def get_uvb(uvb='KS18', Q =18, z = 0.2):
    uvb_file_path =  '/home/vikram/cgm_uvb/cgm_uvb/paper_plots'
    if uvb is 'KS18':
        file_name = uvb_file_path +'/ks19_ebl' + '/KS19_EBL_z_{:.1f}.fits'.format(z)
        q_name = 'Q{}'.format(Q)
        uvb = tab.Table.read(file_name)
        wave = np.array(uvb['Wave'])  # in angstrom
        jnu = np.array(uvb[q_name])  # in standard units

    if uvb is 'FG20':
        file_name = uvb_file_path + '/fg20_fits_files' + '/FG20_EBL_z_{:.2f}.fits'.format(z)
        uvb = tab.Table.read(file_name)
        wave = np.array(uvb['Wave'])  # in angstrom
        jnu = np.array(uvb['Jnu'])  # in standard units

    if uvb is 'P19':
        file_name = uvb_file_path  + '/p19_ebl'+ '/P19_EBL_z_{:.2f}.fits'.format(z)
        uvb = tab.Table.read(file_name)
        wave = np.array(uvb['Wave'])  # in angstrom
        jnu = np.array(uvb['Jnu'])  # in standard units

    if uvb is 'HM12':
        file_name = uvb_file_path  + '/hm12_ebl'+ '/HM12_EBL_z_{:.2f}.fits'.format(z)
        uvb = tab.Table.read(file_name)
        wave = np.array(uvb['Wave'])  # in angstrom
        jnu = np.array(uvb['Jnu'])  # in standard units


    potential_HI = ((13.602*u.eV).to(u.J)).value
    wave_for_ionization = 1e10*(const.h).value * (const.c).value/potential_HI
    print((const.h).value)
    print(wave_for_ionization)
    new_wave = wave[wave<wave_for_ionization]
    new_jnu = jnu[wave<wave_for_ionization]


    return new_wave, new_jnu


w, j= get_uvb()
gamma_HI=get_HI_photoionization_rate(w, j)
print(gamma_HI)
