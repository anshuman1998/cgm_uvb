import numpy as np
import astropy.table as tab
from scipy.interpolate import  interp1D
import astropy.constants as const
import astropy.units as u



def photoionization_cross_section(nu, element = 'HI'):
    """
    :param nu:
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

    sigma_0 = 6.304e18/Z**2
    y = h*nu/IH/Z**2
    print(y)
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

    nu_array =  c*1e10 / wavelength
    integral_array  =  np.pi * 4 * intensity/(nu_array*h)

    integration  = np.trapz (integral_array, nu_array)

    return integration

def get_uvb(uvb, Q):


    uvb = tab.Table.read(file_name)
    wave =  np.array (uvb['wave']) # in angstrom
    jnu = np.array (uvb['Jnu']) # in standard units

    return wave, jnu
