import numpy as np
import astropy.table as tab
from scipy.interpolate import  interp1D
import astropy.constants as const
import astropy.units as u



def hydrogen_photoionization_rate(wave):
     #Z=1.0 for Hydrogen




    return

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
