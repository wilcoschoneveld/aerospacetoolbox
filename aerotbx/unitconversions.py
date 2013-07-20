""" Aerospace Toolbox / unitconversions.py

Convert different types of data structures from one unit to another.

"""

import scipy as sp
from types import ListType
from aerotbx.utils import to_ndarray, from_ndarray

#acceleration
C_ACC = {'m/s^2': 1.0,
         'g': 9.80665,
         'gal': 0.01,
         'km/s^2': 1000.0,
         'in/s^2': 0.0254,
         'ft/s^2': 0.3048,
         'km/h-s': 1 / 3.6,
         'mph/s': 5280 * 0.3048 / 3600}
         
#angle
C_ANGLE = {'rad': 1.0,
         'deg': sp.pi/180,
         'rev': 2*sp.pi}

#angle velocity
C_ANVEL = {'rad/s': 1.0,
         'deg/s': sp.pi/180,
         'rpm': 2*sp.pi/60}

#angle acceleration
C_ANACC = {'rad/s^2': 1.0,
         'deg/s^2': sp.pi/180,
         'rpm/s': 2*sp.pi/60}

#density
C_DENS = {'kg/m^3': 1.0,
         'lbm/ft^3': 0.45359237 / 0.3048**3,
         'lbm/in^3': 0.45359237 / 0.0254**3,
         'slug/ft^3': 0.45359237 * 9.80665 / 0.3048**4}

#force
C_FORCE = {'N': 1.0,
         'kgf': 9.80665,
         'lbf': 0.45359237 * 9.80665,
         'dyne': 0.001 * 0.01}

#length
C_LEN = {'m': 1.0,
         'km': 1000.0,
         'mi': 5280 * 0.3048,
         'nmi': 1852.0,
         'in': 0.0254,
         'ft': 0.3048,
         'yd': 0.9144,
         'pt': 0.0254 / 72}

#mass
C_MASS = {'kg': 1.0,
         'ton': 1000.0,
         'lbm': 0.45359237,
         'slug': 0.45359237 * 9.80665 / 0.3048}

#pressure
C_PRESS = {'Pa': 1.0,
         'atm': 101325.0,
         'bar': 100000.0,
         'psi': 0.45359237 * 9.80665 / 0.0254**2,
         'psf': 0.45359237 * 9.80665 / 0.3048**2}

#temperature
C_TEMP = {'K': 1.0,
         'R': 5 / 9.0,
         'C': [lambda x: x + 273.15,
               lambda x: x - 273.15],
         'F': [lambda x: (x - 32)*(5/9.) + 273.15,
               lambda x: (x - 273.15)*(9/5.) + 32],
         'N': [lambda x: (x * 100/33.) + 273.15,
               lambda x: (x - 273.15) * 33. / 100]}

#velocity
C_VELO = {'m/s': 1.0,
         'ft/s': 0.3048,
         'km/s': 1000.0,
         'in/s': 0.0254,
         'km/h': 1 / 3.6,
         'mph': 5280 * 0.3048 / 3600,
         'kts': 1852.0 / 3600,
         'ft/min': 0.3048 / 60}

#energy
C_ENERGY = {'J': 1.0,
         'kwh': 1000.0 * 3600,
         'cal': 4.184,
         'ftlb': 0.3048 * 0.45359237 * 9.80665}

CNVTS = [C_ACC, C_ANGLE, C_ANVEL, C_ANACC, C_DENS, C_FORCE,
          C_LEN, C_MASS, C_PRESS, C_TEMP, C_VELO, C_ENERGY]

def convert(value, unit_from, unit_to):
    """convert value from one unit to another"""

    itype, value = to_ndarray(value)

    #loop through all converters and search for pairs
    for cnvt in CNVTS:
        if unit_from in cnvt and unit_to in cnvt:
            
            #convert to SI unit
            if type(cnvt[unit_from]) is ListType:
                value_si = cnvt[unit_from][0](value)
            else:
                value_si = value * cnvt[unit_from]

            #convert to target unit
            if type(cnvt[unit_to]) is ListType:
                value = cnvt[unit_to][1](value_si)
            else:
                value = value_si / cnvt[unit_to]

            return from_ndarray(itype, value)

    #no converter found
    raise ValueError('Could not convert between given units.')
