# pylint: disable=C0103,W0142
""" Aerospace Toolbox / matlab.py

Wrappers for the matlab syntax from http://www.mathworks.nl/help/aerotbx/

"""

import aerotbx

#Environment

def atmosisa(height):
    """
    Use International Standard Atmosphere Model.

    http://www.mathworks.nl/help/aerotbx/ug/atmosisa.html
    """
    _, T, P, rho, a = aerotbx.stdatmos(geop=height)
    return T, a, rho, P

def atmoslapse(height, g, heatRatio, characteristicGasConstant,
               lapseRate, heightTroposphere, heightTropopause,
               density0, pressure0, temperature0):
    """
    Use Lapse Rate Atmosphere Model.

    http://www.mathworks.nl/help/aerotbx/ug/atmoslapse.html
    """
    mdl = aerotbx.stdmodel(R=characteristicGasConstant, gamma=heatRatio,
                           g0=g, T0=temperature0, P0=pressure0,
                           lapserate=[-lapseRate, 0],
                           layers=[0, heightTroposphere, heightTropopause])
    _, T, P, rho, a = aerotbx.stdatmos(model=mdl, geop=height)
    return T, a, P, rho

def atmospalt(pressure):
    """
    Calculate pressure altitude based on ambient pressure.

    http://www.mathworks.nl/help/aerotbx/ug/atmospalt.html
    """
    h = aerotbx.stdatmos(P=pressure)[0]
    return h

def geoidegm96(lat, lon):
    """
    Implement the EGM96 Geopotential Model

    http://www.mathworks.nl/help/aerotbx/ug/geoidegm96.html
    """
    N = aerotbx.geoidheight(lat, lon)
    return N

#Gas Dynamics

def flowisentropic(gamma, flow, mtype='mach'):
    """
    Calculate isentropic flow ratios.

    http://www.mathworks.nl/help/aerotbx/ug/flowisentropic.html
    """
    M, T, P, rho, area = aerotbx.flowisentropic(**{'gamma':gamma, mtype:flow})
    return M, T, P, rho, area

def flownormalshock(gamma, normal_shock_relations, mtype='mach'):
    """
    Calculate normal shock relations.

    http://www.mathworks.nl/help/aerotbx/ug/flownormalshock.html
    """
    M, M2, T, P, rho, P0, P1 = aerotbx.flownormalshock(**{'gamma':gamma,
                                                mtype:normal_shock_relations})
    return M, P, T, rho, M2, P0, P1

def flowprandtlmeyer(gamma, prandtlmeyer_array, mtype='mach'):
    """
    Calculate Prandtl-Meyer functions for expansion waves.

    http://www.mathworks.nl/help/aerotbx/ug/flowprandtlmeyer.html
    """
    M, V, U = aerotbx.flowprandtlmeyer(**{'gamma':gamma,
                                       mtype:prandtlmeyer_array})
    return M, V, U

#Unit Conversions

def conv(valuesToConvert, inputUnits, outputUnits):
    """
    Convert from units to desired units.

    http://www.mathworks.nl/help/aerotbx/unit-conversions-1.html
    """
    convertedValues = aerotbx.convert(valuesToConvert, inputUnits,
                                      outputUnits)
    return convertedValues
