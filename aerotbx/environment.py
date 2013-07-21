# pylint: disable=C0103,W0603,E0611,E1103
""" Aerospace Toolbox / environment.py

Analyse the environment with different standard atmosphere models and
evaluation of the EGM96 geopotential model.

"""

import scipy as sp
from scipy.interpolate import RectSphereBivariateSpline
from pkg_resources import resource_string
from types import DictType
from aerotbx.utils import to_ndarray, from_ndarray

_EGM96 = None

def _loadEGM96():
    """load the EGM96 geoid model into a spline object"""
    
    #load the data resource file into a string
    flc = resource_string(__name__, "data/egm96.dac")

    #setup basic coordinates
    lon = sp.linspace(0, 2*sp.pi, 1440, False)
    lat = sp.linspace(0, sp.pi, 721)

    #parse the raw data string
    data = sp.fromstring(flc, sp.dtype(sp.int16).newbyteorder("B"),
        1038240).reshape((lat.size, lon.size)) / 100.0

    #interpolate the bad boy
    lut = RectSphereBivariateSpline(lat[1: -1], lon, data[1: -1],
        pole_values = (sp.mean(data[1]), sp.mean(data[-1])))

    return lut

def stdmodel(**params):
    """
    Define a model for use in evaluating the standard atmosphere.
    
    The user can define as many parameters as needed, all other values
    are default as defined in the International Standard Atmosphere.
    
    Parameters
    ----------
    R : float
        Specific gas constant [J/kg/K].
    gamma : float
        Specific heat ratio.
    radius : float
        Radius of planet [m].
    g0 : float
        Standard gravity at zero altitude [m/s^2].
    T0 : float
        Temperature at zero altitude [K].
    P0 : float
        Pressure at zero altitude [K].
    lapserate : array_like
        Rate of temperature change with height [K/m]. The lapse rate should
        be defined in line with geopotential altitude.
    layers : array_like
        Height layers accompanying the lapse rate input [m]. Should be
        defined as geopotential altitude and contain atleast one more
        value as the given lapse rate to account for layer base and top
        values.
    
    Returns
    -------
    model : dict
        Dictionary containing all given values
    """
    
    return params

def stdatmos(**altitude):
    """
    Evaluate the standard atmosphere at any given altitude.
    
    This function allows input of a single variable to calculate
    the atmospheric properties at different altitudes. The function
    can work with different types of standard models. The default model
    values are set as defined by the International Standard Atmosphere.
    
    Parameters
    ----------
    model : dict, optional
        A standard atmosphere model as obtained from stdmodel. Partial
        models are allowed. The remaining model values default to
        the International Standard Atmosphere.
    h or geom : array_like
        Geometrical altitude [meters].
    geop : array_like
        Geopotential altitude [meters].
    abs : array_like
        Absolute altitude [meters].
    T : array_like
        Temperature altitude [K].
    P : array_like
        Pressure altitude [Pa].
    rho : array_like
        Density altitude [kg/m^3].
        
    Returns
    -------
    out : (h, T, P, rho, a)
        Tuple of geometrical altitude, temperature, pressure, density
        and speed of sound at given altitudes.
    
    Notes
    -----
    This function assumes a continues lapserate below 0 altitude and above
    the top layer, which allows for extrapolation outside the specified
    region (0 to 86km in ISA). Temperature altitude is obtained as the
    first altitude from 0 where the specified temperature exists.
        
    See Also
    --------
    stdmodel
    
    Examples
    --------
    >>> stdatmos(P=[1e5, 1e4, 1e3])[0]
    [110.8864127251899, 16221.007939493587, 31207.084373790043]
    >>> stdatmos(h=sp.linspace(-2000, 81000))
    (array, array, array, array, array)
    """
    
    #pop atmospherical model from input
    model = altitude.pop("model", {})

    #check if model is a dictionary
    if type(model) is not DictType:
        raise Exception("Custom atmosphere model is incompatible.")

    #check if a single remaining input exists
    if len(altitude) is not 1:
        raise Exception("Function needs exactly one altitude input.")

    #pop the altitude input
    mtype, alt = altitude.popitem()

    #check if the altitude input type is valid
    if mtype not in ["h", "geom", "geop", "abs", "T", "P", "rho"]:
        raise Exception("The altitude input should be a valid input type.")

    #convert the input to numpy arrays
    itype, alt = to_ndarray(alt)

    #model values
    R = model.get("R", 287.053) #gas constant [J/kg/K] (air)
    gamma = model.get("gamma", 1.4) #specific heat ratio [-] (air)

    g = model.get("g0", 9.80665) #gravity [m/s^2] (earth)
    radius = model.get("radius", 6356766.0) #earth radius [m] (earth)

    Tb = model.get("T0", 288.15) #base temperature [K]
    Pb = model.get("P0", 101325.0) #base pressure [Pa]
    
    #model lapse rate and height layers
    Hb = sp.array([0, 11, 20, 32, 47, 51, 71, sp.inf], sp.float64) * 1000
    Lr = sp.array([-6.5, 0, 1, 2.8, 0, -2.8, -2], sp.float64) * 0.001

    Hb = model.get("layers", Hb) #layer height [km]
    Lr = model.get("lapserate", Lr) #lapse rate [K/km]

    #preshape solution arrays
    T = sp.ones(alt.shape, sp.float64) * sp.nan
    P = sp.ones(alt.shape, sp.float64) * sp.nan

    #define the height array
    if mtype in ["h", "geom"]:
        h = alt * radius / (radius + alt)
    elif mtype is "geop":
        h = alt
    elif mtype is "abs":
        h = alt - radius
    else:
        h = sp.ones(alt.shape, sp.float64) * sp.nan

    for lr, hb, ht in zip(Lr, Hb[:-1], Hb[1:]):
        #calculate the temperature at layer top
        Tt = Tb + lr*(ht-hb)
        
        if mtype is "T":
            #break the loop if there are no nans in the solution array
            if not sp.isnan(h).any():
                break

            #select all temperatures in current layer
            if lr == 0:
                sel = (alt == Tb)
            else:
                s = sp.sign(lr)
                bot = -sp.inf if hb == 0 else Tb*s
                top = sp.inf if ht == Hb[-1] else Tt*s
                sel = sp.logical_and(alt*s >= bot, alt*s < top)

            #only select when not already solved
            sel = sp.logical_and(sel, sp.isnan(h))

            #temperature is given as input
            T[sel] = alt[sel]

            #solve for height and pressure
            if lr == 0:
                h[sel] = hb
                P[sel] = Pb
            else:
                h[sel] = hb + (1.0/lr)*(T[sel] - Tb)
                P[sel] = Pb * (T[sel] / Tb)**(-g/(lr*R))
                
        elif mtype in ["P", "rho"]:
            #choose base value as pressure or density
            vb = Pb if mtype is "P" else Pb/(R*Tb)

            #select all input values below given pressure or density
            sel = alt <= (sp.inf if hb == 0 else vb)

            #break if nothing is selected
            if not sel.any():
                break

            #solve for temperature and height
            if lr == 0:
                T[sel] = Tb
                h[sel] = hb - sp.log(alt[sel]/vb)*R*Tb/g
            else:
                x = g if mtype is "P" else (lr*R + g)
                T[sel] = Tb * (alt[sel]/vb)**(-lr*R/x)
                h[sel] = hb + (T[sel] - Tb) / lr

            #pressure is given as input
            P[sel] = alt[sel] if mtype is "P" else alt[sel]*R*T[sel]

        else:
            #select all height values above layer base
            sel = h >= (-sp.inf if hb == 0 else hb)

            #break if nothing is selected
            if not sel.any():
                break

            #solve for temperature and pressure
            if lr == 0:
                T[sel] = Tb
                P[sel] = Pb * sp.exp((-g/(R*Tb))*(h[sel]-hb))
            else:
                T[sel] = Tb + lr * (h[sel] - hb)
                P[sel] = Pb * (T[sel] / Tb)**(-g/(lr*R))

        #update pressure base value
        if lr == 0:
            Pb *= sp.exp((-g/(R*Tb))*(ht - hb))
        else:
            Pb *= (Tt / Tb)**(-g/(lr*R))

        #update temperature base value
        Tb = Tt

    #convert geopotential altitude to geometrical altitude
    h *= radius / (radius - h)

    #density
    rho = P / (R*T)

    #speed of sound
    a = sp.sqrt(gamma * R * T)
    
    return from_ndarray(itype, h, T, P, rho, a)

def geoidheight(lat, lon):
    """
    Calculate geoid height using the EGM96 Geopotential Model.

    Parameters
    ----------
    lat : array_like
        Lateral coordinates [degrees]. Values must be -90 <= lat <= 90.
    lon : array_like
        Longitudinal coordinates [degrees]. Values must be 0 <= lon <= 360.

    Returns
    -------
    out : array_like
        Geoidheight [meters]

    Examples
    --------
    >>> geoidheight(30, 20)
    25.829999999999995
    >>> geoidheight([30, 20],[40, 20])
    [9.800000000000002, 14.43]
    """

    global _EGM96

    #convert the input value to array
    itype, lat = to_ndarray(lat)
    itype, lon = to_ndarray(lon)

    if lat.shape != lon.shape:
        raise Exception("Inputs must contain equal number of values.")

    if (lat < -90).any() or (lat > 90).any() or not sp.isreal(lat).all():
        raise Exception("Lateral coordinates must be real numbers" \
            " between -90 and 90 degrees.")

    if (lon < 0).any() or (lon > 360).any() or not sp.isreal(lon).all():
        raise Exception("Longitudinal coordinates must be real numbers" \
            " between 0 and 360 degrees.")

    #if the model is not loaded, do so
    if _EGM96 is None:
        _EGM96 = _loadEGM96()

    #shift lateral values to the right reference and flatten coordinates
    lats = sp.deg2rad(-lat + 90).ravel()
    lons = sp.deg2rad(lon).ravel()

    #evaluate the spline and reshape the result
    evl = _EGM96.ev(lats, lons).reshape(lat.shape)

    return from_ndarray(itype, evl)
