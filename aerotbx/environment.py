import scipy as sp
import scipy.interpolate
from pkg_resources import resource_string

_EGM96 = None

def atmosisa(h, mtype="geom", T0=288.15, P0=101325.0):
    """
    Evaluate the international standard atmosphere (ISA) at a given altitude.
    The function assumes a continued troposphere below 0 meters and an infinite
    mesosphere above 84 kilometers geopotential height. atmosisa returns
    a tuple of temperature T, speed of sound A, pressure P, and a density RHO.

    call as:
        [T, a, P, rho] = atmosisa(h, mtype, T0, P0)

    Parameters
    ----------

        
    Returns
    -------

        
    Examples
    --------
    >>> atmosisa(-2000)
    (301.15409141737081, 347.88799859305686, 127782.84080627175,
        1.4781608008362668)

    References
    ----------
    Definition of the Standard Atmosphere:
        Anderson, John D. (2008). Introduction to Flight. (Sixth Edition
        International). NY: McGraw-Hill.
    """

    #convert the input value to array with ndmin=1
    h = sp.array(h, sp.float64, ndmin=1)

    #check if given input is valud
    if not sp.isreal(h).all():
        raise Exception("Height input must be real numbers.")

    #define constants
    R = 287.053
    g = 9.80665
    Re = 6356766.0

    #convert altitude to geopotential altitude if needed
    if mtype == "geop":
        pass
    elif mtype == "geom":
        h *= Re / (Re + h)
    else:
        raise Exception("mtype input must be an acceptable string to select second input parameter.")

    #define the international standard atmosphere
    Hb = sp.array([0, 11, 20, 32, 47, 51, 71, 84.852], sp.float64) * 1000
    Lr = sp.array([-6.5, 0, 1, 2.8, 0, -2.8, -2.0], sp.float64) * 0.001

    #define base conditions
    Tb = T0
    Pb = P0

    #create solution arrays
    T = sp.ones(h.shape, sp.float64) * Tb
    P = sp.ones(h.shape, sp.float64) * Pb

    #loop through the layers of the international standard atmosphere
    for i in xrange(7):
        #grab a selection with altitudes above current layer
        if i == 0: s = h > -sp.inf
        else: s = h > Hb[i]

        #if no altitudes are selected, stop looping
        if not s.any(): break

        #calculate the standard atmosphere from the hydrostatic equation
        if Lr[i] == 0:
            T[s] = Tb
            P[s] = Pb * sp.exp((-g/(R*Tb))*(h[s]-Hb[i]))

            #update new layer base values
            Pb *= sp.exp((-g/(R*Tb))*(Hb[i+1]-Hb[i]))
        else:
            T[s] = Tb + Lr[i]*(h[s]-Hb[i])
            P[s] = Pb * (T[s] / Tb)**(-g/(Lr[i]*R))

            #update new layer base values
            Tt = Tb + Lr[i] * (Hb[i+1]-Hb[i])
            Pb *= (Tt / Tb)**(-g/(Lr[i]*R))
            Tb = Tt

    #flatten solution if single value was given
    if h.size == 1:
        T = T.flat[0]
        P = P.flat[0]
    
    return T, sp.sqrt(1.4*R*T), P, P/(R*T)

def _loadEGM96():
    #load the data resource file into a string
    flc = resource_string(__name__, 'egm96.dac')

    #parse the raw data string
    data = sp.fromstring(flc, sp.dtype(sp.int16).newbyteorder('B')
                              ).reshape((721, 1440))[1:-1] / 100.0

    #setup basic coordinates
    lon = sp.linspace(0, 2*sp.pi, 1440, False)
    lat = sp.linspace(0, sp.pi, 721)[1: -1]

    #interpolate the bad boy
    lut = sp.interpolate.RectSphereBivariateSpline(lat, lon, data,
                                                   pole_values=(13.61, -29.53))

    return lut

def geoidheight(lat, lon):
    """
    Calculate the geoid height using the EGM96 Geopotential Model.
    """

    global _EGM96

    #if the model is not loaded, do so
    if _EGM96 is None: _EGM96 = _loadEGM96()

    #todo: handle array input
    return _EGM96.ev((-lat + 90.0)*sp.pi/180.0, lon*sp.pi/180.0)
