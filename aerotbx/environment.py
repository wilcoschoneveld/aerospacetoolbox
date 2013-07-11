import scipy as sp
from scipy.interpolate import RectSphereBivariateSpline
from pkg_resources import resource_string

_EGM96 = None

def _loadEGM96():
    #load the data resource file into a string
    flc = resource_string(__name__, 'egm96.dac')

    #setup basic coordinates
    lon = sp.linspace(0, 2*sp.pi, 1440, False)
    lat = sp.linspace(0, sp.pi, 721)

    #parse the raw data string
    data = sp.fromstring(flc, sp.dtype(sp.int16).newbyteorder('B'),
                 1038240, '').reshape((lat.size, lon.size)) / 100.0

    #interpolate the bad boy
    lut = RectSphereBivariateSpline(lat[1: -1], lon, data[1: -1],
               pole_values = (sp.mean(data[1]), sp.mean(data[-1])))

    return lut

def atmosisa(h, geopotential=True, T0=288.15, P0=101325.0):
    """
    Evaluate the international standard atmosphere (ISA) at a given altitude.
    The function assumes a continued troposphere below 0 meters and an infinite
    mesosphere above 84 kilometers geopotential height. atmosisa returns
    a tuple of temperature T, speed of sound A, pressure P, and a density RHO.

    call as:
        [T, a, P, rho] = atmosisa(h, geopotential, T0, P0)

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
    if not geopotential:
        h *= Re / (Re + h)

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

def geoidheight(lat, lon):
    """
    Calculate the geoid height using the EGM96 Geopotential Model.
    """

    global _EGM96

    lat = sp.array(lat, sp.float64, ndmin=1)
    lon = sp.array(lon, sp.float64, ndmin=1)

    if lat.shape != lon.shape:
        raise Exception("Inputs must contain equal number of values.")

    if (lat < -90).any() or (lat > 90).any() or not sp.isreal(lat).all():
        raise Exception("Lateral coordinates must be real numbers between -90 and 90 degrees.")

    if (lon < 0).any() or (lon > 360).any() or not sp.isreal(lon).all():
        raise Exception("Longitudinal coordinates must be real numbers between 0 and 360 degrees.")

    #if the model is not loaded, do so
    if _EGM96 is None: _EGM96 = _loadEGM96()

    #todo: handle array input

    #shift lateral values to the right reference and flatten coordinates
    lats = sp.deg2rad(-lat + 90).ravel()
    lons = sp.deg2rad(lon).ravel()

    #evaluate the spline and reshape the result
    evl = _EGM96.ev(lats, lons).reshape(lat.shape)

    #flatten solution if single value was given
    if evl.size == 1:
        evl = evl.flat[0]

    return evl
