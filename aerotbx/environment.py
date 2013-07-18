import scipy as sp
from scipy.interpolate import RectSphereBivariateSpline
from pkg_resources import resource_string
from aerotbx.utils import to_ndarray, from_ndarray

_EGM96 = None

def _loadEGM96():
    #load the data resource file into a string
    flc = resource_string(__name__, 'data/egm96.dac')

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

def stdatmos(**altitude):
    model = altitude.pop("model", None)

    if len(altitude) is not 1:
        raise Exception("Function needs exactly one altitude input.")

    mtype, alt = altitude.popitem()

    if mtype not in ["geom", "geop", "abs", "T", "P", "rho"]:
        raise Exception("The altitude input should be a valid input type.")
    
    itype, alt = to_ndarray(alt)

    R = 287.053 #gas constant [J/kg/K] (air)
    gamma = 1.4 #specific heat ratio [-] (air)
    
    g = 9.80665 #standard gravity [m/s^2] (earth)
    radius = 6356766.0 #earth radius [m] (earth)
    
    Tb = 288.15 #base temperature [K]
    Pb = 101325.0 #base pressure [Pa]

    Hb = [0, 11, 20, 32, 47, 51, 71, 84.852] #layer height [km]
    Lr = [-6.5, 0, 1, 2.8, 0, -2.8, -2] #lapse rate [K/km]

    if model is not None:
        raise Exception("Custom atmosphere model is incompatible.")

    Hb = sp.array(Hb, sp.float64) * 1000
    Lr = sp.array(Lr, sp.float64) / 1000

    Tb *= sp.ones(Lr.size + 1, sp.float64)
    Tb[1:] += sp.cumsum(Lr*(Hb[1:]-Hb[:-1]))

    T = sp.ones(alt.shape, sp.float64) * sp.nan
    P = sp.ones(alt.shape, sp.float64) * sp.nan

    if mtype is "geom": h = alt * radius / (radius + alt)
    elif mtype is "geop": h = alt
    elif mtype is "abs": h = alt - radius
    else: h = sp.ones(alt.shape, sp.float64) * sp.nan
        
    for i in xrange(len(Lr)):
        if mtype is "T":
            if not sp.isnan(h).any(): break
            if Lr[i] == 0:
                sel = (alt == Tb[i])
            else:
                s = sp.sign(Lr[i])
                bot = -sp.inf if i is 0 else Tb[i]*s
                top = sp.inf if i is len(Lr)-1 else Tb[i+1]*s
                sel = sp.logical_and(alt*s >= bot, alt*s < top)

            sel = sp.logical_and(sel, sp.isnan(h))

            T[sel] = alt[sel]
            h[sel] = Hb[i] + (alt[sel] - Tb[i])/Lr[i]
            
            if Lr[i] == 0:
                P[sel] = Pb * sp.exp((-g/(R*Tb[i]))*(h[sel]-Hb[i]))
            else:
                P[sel] = Pb * (T[sel] / Tb[i])**(-g/(Lr[i]*R))
                
        elif mtype in ["P", "rho"]:
            vb = Pb if mtype is "P" else Pb/(R*Tb[i])
            sel = alt <= (sp.inf if i is 0 else vb)
            if not sel.any(): break

            if Lr[i] == 0:
                T[sel] = Tb[i]
                h[sel] = Hb[i] - sp.log(alt[sel]/vb)*R*Tb[i]/g
            else:
                x = (Lr[i]*R + g) if mtype is "rho" else g
                T[sel] = Tb[i] * (alt[sel]/vb)**(-Lr[i]*R/x)
                h[sel] = Hb[i] + (T[sel] - Tb[i])/Lr[i]

            P[sel] = alt[sel] if mtype is "P" else alt[sel]*R*T[sel]
            
        else:
            sel = h >= (-sp.inf if i is 0 else Hb[i])
            if not sel.any(): break

            if Lr[i] == 0:
                T[sel] = Tb[i]
                P[sel] = Pb * sp.exp((-g/(R*Tb[i]))*(h[sel]-Hb[i]))
            else:
                T[sel] = Tb[i] + Lr[i] * (h[sel]-Hb[i])
                P[sel] = Pb * (T[sel] / Tb[i])**(-g/(Lr[i]*R))

        if Lr[i] == 0:
            Pb *= sp.exp((-g/(R*Tb[i]))*(Hb[i+1] - Hb[i]))
        else:
            Pb *= (Tb[i+1] / Tb[i])**(-g/(Lr[i]*R))

    h *= radius / (radius - h)
    rho = P / (R*T)
    a = sp.sqrt(gamma * R * T)
    
    return from_ndarray(itype, h, T, P, rho, a)

def geoidheight(lat, lon):
    """
    Calculate the geoid height using the EGM96 Geopotential Model.
    """

    global _EGM96

    #convert the input value to array
    itype, lat = to_ndarray(lat)
    itype, lon = to_ndarray(lon)

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

    return from_ndarray(itype, evl)
