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
    #pop atmospherical model from input
    model = altitude.pop("model", None)

    #check if a single remaining input exists
    if len(altitude) is not 1:
        raise Exception("Function needs exactly one altitude input.")

    #pop the altitude input
    mtype, alt = altitude.popitem()

    #check if the altitude input type is valid
    if mtype not in ["geom", "geop", "abs", "T", "P", "rho"]:
        raise Exception("The altitude input should be a valid input type.")

    #convert the input to numpy arrays
    itype, alt = to_ndarray(alt)

    #define default atmospheric values
    R = 287.053 #gas constant [J/kg/K] (air)
    gamma = 1.4 #specific heat ratio [-] (air)
    
    g = 9.80665 #standard gravity [m/s^2] (earth)
    radius = 6356766.0 #earth radius [m] (earth)
    
    Tb = 288.15 #base temperature [K]
    Pb = 101325.0 #base pressure [Pa]

    Hb = [0, 11, 20, 32, 47, 51, 71, 84.852] #layer height [km]
    Lr = [-6.5, 0, 1, 2.8, 0, -2.8, -2] #lapse rate [K/km]

    #TODO: add custom model compatibility (overwrite default values)
    if model is not None:
        raise Exception("Custom atmosphere model is incompatible.")

    #convert height and base layer to correct units
    Hb = sp.array(Hb, sp.float64) * 1000
    Lr = sp.array(Lr, sp.float64) / 1000

    #preshape solution arrays
    T = sp.ones(alt.shape, sp.float64) * sp.nan
    P = sp.ones(alt.shape, sp.float64) * sp.nan

    #define the height array
    if mtype is "geom": h = alt * radius / (radius + alt)
    elif mtype is "geop": h = alt
    elif mtype is "abs": h = alt - radius
    else: h = sp.ones(alt.shape, sp.float64) * sp.nan

    for lr, hb, ht in zip(Lr, Hb[:-1], Hb[1:]):
        #calculate the temperature at layer top
        Tt = Tb + lr*(ht-hb)
        
        if mtype is "T":
            #break the loop if there are no nans in the solution array
            if not sp.isnan(h).any(): break

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
            if not sel.any(): break

            #solve for temperature and height
            if lr == 0:
                T[sel] = Tb
                h[sel] = hb - sp.log(alt[sel]/vb)*R*Tb/g
            else:
                x = g if mtype is "P" else (lr*R + g)
                T[sel] = Tb * (alt[sel]/vb)**(-lr*R/x)
                h[sel] = hb + (T[sel] - Tb) / lr

            #pressure is given as input (or as density, easy convert to pressure)
            P[sel] = alt[sel] if mtype is "P" else alt[sel]*R*T[sel]

        else:
            #select all height values above layer base
            sel = h >= (-sp.inf if hb == 0 else hb)

            #break if nothing is selected
            if not sel.any(): break

            #solve for temperature and pressure
            if lr == 0:
                T[sel] = Tb
                P[sel] = Pb * sp.exp((-g/(R*Tb))*(h[sel]-hb))
            else:
                T[sel] = Tb + lr * (h[sel] - hb)
                P[sel] = Pb * (T[sel] / Tb)**(-g/(lr*R))

        #update pressure base value
        if lr == 0: Pb *= sp.exp((-g/(R*Tb))*(ht - hb))
        else: Pb *= (Tt / Tb)**(-g/(lr*R))

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

    #shift lateral values to the right reference and flatten coordinates
    lats = sp.deg2rad(-lat + 90).ravel()
    lons = sp.deg2rad(lon).ravel()

    #evaluate the spline and reshape the result
    evl = _EGM96.ev(lats, lons).reshape(lat.shape)

    return from_ndarray(itype, evl)
