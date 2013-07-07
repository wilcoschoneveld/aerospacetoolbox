import scipy as sp

#for linear relations:
#- float = unit to SI

#for non-linear relations:
#- list = [unit to SI, SI to unit]

#acceleration
_cacc = 1

#angle (done)
_cang = {'rad': 1.0,
         'deg': sp.pi/180,
         'rev': 2*sp.pi}

#angle velocity
_canv = 1

#angle acceleration
_cana = 1

#density
_cden = 1

#force
_cfor = 1

#length (done)
_clen = {'m': 1.0,
         'km': 1000.0,
         'mi': 5280 * 0.3048,
         'nmi': 1852.0,
         'in': 0.0254,
         'ft': 0.0254 * 12,
         'yd': 0.0254 * 36,
         'pt': 0.0254 / 72}

#mass
_cmas = {'kg': 1.0,
         'g': 1000.0,
         'lb': 0.45359237}

#pressure
_cpre = 1

#temperature (done)
_ctem = {'K': 1.0,
         'R': 5 / 9.0,
         'C': [lambda x: x + 273.15,
               lambda x: x - 273.15],
         'F': [lambda x: (x - 32)*(5/9.) + 273.15,
               lambda x: (x - 273.15)*(9/5.) + 32],
         'N': [lambda x: (x * 100/33.) + 273.15,
               lambda x: (x - 273.15) * 33. / 100]}

#velocity (done)
_cvel = {'m/s': 1.0,
         'ft/s': 0.3048,
         'km/s': 1000.0,
         'in/s': 0.0254,
         'km/h': 1 / 3.6,
         'mph': 5280 * 0.3048 / 3600,
         'kts': 1852.0 / 3600,
         'ft/min': 0.3048 / 60}

_cnvts = [_cacc, _cang, _canv, _cana, _cden,
          _cfor, _clen, _cmas, _cpre, _ctem, _cvel]

def convert(v, frm, to):
    for cnvt in _cnvts:
        if cnvt == 1:
            continue
        if frm in cnvt and to in cnvt:
            if isinstance(cnvt[frm], list):
                bv = cnvt[frm][0](v)
            else:
                bv = v * cnvt[frm]
            if isinstance(cnvt[to], list):
                return cnvt[to][1](bv)
            else:
                return bv / cnvt[to]
    raise ValueError('Could not convert between given units.')
