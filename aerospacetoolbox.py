"""
aerospacetoolbox.py

Functions for aerospace analysis to develop and evaluate your designs.
Currently supports the following functions: flowisentropic, atmosisa

Created to look like the equivalent set of functions in the matlab
aerospace toolbox: http://www.mathworks.nl/help/aerotbx/index.html

Author: Wilco Schoneveld
Date: 7 June 2013
Version: 0.2
"""

import scipy as sp

_AETB_iternum = 10

def flowisentropic(gamma, flow, mtype="mach"):
    """
    Evaluate the isentropic relations with a given set of specific heat ratios
    and any one of the isentropic flow variables. flowisentropic returns a
    tuple of isentropic flow Mach number M, temperature ratio T,
    pressure ratio P, density ratio RHO, and area ratio AREA.

    call as:
        [M, T, P, rho, area] = flowisentropic(gamma, flow, mtype)

    Parameters
    ----------
    gamma : ndarray or scalar
        Array of specific heat ratios. gamma must be a scalar or array of real
        numbers greater than 1.
    flow : ndarray or scalar
        Array of real numerical values for one of the isentropic flow
        relations. This argument can be one of the following:
        - Mach numbers: flow must be a scalar or an array of real numbers
        greater than or equal to 0.
        - Temperature ratios: flow must be a scalar or an array of real numbers
        between and including 0 and 1.
        - Pressure ratios: flow must be a scalar or an array of real numbers
        between and including 0 and 1.
        - Density ratios: flow must be a scalar or an array of real numbers
        between and including 0 and 1.
        - Area ratios: flow must be a scalar or an array of real numbers
        greater than or equal to 1.
        * If flow and gamma are ndarrays, they must be the same size.
    mtype : {'mach', 'temp', 'pres',  'dens', 'sub', 'sup'}, optional
        A string that defines the input mode for the isentropic flow
        in flow. This argument can be one of the following:
        - 'mach' or 'm': Mach number (Default)
        - 'temp' or 't': Temperature ratio.
        - 'pres' or 'p': Pressure ratio.
        - 'dens' or 'd' or 'rho': Density ratio.
        - 'sub': Area ratio. Subsonic solution of the mach area relation.
        - 'sup': Area ratio. Supersonic solution of the mach area relation.
        
    Returns
    -------
    M : ndarray or scalar
        An array of Mach numbers.
    T : ndarray or scalar
        An array of temperature ratios. The temperature ratio is
        defined as the local static temperature over the stagnation
        temperature.
    P : ndarray or scalar
        An array of pressure ratios. The pressure ratio is defined
        as the local static pressure over the stagnation pressure.
    rho : ndarray or scalar
        An array of density ratios. The density ratio is defined as
        the local density over the stagnation density.
    area : ndarray or scalar
        An array of area ratios. The area ratio is defined as the
        local streamtube area over the reference streamtube area for
        sonic conditions.
        
    Examples
    --------
    >>> flowisentropic(1.4, 1.1, 'sub')
    (0.69239913946385423, 0.91250590964351963, 0.72581326332969454,
        0.79540664412052009, 1.0999999999999999)
    >>> [M, T, P, rho, area] = flowisentropic(1.4, datalist, 'sub')
    [array([...]),array([...]),array([...]),array([...]),array([...])]

    References
    ----------
    Isentropic flow equations:
        http://www.grc.nasa.gov/WWW/k-12/airplane/isentrop.html
    """

    #convert the input values to arrays with a minimal dimension of 1
    gamma = sp.array(gamma, sp.float64, ndmin=1)
    flow = sp.array(flow, sp.float64, ndmin=1)

    #check if the given gamma value is valid
    if (gamma <= 1).any() or not sp.isreal(gamma).all():
        raise Exception("Specific heat ratio inputs must be real numbers greater than 1.")

    #if both inputs are non-scalar, they should be equal in shape
    if gamma.size > 1 and flow.size > 1 and gamma.shape != flow.shape:
        raise Exception("Inputs must be same shape or at least one input must be scalar.")

    #if one of the variables is an array, the other should match it size
    n = gamma.shape if gamma.size > flow.size else flow.shape
    if flow.size == 1: flow = sp.ones(n, sp.float64) * flow.flat[0]
    if gamma.size == 1: gamma = sp.ones(n, sp.float64) * gamma.flat[0]

    #calculate gamma-ratios for use in the equations
    a = (gamma+1) / 2
    b = (gamma-1) / 2
    c = a / (gamma-1)

    #check what the input type is, and use the isentropic relations to solve for the mach number
    if mtype == "mach" or mtype == "m":
        if (flow < 0).any() or not sp.isreal(flow).all():
            raise Exception("Mach number inputs must be real numbers greater than 1.")
        M = flow
    elif mtype == "temp" or mtype == "t":
        if (flow < 0).any() or (flow > 1).any() or not sp.isreal(flow).all():
            raise Exception("Temperature ratio inputs must be real numbers 0 <= T <= 1.")
        M = sp.empty(n, sp.float64)
        M[flow == 0] = sp.inf
        M[flow != 0] = sp.sqrt((1/b[flow != 0])*(flow[flow != 0]**(-1) - 1))
    elif mtype == "pres" or mtype == "p":
        if (flow < 0).any() or (flow > 1).any() or not sp.isreal(flow).all():
            raise Exception("Pressure ratio inputs must be real numbers 0 <= P <= 1.")
        M = sp.empty(n, sp.float64)
        M[flow == 0] = sp.inf
        M[flow != 0] = sp.sqrt((1/b[flow != 0])*(flow[flow != 0]**((gamma[flow != 0]-1)/-gamma[flow != 0]) - 1))
    elif mtype == "dens" or mtype == "d" or mtype == "rho":
        if (flow < 0).any() or (flow > 1).any() or not sp.isreal(flow).all():
            raise Exception("Density ratio inputs must be real numbers 0 <= rho <= 1.")
        M = sp.empty(n, sp.float64)
        M[flow == 0] = sp.inf
        M[flow != 0] = sp.sqrt((1/b[flow != 0])*(flow[flow != 0]**((gamma[flow != 0]-1)/-1) - 1))
    elif mtype == "sub" or mtype == "sup":
        if (flow < 1).any() or not sp.isreal(flow).all():
            raise Exception("Area ratio inputs must be real numbers greater than or equal to 1.")
        M = sp.empty(n, sp.float64)
        if mtype == "sub": M[:] = 0.2 #initial guess for the subsonic solution
        if mtype == "sup": M[:] = 1.8 #initial guess for the supersonic solution
        for i in xrange(_AETB_iternum):
            K = M ** 2
            f = -flow + a**(-c) * ((1+b*K)**c) / M #mach-area relation
            g = a**(-c) * ((1+b*K)**(c-1)) * (b*(2*c - 1)*K - 1) / K #derivative
            M = M - (f / g) #Newton-Raphson
        M[flow == 1] = 1
    else:
        raise Exception("Third input must be an acceptable string to select second input parameter.")

    #if single mach number is calculated
    if M.size == 1:
        #flatten the values to a scalar
        M = M.flat[0]; gamma = gamma.flat[0]

        #recalculate gamma-ratios as scalar values
        a = (gamma+1) / 2
        b = (gamma-1) / 2
        c = a / (gamma-1)

        #insert values in the isentropic relations
        d = 1 + b*M**2
        T = d**(-1)
        P = d**(-gamma/(gamma-1))
        rho = d**(-1/(gamma-1))

        #the mach-area relation has limits 0 and infinite
        if sp.isinf(M): area = 0
        elif M == 0: area = sp.inf
        else: area = a**(-c) * d**c / M

        return M, T, P, rho, area

    #if an ndarray M is calculated
    else:
        #insert values in the isentropic relations
        d = 1 + b*M**2
        T = d**(-1)
        P = d**(-gamma/(gamma-1))
        rho = d**(-1/(gamma-1))

        #start with the mach-area limits
        area = sp.zeros(n, sp.float64)
        area[M == 0] = sp.inf

        #calculate the mach-area relation only when non-zero and non-infinite
        r = sp.logical_and(sp.logical_not(M==0),sp.logical_not(sp.isinf(M))) 
        area[r] = a[r]**(-c[r]) * d[r]**c[r] / M[r]

        return M, T, P, rho, area

def atmosisa(h, mtype="", Toffset=0):
    """
    Evaluate the international standard atmosphere (ISA) at a given altitude.
    The function assumes a continued troposphere below 0 meters and an infinite
    mesosphere above 71 kilometers geopotential height. atmosisa returns
    a tuple of temperature T, speed of sound A, pressure P, and a density RHO.

    call as:
        [T, a, P, rho] = atmosisa(h, mtype, Toffset)

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
    Pb = 101325.0
    R = 287.053
    g = 9.80665
    Re = 6356766.0

    #convert altitude to geopotential altitude if needed
    if not mtype == "geop":
        h *= Re / (Re + h)

    #define the international standard atmosphere
    Hb = sp.array([0, 11, 20, 32, 47, 51, 71, 84.852], sp.float64) * 1000
    Lr = sp.array([-6.5, 0, 1, 2.8, 0, -2.8, -2.0], sp.float64) * 0.001
    Tb = sp.array([15, -56.5, -56.5, -44.5, -2.5, -2.5, -58.5, -86.28])
    Tb += 273.15 + Toffset

    #create solution arrays
    T = sp.ones(h.shape, sp.float64) * Tb[0]
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
            T[s] = Tb[i]
            P[s] = Pb * sp.exp((-g/(R*Tb[i]))*(h[s]-Hb[i]))
            Pb *= sp.exp((-g/(R*Tb[i]))*(Hb[i+1]-Hb[i]))
        else:
            T[s] = Tb[i] + Lr[i]*(h[s]-Hb[i])
            P[s] = Pb * (T[s] / Tb[i])**(-g/(Lr[i]*R))
            Pb *= (Tb[i+1] / Tb[i])**(-g/(Lr[i]*R))

    #flatten solution if single value was given
    if h.size == 1:
        T = T.flat[0]
        P = P.flat[0]
    
    return T, sp.sqrt(1.4*R*T), P, P/(R*T)
