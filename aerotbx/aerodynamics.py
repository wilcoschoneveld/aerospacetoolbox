# pylint: disable=C0103,E1103
""" Aerospace Toolbox / aerodynamics.py

Functions for calculating aerodynamic properties of different flows.

"""

import scipy as sp
from aerotbx.utils import to_ndarray, from_ndarray

_AETB_iternum = 10

def _flowinput(flow):
    """Parse the flow input used in the aerodynamics module"""
    
    #pop flow from the input
    gamma = flow.pop("gamma", 1.4)

    #check if single input is given
    if len(flow) != 1:
        raise Exception("Function needs exactly one flow input.")

    #pop the flow variable and type
    mtype, flow = flow.popitem()

    if not sp.isreal(gamma).all() or not sp.isreal(flow).all():
        raise Exception("Flow input variables must be real numbers.")

    #convert the input values to arrays
    ftype, flow = to_ndarray(flow)
    gtype, gamma = to_ndarray(gamma)

    #check if the given gamma value is valid
    if (gamma <= 1).any():
        raise Exception("Specific heat ratio inputs must be real numbers" \
            " greater than 1.")

    #if both inputs are non-scalar, they should be equal in shape
    if gamma.size > 1 and flow.size > 1 and gamma.shape != flow.shape:
        raise Exception("Inputs must be same shape or at least one input" \
            " must be scalar.")

    #if one of the variables is an array, the other should match it size
    if gamma.size > flow.size:
        n = gamma.shape
        itype = gtype
        flow = sp.ones(n, sp.float64) * flow[0]
    else:
        n = flow.shape
        itype = ftype
        gamma = sp.ones(n, sp.float64) * gamma[0]

    return gamma, flow, mtype.lower(), itype

def flowisentropic(**flow):
    """
    Evaluate the isentropic relations with a given set of specific heat
    ratios and any one of the isentropic flow variables. flowisentropic
    returns a tuple of isentropic flow Mach number M, temperature ratio T,
    pressure ratio P, density ratio RHO, and area ratio AREA.

    call as:
        M, T, P, rho, area = flowisentropic(**flow)
    """

    #parse the input
    gamma, flow, mtype, itype = _flowinput(flow)

    #calculate gamma-ratios for use in the equations
    a = (gamma+1) / 2
    b = (gamma-1) / 2
    c = a / (gamma-1)

    #preshape mach array
    M = sp.empty(flow.shape, sp.float64)

    #use the isentropic relations to solve for the mach number
    if mtype in ["mach", "m"]:
        if (flow < 0).any() or not sp.isreal(flow).all():
            raise Exception("Mach number inputs must be real numbers" \
                " greater than or equal to 0.")
        M = flow
    elif mtype in ["temp", "t"]:
        if (flow < 0).any() or (flow > 1).any():
            raise Exception("Temperature ratio inputs must be real numbers" \
                " 0 <= T <= 1.")
        M[flow == 0] = sp.inf
        M[flow != 0] = sp.sqrt((1/b[flow != 0])*(flow[flow != 0]**(-1) - 1))
    elif mtype in ["pres", "p"]:
        if (flow < 0).any() or (flow > 1).any():
            raise Exception("Pressure ratio inputs must be real numbers" \
                " 0 <= P <= 1.")
        M[flow == 0] = sp.inf
        M[flow != 0] = sp.sqrt((1/b[flow != 0]) * \
            (flow[flow != 0]**((gamma[flow != 0]-1)/-gamma[flow != 0]) - 1))
    elif mtype in ["dens", "d", "rho"]:
        if (flow < 0).any() or (flow > 1).any():
            raise Exception("Density ratio inputs must be real numbers" \
                " 0 <= rho <= 1.")
        M[flow == 0] = sp.inf
        M[flow != 0] = sp.sqrt((1/b[flow != 0]) * \
            (flow[flow != 0]**((gamma[flow != 0]-1)/-1) - 1))
    elif mtype in ["sub", "sup"]:
        if (flow < 1).any():
            raise Exception("Area ratio inputs must be real numbers greater" \
                " than or equal to 1.")
        M[:] = 0.2 if mtype == "sub" else 1.8
        for _ in xrange(_AETB_iternum):
            K = M ** 2
            f = -flow + a**(-c) * ((1+b*K)**c) / M #mach-area relation
            g = a**(-c) * ((1+b*K)**(c-1)) * (b*(2*c - 1)*K - 1) / K #deriv
            M = M - (f / g) #Newton-Raphson
        M[flow == 1] = 1
        M[sp.isinf(flow)] = sp.inf
    else:
        raise Exception("Keyword input must be an acceptable string to" \
            " select input parameter.")

    d = 1 + b*M**2
    
    T = d**(-1)
    P = d**(-gamma/(gamma-1))
    rho = d**(-1/(gamma-1))

    area = sp.empty(M.shape, sp.float64)
    r = sp.logical_and(M != 0, sp.isfinite(M))
    area[r] = a[r]**(-c[r]) * d[r]**c[r] / M[r]
    area[sp.logical_not(r)] = sp.inf

    return from_ndarray(itype, M, T, P, rho, area)

def flownormalshock(**flow):
    """Normal shock relations"""

    #parse the input
    gamma, flow, mtype, itype = _flowinput(flow)

    #calculate gamma-ratios for use in the equations
    a = (gamma+1) / 2
    b = (gamma-1) / 2
    c = gamma / (gamma-1)

    #preshape mach array
    M = sp.empty(flow.shape, sp.float64)

    #use the normal shock relations to solve for the mach number
    if mtype in ["mach", "m1", "m"]:
        if (flow < 1).any():
            raise Exception("Mach number inputs must be real numbers" \
                " greater than or equal to 1.")
        M = flow
    elif mtype in ["down", "mach2", "m2", "md"]:
        lowerbound = sp.sqrt((gamma - 1)/(2*gamma))
        if (flow < lowerbound).any() or (flow > 1).any():
            raise Exception("Mach number downstream inputs must be real" \
                " numbers SQRT((GAMMA-1)/(2*GAMMA)) <= M <= 1.")
        M[flow <= lowerbound] = sp.inf
        M[flow > lowerbound] = sp.sqrt((1 + b*flow**2) / (gamma*flow**2 - b))
    elif mtype in ["pres", "p"]:
        if (flow < 1).any() or not sp.isreal(flow).all():
            raise Exception("Pressure ratio inputs must be real numbers" \
                " greater than or equal to 1.")
        M = sp.sqrt(((flow-1)*(gamma+1)/(2*gamma)) + 1)
    elif mtype in ["dens", "d", "rho"]:
        upperbound = (gamma+1) / (gamma-1)
        if (flow < 1).any() or (flow > upperbound).any():
            raise Exception("Density ratio inputs must be real numbers" \
                " 1 <= M <= (GAMMA+1)/(GAMMA-1).")
        M[flow >= upperbound] = sp.inf
        M[flow < upperbound] = sp.sqrt(2*flow / (1 + gamma + flow - flow*gamma))
    elif mtype in ["temp", "t"]:
        if (flow < 1).any():
            raise Exception("Temperature ratio inputs must be real numbers" \
                " greater than or equal to 1.")
        B = b + gamma/a - gamma*b/a - flow*a
        M = sp.sqrt((-B + sp.sqrt(B**2 - \
            4*b*gamma*(1-gamma/a)/a)) / (2*gamma*b/a))
    elif mtype in ["totalp", "p0"]:
        if (flow < 0).any() or (flow > 1).any():
            raise Exception("Total pressure ratio inputs must be real" \
                " numbers 0 <= P0 <= 1.")
        M[:] = 2.0 #initial guess for the solution
        for _ in xrange(_AETB_iternum):
            f = -flow + (1 + (gamma/a)*(M**2 - 1))**(1-c) \
                * (a*M**2 / (1 + b*M**2))**c
            g = 2*M*(a*M**2 / (b*M**2 + 1))**(c-1)*((gamma* \
                (M**2-1))/a + 1)**(-c)*(a*c + gamma*(-c*(b*M**4 + 1) \
                + b*M**4 + M**2)) / (b*M**2 + 1)**2
            M = M - (f / g) #Newton-Raphson
        M[flow == 0] = sp.inf
    elif mtype in ["pito", "pitot", "rpr"]:
        lowerbound = a**c
        if (flow < lowerbound).any():
            raise Exception("Rayleigh-Pitot ratio inputs must be real" \
                " numbers greater than or equal to ((G+1)/2)**(-G/(G+1)).")
        M[:] = 5.0 #initial guess for the solution
        K = a**(2*c - 1)
        for _ in xrange(_AETB_iternum):
            f = -flow + K * M**(2*c) / (gamma*M**2 - b)**(c-1) #Rayleigh-Pitot
            g = 2*K*M**(2*c - 1) * (gamma*M**2 - b)**(-c) * (gamma*M**2 - b*c)
            M = M - (f / g) #Newton-Raphson
    else:
        raise Exception("Keyword input must be an acceptable string to" \
            " select input parameter.")

    #normal shock relations
    M2 = sp.sqrt((1 + b*M**2) / (gamma*M**2 - b))
    rho = ((gamma+1)*M**2) / (2 + (gamma-1)*M**2)
    P = 1 + (M**2 - 1)*gamma / a
    T = P / rho
    P0 = P**(1-c) * rho**(c)
    P1 = a**(2*c - 1) * M**(2*c) / (gamma*M**2 - b)**(c-1)

    #handle infinite mach
    M2[M == sp.inf] = sp.sqrt((gamma[M == sp.inf] - 1)/(2*gamma[M == sp.inf]))
    T[M == sp.inf] = sp.inf
    P[M == sp.inf] = sp.inf
    rho[M == sp.inf] = (gamma[M == sp.inf]+1) / (gamma[M == sp.inf]-1)
    P0[M == sp.inf] = 0
    P1[M == sp.inf] = sp.inf    

    return from_ndarray(itype, M, M2, T, P, rho, P0, P1)

def flowprandtlmeyer(**flow):
    """Prandtl-Meyer function for expansion waves"""

    #parse the input
    gamma, flow, mtype, itype = _flowinput(flow)

    #calculate gamma-ratios for use in the equations
    l = sp.sqrt((gamma-1)/(gamma+1))

    #preshape mach array
    M = sp.empty(flow.shape, sp.float64)

    #use prandtl-meyer relation to solve for the mach number
    if mtype in ["mach", "m"]:
        if (flow < 1).any():
            raise Exception("Mach number inputs must be real numbers greater" \
                " than or equal to 1.")
        M = flow
    elif mtype in ["mu", "machangle"]:
        if (flow < 0).any() or  (flow > 90).any():
            raise Exception("Mach angle inputs must be real numbers" \
                " 0 <= M <= 90.")
        M = 1 / sp.sin(sp.deg2rad(flow))
    elif mtype in ["nu", "pm", "pmangle"]:
        if (flow < 0).any() or  (flow > 90*((1/l)-1)).any():
            raise Exception("Prandtl-Meyer angle inputs must be real" \
                " numbers 0 <= M <= 90*(sqrt((g+1)/(g-1))-1).")
        M[:] = 2 #initial guess for the solution
        for _ in xrange(_AETB_iternum):
            b = sp.sqrt(M**2 - 1)
            f = -sp.deg2rad(flow) + (1/l) * sp.arctan(l*b) - sp.arctan(b)
            g = b*(1 - l**2) / (M*(1 + (l**2)*(b**2))) #derivative
            M = M - (f / g) #Newton-Raphson
    else:
        raise Exception("Keyword input must be an acceptable string to" \
            " select input parameter.")

    #normal shock relations
    b = sp.sqrt(M**2 - 1)
    V = (1/l) * sp.arctan(l*b) - sp.arctan(b)
    U = sp.arcsin(1 / M)

    return from_ndarray(itype, M, sp.rad2deg(V), sp.rad2deg(U))
