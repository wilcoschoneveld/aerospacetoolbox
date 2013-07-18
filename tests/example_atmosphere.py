import scipy as sp
import matplotlib.pyplot as plt
from aerotbx import stdatmos

#create a linearly spaced height array
h = sp.linspace(-1000, 80000, 1000)

#use stdatmos to evaluate the standard atmosphere
h, T, P, rho, a = stdatmos(geom=h)

#create a logspaced pressure array (to simulate pressure measurements)
pstatic = sp.logspace(1, 4, 30)

#use stdatmos to obtain geometrical height from pressures
h2 = stdatmos(P=pstatic)[0]

#plot the temperature gradient
plt.figure()
plt.plot(T, h)
plt.title('temperature')

#plot the pressure curve and simulated measurements
plt.figure()
plt.semilogx(P, h)
plt.plot(pstatic, h2, '+')
plt.title('pressure')

#plot the density curve
plt.figure()
plt.semilogx(rho, h)
plt.title('density')

plt.show()

