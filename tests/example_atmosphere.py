import scipy as sp
import matplotlib.pyplot as plt
from aerotbx import stdatmos

#create a linearly spaced height array
h = sp.linspace(-2, 82, 1000)

#use stdatmos to evaluate the standard atmosphere
_, T, P, rho, _ = stdatmos(geom=h*1000)

#use stdatmos to evaluate temperature altitude
Tt = sp.array([280, 250, 220, 200])
Th = stdatmos(T=Tt)[0] * 0.001

#use stdatmos to evaluate pressure altitude
Pp = sp.logspace(1, 4, 30)
Ph = stdatmos(P=Pp)[0] * 0.001

#use stdatmos to evaluate density altitude
rr = sp.logspace(-1, -4, 20)
rh = stdatmos(rho=rr)[0] * 0.001

#create a new plot
plt.figure(figsize=(13, 6))
plt.suptitle('International Standard Atmosphere')

#Temperature
plt.subplot(131)
plt.plot(T, h, 'r')
plt.plot(Tt, Th, 'kv')
plt.xlabel('Temperature [$K$]')
plt.ylabel('Height (km)')
plt.ylim(h[0], h[-1])

#Pressure
plt.subplot(132)
plt.semilogx(P, h, 'g')
plt.plot(Pp, Ph, 'k+')
plt.xlabel('Pressure [$Pa$]')
plt.ylim(h[0], h[-1])

#Density
plt.subplot(133)
plt.semilogx(rho, h, 'b')
plt.plot(rr, rh, 'k+')
plt.xlabel(r'Density [$N/m^2$]')
plt.ylim(h[0], h[-1])

plt.show()

