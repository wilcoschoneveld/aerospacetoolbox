import scipy as sp
import matplotlib.pyplot as plt
from aerotbx import geoidheight

lats = sp.linspace(90, -90, 180)
lons = sp.linspace(-180, 180, 360) % 360

lons, lats = sp.meshgrid(lons, lats)

h = geoidheight(lats, lons)

plt.imshow(h, extent=[-180, 180, -90, 90])
plt.colorbar()
plt.show()
