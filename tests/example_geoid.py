import scipy as sp
import matplotlib.pyplot as plt
from aerotbx import geoidheight

lats = sp.linspace(90, -90, 180)
lons = sp.linspace(0, 360, 360)

lons, lats = sp.meshgrid(lons, lats)

h = geoidheight(lats, lons)

plt.imshow(h, extent=[0, 360, -90, 90])
plt.show()
