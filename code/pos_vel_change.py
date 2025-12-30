from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from earth_movement_3 import ICRS_to_Galactocentric


def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    max_range = max([x_range, y_range, z_range]) / 2

    mid_x = np.mean(x_limits)
    mid_y = np.mean(y_limits)
    mid_z = np.mean(z_limits)

    ax.set_xlim3d(mid_x - max_range, mid_x + max_range)
    ax.set_ylim3d(mid_y - max_range, mid_y + max_range)
    ax.set_zlim3d(mid_z - max_range, mid_z + max_range)

#ax.plot(x, y, z)
#set_axes_equal(ax)



times = Time(59000.0, format="mjd") + np.arange(0, 365*2, 25) * u.day

x, y, z = [], [], []

for t in times:
    pos_gal, _ = ICRS_to_Galactocentric(t)
    x.append(pos_gal[0].to_value(u.au))
    y.append(pos_gal[1].to_value(u.au))
    z.append(pos_gal[2].to_value(u.au))

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot(x, y, z)
set_axes_equal(ax)
plt.show()


