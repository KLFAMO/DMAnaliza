from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from earth_movement_3 import ICRS_to_Galactocentric
from empm import earth_velocity_itrf_from_mjd_astropy


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




times = Time(59000.0, format="mjd") + np.arange(0, 365*40, 100) * u.day
# ↑ np. 5 lat, co 30 dni (więcej = gęściej)

x, y, z = [], [], []
vx, vy, vz = [], [], []
x0, y0, z0 = [], [], []

for t in times:
    pos_gal, vel_gal = ICRS_to_Galactocentric(t)
    pos_0, vel_0 = ICRS_to_Galactocentric(Time(59000.0, format="mjd"))
    #pos_gal, vel_gal = earth_velocity_itrf_from_mjd_astropy(t)

    x.append(pos_gal[0].to_value(u.au))
    y.append(pos_gal[1].to_value(u.au))
    z.append(pos_gal[2].to_value(u.au))

    vx.append(vel_gal[0].to_value(u.km/u.s))
    vy.append(vel_gal[1].to_value(u.km/u.s))
    vz.append(vel_gal[2].to_value(u.km/u.s))

    x0.append(pos_0[0].to_value(u.au))
    y0.append(pos_0[1].to_value(u.au))
    z0.append(pos_0[2].to_value(u.au))

    #x_plot, y_plot, z_plot = [], [], []

    #x_plot = x - x0
    #y_plot = y - y0
    #z_plot = z - z0



x = np.array(x)
y = np.array(y)
z = np.array(z)
vx = np.array(vx)
vy = np.array(vy)
vz = np.array(vz)
x0 = np.array(x0)
y0 = np.array(y0)
z0 = np.array(z0)

x_plot = x - x0
y_plot = y - y0
z_plot = z - z0

print()
print(x_plot)
print()
print(y_plot)
print()
print(z_plot)




fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection="3d")

ax.plot(x_plot, y_plot, z_plot, lw=1, color="black", label="Tor punktu")

# współczynnik wizualny (DOBIERZ)
velocity_scale = 0.1  # au / (km/s)


ax.quiver(
    x_plot, y_plot, z_plot,          # początek strzałek
    vx, vy, vz,       # kierunek (prędkość)
    length=velocity_scale,
    normalize=True,
    color="red",
    linewidth=1,
    label="Kierunek prędkości"
)


ax.set_xlabel("X [au]")
ax.set_ylabel("Y [au]")
ax.set_zlabel("Z [au]")

ax.set_title("Ruch punktu na Ziemi w układzie Galactocentric")
ax.legend()

set_axes_equal(ax)


'''x0, y0, z0 = x[0], y[0], z[0]

x_plot = x - x0
y_plot = y - y0
z_plot = z - z0'''


#ax.plot(x, y, z, lw=1)
#ax.quiver(x, y, z, vx, vy, vz,
#          length=0.002, normalize=True, color="red")
ax.plot(x_plot, y_plot, z_plot, lw=1)
ax.quiver(x_plot, y_plot, z_plot, vx, vy, vz,
          length=0.002, normalize=True, color="red")
plt.show()



