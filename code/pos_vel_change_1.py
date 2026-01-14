from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from empm_1 import earth_velocity_itrf_from_mjd_astropy

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


t = Time(59000.0, format="mjd")

#print("Prędkość Ziemi:", earth_velocity_itrf_from_mjd_astropy(t.mjd))
#print("Norma prędkości Ziemi:", np.linalg.norm(earth_velocity_itrf_from_mjd_astropy(t.mjd)))


times = Time(59000.0, format="mjd") + np.arange(0, 365*0.2, 1) * u.day
# ↑ np. 5 lat, co 30 dni (więcej = gęściej)
print("Czas poczatkowy:", Time(59000.0, format="mjd"))

x, y, z = [], [], []
vx, vy, vz = [], [], []

for t in times:
    pos_gal, vel_gal = earth_velocity_itrf_from_mjd_astropy(t.mjd)

    x.append(pos_gal[0])
    y.append(pos_gal[1])
    z.append(pos_gal[2])

    vx.append(vel_gal[0])
    vy.append(vel_gal[1])
    vz.append(vel_gal[2])

    print("Norma prędkości [km/s]:", (vx[-1]**2 + vy[-1]**2 + vz[-1]**2)**0.5)


x = np.array(x)
y = np.array(y)
z = np.array(z)
vx = np.array(vx)
vy = np.array(vy)
vz = np.array(vz)

#print(x, y, z)
#print(vx, vy, vz)



fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection="3d")
ax.plot(x, y, z, lw=1, color="black", label="Tor punktu")

# współczynnik wizualny (DOBIERZ)
velocity_scale = 100  # au / (km/s)

ax.quiver(
    x, y, z,
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

ax.plot(x, y, z, lw=1)
ax.quiver(x, y, z, vx, vy, vz, length=0.002, normalize=True, color="red")
plt.show()