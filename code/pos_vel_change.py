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


def earth_vel_change():

    times = Time(59000.0, format="mjd") + np.arange(0, 365*1, 30) * u.day
    # ↑ np. 5 lat, co 30 dni (więcej = gęściej)
    #print("Czas poczatkowy:", Time(59000.0, format="mjd"))

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

        #print("Norma prędkości [km/s]:", (vx[-1]**2 + vy[-1]**2 + vz[-1]**2)**0.5)
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

    return x, y, z, vx, vy, vz, x0, y0, z0


#rint((earth_vel_change()[3]**2 + earth_vel_change()[4]**2 + earth_vel_change()[5]**2)**0.5)





def plot_earth_velocity_change():
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(earth_vel_change()[0], earth_vel_change()[1], earth_vel_change()[2], lw=1, color="black", label="Tor punktu")
    #ax.plot(x_plot, y_plot, z_plot, lw=1, color="black", label="Tor punktu")

    # współczynnik wizualny (DOBIERZ)
    velocity_scale = 0.1  # au / (km/s)

    ax.quiver(
        earth_vel_change()[0], earth_vel_change()[1], earth_vel_change()[2],          # początek strzałek
        earth_vel_change()[3], earth_vel_change()[4], earth_vel_change()[5],       # kierunek (prędkość)
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

    ax.plot(earth_vel_change()[0], earth_vel_change()[1], earth_vel_change()[2], lw=1)
    ax.quiver(earth_vel_change()[0], earth_vel_change()[1], earth_vel_change()[2],
            earth_vel_change()[3], earth_vel_change()[4], earth_vel_change()[5],
            length=0.002, normalize=True, color="red")
    #ax.plot(x_plot, y_plot, z_plot, lw=1)
    #ax.quiver(x_plot, y_plot, z_plot, vx, vy, vz, length=0.002, normalize=True, color="red")
    plt.show()


#lot_earth_velocity_change()


def pos_from_vel():
    t = 30 * u.day
    #t = np.arange(0, 365*1, 30) * u.day

    x = 0
    y = 0
    z = 0

    x = x + t*earth_vel_change()[3]  # vx
    y = y + t*earth_vel_change()[4]  # vy
    z = z + t*earth_vel_change()[5]  # vz

    return x, y, z


print(pos_from_vel())

def plot_pos_from_vel():
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(pos_from_vel()[0], pos_from_vel()[1], pos_from_vel()[2], lw=1, color="blue", label="Pozycja z prędkości")

    ax.set_xlabel("X [au]")
    ax.set_ylabel("Y [au]")
    ax.set_zlabel("Z [au]")

    ax.set_title("Pozycja punktu na Ziemi wyliczona z prędkości w układzie Galactocentric")
    ax.legend()

    set_axes_equal(ax)

    plt.show()

plot_pos_from_vel()