from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from earth_movement_clean import ICRS_to_Galactocentric


# ---- Zmiana położenia i prędkości z funkcji w earth_movement.py ----
def earth_vel_change(krok, mnoznik):

    #times = Time(59000.0, format="mjd") + np.arange(0, 365*1, 30) * u.day
    times = Time(59000.0, format="mjd") + np.arange(0, krok*mnoznik, krok) * u.day # ilość dni * ilość kroków, co tyle samo dni (krok)

    x, y, z = [], [], []
    vx, vy, vz = [], [], []
    x0, y0, z0 = [], [], []

    for t in times:
        pos_gal, vel_gal = ICRS_to_Galactocentric(t)
        pos_0, vel_0 = ICRS_to_Galactocentric(Time(59000.0, format="mjd"))

        x.append(pos_gal[0].to_value(u.au))
        y.append(pos_gal[1].to_value(u.au))
        z.append(pos_gal[2].to_value(u.au))

        vx.append(vel_gal[0].to_value(u.km/u.s))
        vy.append(vel_gal[1].to_value(u.km/u.s))
        vz.append(vel_gal[2].to_value(u.km/u.s))

        x0.append(pos_0[0].to_value(u.au))
        y0.append(pos_0[1].to_value(u.au))
        z0.append(pos_0[2].to_value(u.au))


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



# --- Zmiana pozycji wyliczona z prędkości ---
def pos_from_vel(krok, mnoznik):
    t = krok * 86400 # dni w sekundach (krok)
    
    x = [0]
    y = [0] 
    z = [0]
    tt = [0]

    for i in range(earth_vel_change(krok, mnoznik)[0].shape[0]):
        x.append(x[-1] + t*earth_vel_change(krok, mnoznik)[3][i])  # vx w km/s
        y.append(y[-1] + t*earth_vel_change(krok, mnoznik)[4][i])  # vy
        z.append(z[-1] + t*earth_vel_change(krok, mnoznik)[5][i])  # vz
        tt.append(tt[-1] + t)

        print(t*(i+1)/86400)


    return x, y, z, tt # w km



# ---- Zapis pozycji do pliku ----
def save_pos_to_file(krok, mnoznik):
    tt = np.array(pos_from_vel(krok, mnoznik)[3])/86400
    x, y, z = np.array(pos_from_vel(krok, mnoznik)[0]) / 1.496e+8, np.array(pos_from_vel(krok, mnoznik)[1]) / 1.496e+8, np.array(pos_from_vel(krok, mnoznik)[2]) / 1.496e+8    

    data = np.column_stack((tt, x, y, z))

    with open("figures/earth_movement/pos_75_10.txt", "w") as f:
        f.write("#t[dni]   x [AU]    y [AU]    z [AU]\n")
        np.savetxt(f, data, fmt="%.4e")

#save_pos_to_file(75,35)



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

print("========================================================================================")



# ---- Wykres pozycji z pliku ----
def plot_from_file(krok, mnoznik):
    data = np.loadtxt("figures/earth_movement/pos_25_40.txt", comments="#")
    #tt = data[:,0]
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]

    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(x, y, z, lw=1, color="green", label=f"Pozycja z pliku po {krok*mnoznik} dniach")

    ax.set_xlabel("X [au]")
    ax.set_ylabel("Y [au]")
    ax.set_zlabel("Z [au]")

    ax.legend()

    #set_axes_equal(ax)

    plt.show()


plot_from_file(25,40) #krok i mnożnik