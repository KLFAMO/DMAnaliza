from astropy.time import Time
import astropy.units as u
from astropy.coordinates import CartesianRepresentation
import numpy as np
import matplotlib.pyplot as plt
from earth_movement import ITRS_to_Galactocentric

xyz_list=[6378137, 0, 0]
point_on_earth = CartesianRepresentation(*xyz_list) * u.m

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


# ---- Zmiana położenia i prędkości z funkcji w earth_movement_3.py ----
def earth_vel_change(krok, mnoznik):
    
    #times = Time(59000.0, format="mjd") + np.arange(0, 365*1, 30) * u.day
    #times = Time(59000.0, format="mjd") + np.arange(0, 10*30, 10) * u.day # ilość dni * ilość kroków, co tyle samo dni (krok)
    times = Time(59000.0, format="mjd") + np.arange(0, krok*mnoznik, krok) * u.day # ilość dni * ilość kroków, co tyle samo dni (krok)
    #print("Czas poczatkowy:", Time(59000.0, format="mjd"))

    x, y, z = [], [], []
    vx, vy, vz = [], [], []
    x0, y0, z0 = [], [], []

    for t in times:
        pos_gal, vel_gal = ITRS_to_Galactocentric(t, point_on_earth)
        pos_0, vel_0 = ITRS_to_Galactocentric(Time(59000.0, format="mjd"), point_on_earth)
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

#print((earth_vel_change(10,20)[3]**2 + earth_vel_change(10,20)[4]**2 + earth_vel_change(10,20)[5]**2)**0.5)




# ---- Wykres zmiany położenia z funkcji earth_vel_change ----
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

#pos_from_vel(10,10)



# ---- Zapis pozycji do pliku ----
def save_pos_to_file(krok, mnoznik):
    tt = np.array(pos_from_vel(krok, mnoznik)[3])/86400
    x, y, z = np.array(pos_from_vel(krok, mnoznik)[0]) / 1.496e+8, np.array(pos_from_vel(krok, mnoznik)[1]) / 1.496e+8, np.array(pos_from_vel(krok, mnoznik)[2]) / 1.496e+8    

    data = np.column_stack((tt, x, y, z))

    with open("figures/earth_movement/new_pos_25_40.txt", "w") as f:
        f.write("#t[dni]   x [AU]    y [AU]    z [AU]\n")
        np.savetxt(f, data, fmt="%.4e")

#save_pos_to_file(25,40)



# ---- Wykres pozycji wyliczonej z prędkości ----
'''def plot_pos_from_vel():
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection="3d")

    #ax.plot(pos_from_vel()[0], pos_from_vel()[1], pos_from_vel()[2], lw=1, color="blue", label="Pozycja z prędkości")
    ax.plot(x, y, z, lw=1, color="blue", label="Pozycja z prędkości")

    ax.set_xlabel("X [au]")
    ax.set_ylabel("Y [au]")
    ax.set_zlabel("Z [au]")

    #ax.set_title(f"Pozycja punktu na Ziemi po {t*5} dniach obliczona z prędkości")
    ax.legend()

    set_axes_equal(ax)

    plt.show()

#plot_pos_from_vel()'''



# ---- Wykres pozycji z pliku ----
def plot_from_file(krok, mnoznik):
    data = np.loadtxt("figures/earth_movement/new_pos_25_40.txt", comments="#")
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

    set_axes_equal(ax)

    plt.show()


#plot_from_file(25,40) #krok i mnożnik

def plot_scatter():
    data = np.loadtxt("figures/tabela_v/EV_xyz_pelne.txt", skiprows=1)

    t = data[:,0]
    x = data[:,1]
    xf = data[:,2]
    y = data[:,4]
    yf = data[:,5]
    z = data[:,7]
    zf = data[:,8]

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')

    # Punkty pokolorowane czasem
    sc1 = ax.scatter(x, y, z, c=t, cmap='winter', s=8)
    sc2 = ax.scatter(xf, yf, zf, c=t, cmap='autumn', s=8)
    
    cbar1 = plt.colorbar(sc1, ax=ax, pad=0.02)
    cbar1.set_label("mjd")
    cbar1.ax.set_title("velocity", fontsize=10, pad=8)
    cbar2 = plt.colorbar(sc2, ax=ax, pad=0.12)
    cbar2.set_label("mjd")
    cbar2.ax.set_title("velocity fast", fontsize=10, pad=8)

    # Opcjonalnie połącz punkty linią
    #ax.plot(x, y, z, color='black', linewidth=0.5)

    #ax.scatter(x, y, z, s=8, label = 'ev')
    #ax.scatter(xf, yf, zf, s=8, label= 'ev fast')
    set_axes_equal(ax)

    ax.set_xlabel("X [m/s]")
    ax.set_ylabel("Y [m/s]")
    ax.set_zlabel("Z [m/s]")

    plt.show()

plot_scatter()