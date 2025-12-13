from astropy.coordinates import get_body_barycentric_posvel, ICRS, get_body, SkyCoord, ITRS, Galactic,Galactocentric, CartesianRepresentation, CartesianDifferential
from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

print("=============================================")
t = Time(Time.now(), format="mjd")
#print(t)

def ITRS_sum_to_Galactocentric(time, w):
    
    #=== Koordynaty i prędkości Ziemi i Słońca w układzie barycentrycznym ICRS ===
    earth_pv = get_body_barycentric_posvel("earth", time)
    sun_pv = get_body_barycentric_posvel("sun", time)

    p_earth_icrs = earth_pv[0].xyz.to(u.m)
    p_sun_icrs = sun_pv[0].xyz.to(u.m)

    v_earth_icrs = earth_pv[1].xyz.to(u.m/u.s)
    v_sun_icrs = sun_pv[1].xyz.to(u.m/u.s)

    pos_icrs_earth = CartesianRepresentation(p_earth_icrs[0], p_earth_icrs[1], p_earth_icrs[2] )
    pos_icrs_sun = CartesianRepresentation(p_sun_icrs[0], p_sun_icrs[1], p_sun_icrs[2] )

    vel_icrs_earth = CartesianDifferential(v_earth_icrs[0], v_earth_icrs[1], v_earth_icrs[2])
    vel_icrs_sun = CartesianDifferential(v_sun_icrs[0], v_sun_icrs[1], v_sun_icrs[2])

    #=== Sum w układzie ICRS ===
    sum_pos_icrs = pos_icrs_earth + pos_icrs_sun
    sum_vel_icrs = vel_icrs_earth + vel_icrs_sun

    pos_coordinates = SkyCoord(sum_pos_icrs.with_differentials(sum_vel_icrs), frame=ICRS(), obstime=time)
    pos_coordinates_itrs = pos_coordinates.transform_to(ITRS(obstime=time))
    pos_coordinates_galactocentric = pos_coordinates_itrs.transform_to(Galactocentric())
    pos_galactocentric = pos_coordinates_galactocentric.cartesian.xyz

    sum_coordinates = SkyCoord(sum_pos_icrs.with_differentials(sum_vel_icrs), frame=ICRS(), obstime=time)
    sum_coordinates_itrs = sum_coordinates.transform_to(ITRS(obstime=time))
    sum_coordinates_galactocentric = sum_coordinates_itrs.transform_to(Galactocentric())
    vel_sum_galactocentric = sum_coordinates_galactocentric.velocity.d_xyz

    
    print("mjd= "+ str(time.value),"|  pos=" , pos_galactocentric.to(u.km*10**6),"|  vel=", vel_sum_galactocentric.to(u.km/u.s))

    #print("mjd="+ str(time.value), pos_galactocentric.to(u.km))
    #print("mjd="+ str(time.value), vel_sum_galactocentric.to(u.km/u.s))

    #print("mjd="+ str(time.value), "v={:.2f}".format(np.sqrt(vel_sum_galactocentric[0]**2 + vel_sum_galactocentric[1]**2 + vel_sum_galactocentric[2]**2).to(u.km/u.s)))

    if w==0:
        return pos_galactocentric
    if w==1:
        return vel_sum_galactocentric


#ITRS_sum_to_Galactocentric(t, 0)


def plot_3d_position():
    # listy na dane
    times = []
    x_list = []
    y_list = []
    z_list = []

    # obliczenia dla kolejnych dni
    for i in range(1, 60):
        time_future = Time(t.mjd + i*10, format="mjd")

        p = ITRS_sum_to_Galactocentric(time_future,0)
        x, y, z = p.to_value("10^6 km")

        times.append(time_future.mjd)
        x_list.append(x)
        y_list.append(y)
        z_list.append(z)

    # Rysowanie wykresu 3D
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(x_list, y_list, z_list)
    #ax.set_xlabel("X [km]")
    #ax.set_ylabel("Y [km]")
    #ax.set_zlabel("Z [km]")
    ax.set_title("Position vector trajectory")

    plt.show()


plot_3d_position()

#robi okrąg a powinna być spirala ?



def plot_3d_velocity_trajectory():

    # listy na dane
    times = []
    vx_list = []
    vy_list = []
    vz_list = []
    v_norm = []

    # obliczenia dla kolejnych dni
    for i in range(1, 60):
        time_future = Time(t.mjd + i, format="mjd")

        v = ITRS_sum_to_Galactocentric(time_future,1)
        vx, vy, vz = v.to_value("km/s")

        times.append(time_future.mjd)
        vx_list.append(vx)
        vy_list.append(vy)
        vz_list.append(vz)
        v_norm.append(np.sqrt(vx**2 + vy**2 + vz**2))

    # Rysowanie wykresu 3D
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(vx_list, vy_list, vz_list)
    ax.set_xlabel("v_x [km/s]")
    ax.set_ylabel("v_y [km/s]")
    ax.set_zlabel("v_z [km/s]")
    ax.set_title("Velocity vector trajectory")

    plt.show()

#plot_3d_velocity_trajectory()



