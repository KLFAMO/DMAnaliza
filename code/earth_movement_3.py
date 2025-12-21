from astropy.coordinates import get_body_barycentric_posvel, ICRS, get_body, SkyCoord, ITRS, CIRS, Galactocentric, CartesianRepresentation, CartesianDifferential
from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import parameters

print("=============================================")
#t = Time(Time.now(), format="mjd")
t = Time(59000.0, format="mjd")
#print(t)
omega_vec = [0, 0, 7.2921150e-5] * u.rad/u.s # Earth's angular velocity

point_on_earth = CartesianRepresentation(parameters.inf['UMK1']['X']*u.m,
                                                        parameters.inf['UMK1']['Y']*u.m,
                                                        parameters.inf['UMK1']['Z']*u.m)


def ITRS_to_CIRS(t):
    #=== Koordynaty miejsca na Ziemi w układzie ITRS (to samo co ITRF) ===
    on_earth_itrs = SkyCoord(point_on_earth, frame=ITRS(obstime=t))

    on_earth_velocity = CartesianDifferential(np.cross(omega_vec.to(u.rad/u.s).value, on_earth_itrs.cartesian.xyz.to(u.m).value)*(u.m/u.s))

    print("On Earth ITRS position (m):", on_earth_itrs.cartesian.xyz.to(u.m))
    print("On Earth ITRS velocity (m/s):", on_earth_velocity.d_xyz.to(u.m/u.s))
    print("On Earth ITRS velocity norm (m/s):", np.sqrt(on_earth_velocity.d_x**2 + on_earth_velocity.d_y**2 + on_earth_velocity.d_z**2).to(u.m/u.s))
    print()
    
    #=== Koordynaty i prędkości Ziemi w układzie CIRS ===
    earth_cirs = on_earth_itrs.transform_to(CIRS())

    data_with_vel = SkyCoord(point_on_earth.with_differentials(on_earth_velocity), frame=ITRS(obstime=t))
    
    earth_cirs_velocity = data_with_vel.transform_to(CIRS())

    #earth_cirs_velocity = SkyCoord(on_earth_itrs.point_on_earth.with_differentials(on_earth_velocity), frame=CIRS())

    print("On Earth CIRS position (m):", earth_cirs.cartesian.xyz.to(u.m))
    print("Earth CIRS velocity (m/s):", earth_cirs_velocity.velocity.d_xyz.to(u.m/u.s))
    print("Earth CIRS velocity norm (m/s):", np.sqrt(earth_cirs.velocity.d_x**2 + earth_cirs.velocity.d_y**2 + earth_cirs.velocity.d_z**2).to(u.m/u.s))
    print()

    return earth_cirs, earth_cirs_velocity

ITRS_to_CIRS(t)




def CIRS_to_ICRS(t):
    earth_cirs = ITRS_to_CIRS(t)
    
    #=== Koordynaty Ziemi w układzie ICRS ===
    earth_icrs = earth_cirs.transform_to(ICRS())

    print("Earth ICRS position (au):", earth_icrs.cartesian.xyz.to(u.au))
    #print("Earth ICRS velocity (m/s):", earth_icrs.velocity.d_xyz.to(u.km/u.s))
    #print("Earth ICRS velocity norm (km/s):", np.sqrt(earth_icrs.velocity.d_x**2 + earth_icrs.velocity.d_y**2 + earth_icrs.velocity.d_z**2).to(u.km/u.s))
    print()

    return earth_icrs

#CIRS_to_ICRS(t)

def ICRS_to_Galactocentric(t):
    earth_icrs = CIRS_to_ICRS(t)
    
    #=== Koordynaty i prędkości Ziemi w układzie Galactocentric ===
    earth_galactocentric = earth_icrs.transform_to(Galactocentric())

    #print("Earth Galactocentric position (kpc):", earth_galactocentric.cartesian.xyz.to(u.kpc))
    print("Earth Galactocentric position (au):", earth_galactocentric.cartesian.xyz.to(u.au))
    print("Earth Galactocentric velocity (km/s):", earth_galactocentric.velocity.d_xyz.to(u.km/u.s))
    print("Earth Galactocentric velocity norm (km/s):", np.sqrt(earth_galactocentric.velocity.d_x**2 + earth_galactocentric.velocity.d_y**2 + earth_galactocentric.velocity.d_z**2).to(u.km/u.s))

    return earth_galactocentric

#ICRS_to_Galactocentric(t)