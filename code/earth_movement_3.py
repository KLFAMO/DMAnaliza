from astropy.coordinates import ICRS, SkyCoord, ITRS, Galactocentric, CartesianRepresentation, CartesianDifferential, get_body_barycentric_posvel
from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import parameters

print("=============================================")
#t = Time(Time.now(), format="mjd")
t = Time(59000.0, format="mjd")
#print(t)
omega_vec = [0, 0, 7.2921150e-5] /u.s # Earth's angular velocity

point_on_earth = CartesianRepresentation(parameters.inf['UMK1']['X']*u.m,
                                                        parameters.inf['UMK1']['Y']*u.m,
                                                        parameters.inf['UMK1']['Z']*u.m)


def ITRS_to_ICRS(t):
    #=== Koordynaty miejsca na Ziemi w układzie ITRS (to samo co ITRF) ===
    on_earth_itrs = SkyCoord(point_on_earth, frame=ITRS(obstime=t))
    
    #=== Sprawdzanie czy predkość obrotowa Ziemi jest dobrze liczona ===
    #on_earth_velocity = CartesianDifferential(np.cross(omega_vec.to(1/u.s).value, on_earth_itrs.cartesian.xyz.to(u.m).value)*(u.m/u.s))
    #print("On Earth ITRS position:", on_earth_itrs.cartesian.xyz.to(u.km))
    #print("On Earth velocity:", on_earth_velocity.d_xyz.to(u.km/u.s))
    #print("On Earth velocity norm:", np.sqrt(on_earth_velocity.d_x**2 + on_earth_velocity.d_y**2 + on_earth_velocity.d_z**2).to(u.m/u.s))
    #print()


    #=== Transformacja położenia do układu ICRS ===
    on_earth_icrs = on_earth_itrs.transform_to(ICRS())

    #=== Koordynaty i prędkości Ziemi w układzie ICRS ===
    pos_e, vel_e = get_body_barycentric_posvel("earth", t)
    earth_icrs = SkyCoord(CartesianRepresentation(pos_e.xyz).with_differentials(CartesianDifferential(vel_e.xyz)),frame=ICRS(),obstime=t)

    #=== Obliczanie prędkości obrotowej w układzie ICRS ===
    r_icrs = (on_earth_icrs.cartesian.xyz - earth_icrs.cartesian.xyz).to(u.m)
    v_rot_icrs = np.cross(omega_vec, r_icrs).to(u.m/u.s)

    print("On Earth ICRS - Earth ICRS position:", r_icrs.to(u.km))
    #print("On Earth rotational velocity in ICRS:", v_rot_icrs.to(u.km/u.s))
    print("On Earth rotational velocity norm in ICRS:", np.sqrt(v_rot_icrs[0]**2 + v_rot_icrs[1]**2 + v_rot_icrs[2]**2).to(u.m/u.s))
    print()

    #=== Całkowita prędkość miejsca na Ziemi w układzie ICRS ===
    v_earth_icrs = vel_e.xyz.to(u.m/u.s)
    v_total_icrs = v_earth_icrs + v_rot_icrs

    print("Prędkość orbitalna w ICRS:", np.linalg.norm(v_earth_icrs).to(u.km/u.s))
    print("Prędkość obrotowa w ICRS:", np.linalg.norm(v_rot_icrs).to(u.km/u.s))
    print("Całkowita prędkość w ICRS:", np.linalg.norm(v_total_icrs).to(u.km/u.s))
    print()

    return on_earth_icrs, r_icrs, v_total_icrs


#ITRS_to_ICRS(t)



def ICRS_to_Galactocentric(t):
    # --- bierzemy wynik z pierwszej funkcji ---
    on_earth_icrs, r_icrs, v_total_icrs = ITRS_to_ICRS(t)

    # --- tworzymy SkyCoord w ICRS z PRĘDKOŚCIĄ ---
    coord_icrs = SkyCoord(
        CartesianRepresentation(
            on_earth_icrs.cartesian.xyz
        ).with_differentials(
            CartesianDifferential(v_total_icrs)
        ),
        frame=ICRS(),
        obstime=t
    )

    # --- transformacja do Galactocentric ---
    coord_gal = coord_icrs.transform_to(Galactocentric())

    # --- wyciągamy wektory ---
    pos_gal = coord_gal.cartesian.xyz
    vel_gal = coord_gal.velocity.d_xyz

    # --- wypisywanie kontrolne ---
    print("Pozycja w Galactocentric:", pos_gal.to(u.kpc))
    print("Prędkość w Galactocentric:", vel_gal.to(u.km/u.s))
    print("Prędkość w Galactocentric (norma):", np.linalg.norm(vel_gal.to_value(u.km/u.s)) * u.km/u.s)

    return pos_gal, vel_gal

ICRS_to_Galactocentric(t)
