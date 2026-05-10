from astropy.coordinates import ICRS, SkyCoord, ITRS, Galactocentric, CartesianRepresentation, CartesianDifferential, get_body_barycentric_posvel, galactocentric_frame_defaults
from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

print("========================================================================================")

t = Time(Time.now(), format="mjd")
omega_vec = [0, 0, 7.2921150e-5] /u.s # Earth's angular velocity
point_on_earth = CartesianRepresentation(0*u.m,  0*u.m,  0*u.m)

#t = Time(59000.0, format="mjd")

def ITRS_to_ICRS(t):
    #=== Koordynaty miejsca na Ziemi w układzie ITRS ===
    t = Time(t, format="mjd")
    on_earth_itrs = SkyCoord(point_on_earth, frame=ITRS(obstime=t))
    
    #=== Sprawdzanie czy predkość obrotowa Ziemi jest dobrze liczona ===
    on_earth_velocity = CartesianDifferential(np.cross(omega_vec.to(1/u.s).value, on_earth_itrs.cartesian.xyz.to(u.m).value)*(u.m/u.s))

    #=== Transformacja położenia do układu ICRS ===
    on_earth_icrs = on_earth_itrs.transform_to(ICRS())

    #=== Koordynaty i prędkości Ziemi w układzie ICRS ===
    pos_e, vel_e = get_body_barycentric_posvel("earth", t)
    earth_icrs = SkyCoord(CartesianRepresentation(pos_e.xyz).with_differentials(CartesianDifferential(vel_e.xyz)),frame=ICRS(),obstime=t)
    #cr = CartesianRepresentation(pos_e.xyz).with_differentials(CartesianDifferential(vel_e.xyz))
    #earth_icrs = SkyCoord(cr, frame=ICRS(), obstime=t)

    #=== Obliczanie prędkości obrotowej w układzie ICRS ===
    r_icrs = (on_earth_icrs.cartesian.xyz - earth_icrs.cartesian.xyz).to(u.m)
    v_rot_icrs = np.cross(omega_vec, r_icrs).to(u.m/u.s)

    #=== Całkowita prędkość miejsca na Ziemi w układzie ICRS ===
    v_earth_icrs = vel_e.xyz.to(u.m/u.s)
    v_total_icrs = v_earth_icrs + v_rot_icrs


    return on_earth_icrs, r_icrs, v_total_icrs



def ICRS_to_Galactocentric(t):
    # --- bierzemy wynik z pierwszej funkcji ---
    t = Time(t, format="mjd")
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
    

    return pos_gal, vel_gal


# def earth_vel_vector(mjd):
#     """
#     Calculate Earth's velocity in Galactocentric frame for given MJD.
#     """
#     return ICRS_to_Galactocentric(mjd)[0].to(u.m/u.s)


def earth_velocity(mjd):
    """
    Calculate Earth's velocity in Galactocentric frame for given MJD.
    """
    #return ICRS_to_Galactocentric(mjd)[1].to(u.m/u.s)
    return ICRS_to_Galactocentric(mjd)[1].to(u.m/u.s)


if __name__ == "__main__":
    print("Położenie i prędkość w Galactocentric:", ICRS_to_Galactocentric(t)[0].to(u.pc), ICRS_to_Galactocentric(t)[1].to(u.km/u.s))
    print("Norma prędkości w Galactocentric:", np.linalg.norm(ICRS_to_Galactocentric(t)[1].to_value(u.km/u.s)) * u.km/u.s)
    print(earth_velocity(t))
    print(np.linalg.norm(earth_velocity(t)))