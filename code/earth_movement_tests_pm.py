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
    on_earth_itrs = SkyCoord(point_on_earth, frame=ITRS(obstime=t))
    # print("on_earth_itrs:", on_earth_itrs)
    
    #=== Sprawdzanie czy predkość obrotowa Ziemi jest dobrze liczona ===
    # ta zmienna nigdzie nie jest używana - ten fragment można chyba usunąć
    on_earth_velocity = CartesianDifferential(np.cross(omega_vec.to(1/u.s).value, on_earth_itrs.cartesian.xyz.to(u.m).value)*(u.m/u.s))
    # print("on_earth_velocity:", on_earth_velocity)

    #=== Transformacja położenia do układu ICRS ===
    on_earth_icrs = on_earth_itrs.transform_to(ICRS())
    # print("on_earth_icrs:", on_earth_icrs)

    #=== Koordynaty i prędkości Ziemi w układzie ICRS ===
    pos_e, vel_e = get_body_barycentric_posvel("earth", t)
    print("---")
    print("pos_e:", pos_e)
    print("vel_e:", vel_e)
    earth_icrs = SkyCoord(CartesianRepresentation(pos_e.xyz).with_differentials(CartesianDifferential(vel_e.xyz)),frame=ICRS(),obstime=t)
    print("^^^")
    print("earth_icrs:", earth_icrs)
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



def earth_velocity(mjd):
    """
    Calculate Earth's velocity in Galactocentric frame for given MJD.
    """
    _mjd = Time(mjd, format="mjd")
    return ICRS_to_Galactocentric(_mjd)[1].to(u.m/u.s)

def earth_velocity_pm(mjd):
    _mjd = Time(mjd, format="mjd")
    return earth_icrs_velocity_pm(_mjd).velocity.d_xyz.to(u.m/u.s)


def earth_icrs_velocity_pm(t):

    # położenie i prędkość środka Ziemi
    pos_e, vel_e = get_body_barycentric_posvel("earth", t)

    earth_icrs = SkyCoord(
        CartesianRepresentation(pos_e.xyz).with_differentials(
            CartesianDifferential(vel_e.xyz)
        ),
        frame=ICRS(),
        obstime=t
    )
    print("earth_icrs:", earth_icrs)

    return earth_icrs

from astropy.coordinates import (
    SkyCoord,
    Galactocentric,
    ITRS,
    CartesianRepresentation
)

import astropy.units as u
import numpy as np

def gal_velocity_direction_to_icrs(pos_gal, vel_gal):
    """
    Zamienia kierunek prędkości z Galactocentric na kierunek w ICRS.
    """

    speed = np.linalg.norm(vel_gal.to_value(u.km/u.s)) * u.km/u.s
    n_gal = vel_gal / speed

    eps = 1.0 * u.pc

    p0_gal = SkyCoord(
        CartesianRepresentation(pos_gal),
        frame=Galactocentric()
    )

    p1_gal = SkyCoord(
        CartesianRepresentation(pos_gal + eps * n_gal),
        frame=Galactocentric()
    )

    p0_icrs = p0_gal.transform_to(ICRS()).cartesian.xyz
    p1_icrs = p1_gal.transform_to(ICRS()).cartesian.xyz

    n_icrs = p1_icrs - p0_icrs
    n_icrs = n_icrs / np.linalg.norm(n_icrs.to_value())

    return n_icrs, speed


def normalize_quantity_vector(v):
    """
    Zwraca bezwymiarowy jednostkowy wektor numpy.
    """
    v_val = v.to_value(v.unit)
    return v_val / np.linalg.norm(v_val)


def gal_velocity_direction_to_icrs(pos_gal, vel_gal):
    speed = np.linalg.norm(vel_gal.to_value(u.km/u.s)) * u.km/u.s

    n_gal = vel_gal / speed

    eps = 1.0 * u.pc

    p0_gal = SkyCoord(
        CartesianRepresentation(pos_gal),
        frame=Galactocentric()
    )

    p1_gal = SkyCoord(
        CartesianRepresentation(pos_gal + eps * n_gal),
        frame=Galactocentric()
    )

    p0_icrs = p0_gal.transform_to(ICRS()).cartesian.xyz
    p1_icrs = p1_gal.transform_to(ICRS()).cartesian.xyz

    n_icrs_vec = p1_icrs - p0_icrs

    # ważne: bezwymiarowy numpy array
    n_icrs = normalize_quantity_vector(n_icrs_vec)

    return n_icrs, speed


def icrs_direction_to_itrs(n_icrs, t):
    """
    n_icrs: bezwymiarowy wektor numpy, np. array([nx, ny, nz])
    """
    if not isinstance(t, Time):
        t = Time(t, format="mjd")

    eps = 1.0 * u.pc

    p0_icrs = SkyCoord(
        CartesianRepresentation([0, 0, 0] * u.pc),
        frame=ICRS()
    )

    p1_icrs = SkyCoord(
        CartesianRepresentation(eps * n_icrs),
        frame=ICRS()
    )

    p0_itrs = p0_icrs.transform_to(ITRS(obstime=t)).cartesian.xyz
    p1_itrs = p1_icrs.transform_to(ITRS(obstime=t)).cartesian.xyz

    n_itrs_vec = p1_itrs - p0_itrs

    # znowu robimy bezwymiarowy kierunek
    n_itrs = normalize_quantity_vector(n_itrs_vec)

    return n_itrs

def vel(mjd):
    _mjd = Time(mjd, format="mjd")
    pos_gal, vel_gal = ICRS_to_Galactocentric(_mjd)
    n_icrs, speed = gal_velocity_direction_to_icrs(pos_gal, vel_gal)
    n_itrs = icrs_direction_to_itrs(n_icrs, _mjd)
    print("n_itrs:", n_itrs)
    print("speed:", speed)
    return n_itrs * speed



OMEGA_EARTH = 7.2921150e-5  # rad/s


def rotate_z(v, theta):
    """
    Obrót wektora wokół osi Z.
    """

    c = np.cos(theta)
    s = np.sin(theta)

    x, y, z = v

    return np.array([
        c*x + s*y,
        -s*x + c*y,
        z
    ])


def fast_itrs_direction(mjd, n0_itrs):
    """
    mjd:
        dowolny MJD

    n0_itrs:
        kierunek w ITRS policzony dla pełnego MJD
        np. dla floor(mjd)

    Zwraca:
        przybliżony kierunek w ITRS
        uwzględniający tylko obrót Ziemi.
    """

    mjd0 = np.floor(mjd)

    # część ułamkowa dnia
    frac_day = mjd - mjd0

    # sekundy od początku dnia
    dt_seconds = frac_day * 86400.0

    # obrót Ziemi
    theta = OMEGA_EARTH * dt_seconds

    n_rot = rotate_z(n0_itrs, theta)

    return n_rot / np.linalg.norm(n_rot)

if __name__ == "__main__":
    print("velocity:", vel(59000))
    print("velocity:", vel(59000.5))
    