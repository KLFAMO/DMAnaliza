import math
import numpy as np
import re

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord, CartesianDifferential
import matplotlib.pyplot as plt
#!pip install --upgrade numpy==1.20.0
import numpy
print(numpy.__version__)

import gala
import gala.dynamics as gd
import gala.potential as gp

import spiceypy
from skyfield.api import load


def sun_ecliptic_position(julian_date):
    """
    Returns the position of the Sun in ecliptic coordinates for given Julian Date 
    """
    n = julian_date - 2451545 # number of days since J2000
    
    # mean longitude of the Sun 
    L = 280.459 + 0.9856474 * n 
    
    # mean anomaly of the Sun 
    g = 357.529 + 0.9856003 * n 
    
    # get the values of L and g to be between 0 and 360 degrees
    L = L % 360
    g = g % 360
    
    # ecliptic longitude 
    l = L + 1.915 * np.sin(np.deg2rad(g)) + 0.020 * np.sin(np.deg2rad(2 * g))
    
    # ecliptic latitude 
    b = 0
    
    # sun-earth distance in au 
    r = 1.00014 - 0.01671 * np.cos(np.deg2rad(g)) - 0.00014 * np.cos(np.deg2rad(2 * g))
    
    r_km = r*1.496e8
    print("Sun-Earth distance in km:", r_km) # sanity check
    
    # obliquity of the ecliptic 
    e = 23.439 - 0.0000003 * n 
    
    return l, b, r, e, r_km
 
def sun_equatorial_position(l_deg, b_deg, r, e_deg):
    """
    Returns position of the Sun in equatorial coordinates for a given set of ecliptic coordinates
    """
    
    # Convert input values to radians
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    e = np.deg2rad(e_deg)
    
    # Calculate the equatorial coordinates
    ra = np.arctan2(np.cos(e) * np.sin(l), np.cos(l))
    dec = np.arcsin(np.sin(e) * np.sin(l))
    
    return np.rad2deg(ra), np.rad2deg(dec)  # Convert back to degrees

def earth_velocity_vector(julian_date):
    """
    Computes Earth's Vx, Vy, Vz in the Earth-Centered, Earth-Fixed coordinates for every given Julian Date
    """
    # VELOCITY OF SOLAR SYSTEM AROUND GALAXY ----------------------------------------------------
    
    # Get the ecliptic coordinates of the Sun for this date 
    l_deg, b_deg, r, e_deg, r_km = sun_ecliptic_position(julian_date)

    # Calculate equatorial position (in DD) for this ecliptic position
    ra_deg, dec_deg = sun_equatorial_position(l_deg, b_deg, r, e_deg)
    
    # Position of the Sun on this Julian date in IRCS coordinates 
    ra = "9 h 16 m 27.73590156690716 s"
    dec = dd2dms(dec_deg)
    icrs = coord.SkyCoord(ra=coord.Angle(f'{ra}'), dec=coord.Angle(dec), 
                          distance = r_km*u.km) 

    # Converting to hour angle format for Right Ascension 
    ra_h = int(ra_deg // 15)  # Convert degrees to hours (1 hour = 15 degrees)
    ra_m = int((ra_deg % 15) * 4)  # Convert remaining degrees to minutes (1 min = 0.25 degrees)
    ra_s = (ra_deg % 15 % 0.25) * 240  # Convert remaining degrees to seconds (1 sec = 0.00416667 degrees)
    
    ra = f"{ra_h} h {ra_m} m {ra_s} s"

    # Get position of the Sun for this Julian Date in galactocentric coordinates 
    galcen_frame = coord.Galactocentric()
    galcen = icrs.transform_to(galcen_frame)
    
   # Coordinates and velocity of the Sun in galactocentric coordinates 
    galactocentric = coord.Galactocentric(x = galcen.x, 
                                      y = galcen.y, 
                                      z = galcen.z, 
                                      v_x = 12.9 * u.km/u.s, 
                                      v_y = 245.6 * u.km/u.s, 
                                      v_z = 7.78 * u.km/u.s)
    # Transform the velocities to ecliptic coordinates 
    heliocentric_ecliptic = galactocentric.transform_to(coord.HeliocentricTrueEcliptic)

    # Extract heliocentric velocities in ecliptic coordinates
    
    vx__gal_e = heliocentric_ecliptic.velocity.d_x
    vy_gal_e = heliocentric_ecliptic.velocity.d_y
    vz_gal_e = heliocentric_ecliptic.velocity.d_z
    
    # Convert heliocentric ecliptic velocities to IRCS 
    icrs = heliocentric_ecliptic.transform_to(coord.ICRS)

    # Convert ICRS Velocities to Geocentric ECEF Velocities
    ecef_velocities = icrs.velocity.d_xyz

    vx_gal_ECEF = ecef_velocities[0]
    vy_gal_ECEF = ecef_velocities[1] 
    vz_gal_ECEF = ecef_velocities[2]
    
    #  VELOCITY OF EARTH AROUND SUN --------------------------------------------------------------
    
    # Obtain velocity of Earth around the Sun in ecliptic coordinates for a given julian date
    spiceypy.furnsh('de421.bsp')

    eart_state_wrt_sun, earth_sun_lt = spiceypy.spkgeo(targ=399, \
                                                    et=julian_date, \
                                                    ref='ECLIPJ2000', obs=10)

    # Transform the velocities to IRCS velocities 
    heliocentric_ecliptic = coord.HeliocentricTrueEcliptic(representation_type='cartesian', 
                                      x = eart_state_wrt_sun[0] * u.km, 
                                      y = eart_state_wrt_sun[1] * u.km, 
                                      z = eart_state_wrt_sun[2] * u.km, 
                                      differential_type='cartesian',
                                      v_x = eart_state_wrt_sun[3] * u.km/u.s, 
                                      v_y = eart_state_wrt_sun[4] * u.km/u.s, 
                                      v_z = eart_state_wrt_sun[5] * u.km/u.s 
                                      )
    icrs = heliocentric_ecliptic.transform_to(coord.ICRS)

    # Convert ICRS Velocities to Geocentric ECEF Velocities
    ecef_velocities_sun = icrs.velocity.d_xyz
    
    vx_sun_ECEF = ecef_velocities_sun[0]
    vy_sun_ECEF = ecef_velocities_sun[1] 
    vz_sun_ECEF = ecef_velocities_sun[2]
    
    # COMPILING THE VELOCITIES ---------------------------------------------------------------------
    
    vx_ECEF = vx_gal_ECEF + vx_sun_ECEF 
    vy_ECEF = vy_gal_ECEF + vy_sun_ECEF
    vz_ECEF = vz_gal_ECEF + vz_sun_ECEF
    
    return vx_ECEF, vy_ECEF, vz_ECEF

# Test 

julian_date = 2460166

vx, vy, vz = earth_velocity_vector(julian_date)
print("Geocentric ECEF Velocities:")
print("vx :", vx)
print("vy :", vy)
print("vz :", vz)