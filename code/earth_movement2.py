import numpy as np 
import astropy

from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris, SkyCoord, CartesianDifferential, CartesianRepresentation
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_body_barycentric_posvel, solar_system_ephemeris, ITRS
from astropy import units as u
from astropy.coordinates import ( 
    AltAz,
    BarycentricTrueEcliptic,
    EarthLocation,
    Galactic,
    GCRS,
    GeocentricMeanEcliptic,
    get_body,
    get_sun,
    SkyCoord,
    ICRS,
    get_body_barycentric,
)

def sun_speed_astropy(mjd):

    """
    Gives a Cartesian velocity vector of the speed of the Sun around the Galaxy in Earth-Centered Earth-Fixed coordinates (ECEF)
    :param mjd - Modified Julian Date
    :return: [vx, vy, vz] in km/s

    """
    time = Time(mjd, format='mjd')

    # Location of the observer on the Earth's surface (can be customized)
    location = EarthLocation(lat=0*u.deg, lon=0*u.deg, height=0*u.m)  # At the equator, the zero meridian

    # We create an AltAz coordinate system for this location and time = how the sun appears in the sky from a particular location
    altaz = AltAz(obstime=time, location=location)

    # Transform to galactic coordinates
    sun_galactic = SkyCoord(l=0*u.deg, b=0*u.deg, frame='galactic')

    # We convert to the AltAz layout (this is our approximation for the ECEF layout)
    sun_ecef = sun_galactic.transform_to(altaz)
    print(sun_ecef)

    # The velocity vector of the Sun in the Galaxy (approximately)
    v_sun_galactic = 230 * u.km/u.s  # 230 km/s w kierunku l=90°

    # Now we need to transform the velocity vector to the ECEF (AltAz) system
    # We will do this in an approximate way, assuming that the azimuthal component is the X component and the elevation component is the Z component.
    v_sun_ecef = [
        v_sun_galactic * np.cos(sun_ecef.az.to(u.rad).value) * np.cos(sun_ecef.alt.to(u.rad).value),
        v_sun_galactic * np.sin(sun_ecef.az.to(u.rad).value) * np.cos(sun_ecef.alt.to(u.rad).value),
        v_sun_galactic * np.sin(sun_ecef.alt.to(u.rad).value)
    ]

    return [v.value for v in v_sun_ecef]

# Test

def earth_velocity_xyz_astropy(mjd):

    """
    Gives a Cartesian velocity vector of the speed of the Earth around the Galaxy and the Sun in Earth-Centered Earth-Fixed coordinates (ECEF)
    :param mjd - Modified Julian Date
    :return: [vx, vy, vz] in km/s
    
    """
    
    # We define the time in MJD format
    time = Time(mjd, format='mjd')
    
    # Location of the observer on the Earth's surface (can be customized)
    location = EarthLocation(lat=0*u.deg, lon=0*u.deg, height=0*u.m)  # At the equator, the zero meridian

    # We create an AltAz coordinate system for this location and time = how the sun appears in the sky from a particular location
    altaz = AltAz(obstime=time, location=location)

    # Get speed of the Sun around Galaxy in Earth-Centered coordinates
    v_sun_galaxy = sun_speed_astropy(mjd)

    print("Sun vel:", v_sun_galaxy)
    
    # Get the barycentric position and velocity of the Earth
    earth_barycentric_pos, earth_barycentric_vel = get_body_barycentric_posvel('earth', time)

    print("get body:", earth_barycentric_pos, earth_barycentric_vel)
    
    # Transform barycentric position to AltAz (ECEF)
    earth_barycentric_pos_obj = SkyCoord(earth_barycentric_pos, frame='icrs')
    earth_altaz = earth_barycentric_pos_obj.transform_to(altaz)
    print("Earth pos", earth_altaz)

    # Extract velocity components in ITRS frame (ECEF)
    earth_barycentric_vel_obj = SkyCoord(earth_barycentric_vel, frame='icrs')
    earth_vel_altaz = earth_barycentric_vel_obj.transform_to(altaz)
    print("Earth vel:", earth_vel_altaz)

    # Total velocity in ECEF by adding Sun's velocity around Galaxy in ECEF
    total_velocity = earth_vel_altaz + v_sun_galaxy

    return total_velocity

if __name__ == "__main__":

    mjd = 60493.60420139
    sun_speed = sun_speed_astropy(mjd)
    print("sun speed", sun_speed)

    earth_speed = earth_velocity_xyz_astropy(mjd)
    print("earth speed", earth_speed)
