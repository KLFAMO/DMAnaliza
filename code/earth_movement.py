import math
import numpy as np
import re

from astropy.time import Time
import astropy.units as u
import astropy.coordinates as coord

from astropy.coordinates import ( 
    AltAz,
    BarycentricMeanEcliptic,
    BarycentricTrueEcliptic,
    EarthLocation,
    Galactic,
    GCRS,
    GeocentricMeanEcliptic,
    get_body,
    get_sun,
    SkyCoord, 
)


import spiceypy


def dms2dd(lat):
    """
    Converts GPS coordinates from hours, minutes, second to decimal degrees format
    
    :Changes:
        2023-08-20 by Joséphine Strübing-Tardy: First version
    """
    deg, minutes, seconds, direction =  re.split('[°\'"]', lat)
    newlat = (float(deg) + float(minutes)/60 + float(seconds)/(60*60)) * (-1 if direction in ['W', 'S'] else 1)
    return newlat

def dd2dms(dd):
    """
    Converts coordinates from decimal degrees to hours, minute, seconds
    
    :Changes:
        2023-08-20 by Joséphine Strübing-Tardy: First version
    """
    d = int(dd)
    m = int((dd - d) * 60)
    s = (dd - d - m/60) * 3600.00
    z= round(s, 2)
    if d >= 0:
        dms = (f"+ {abs(d)}° {abs(m)}' {abs(z)}\"")
    else:
        dms = (f"- {abs(d)}° {abs(m)}' {abs(z)}\"")
    return dms

def gps2ecef(long, lat, h):
    """
    Converts GPS coordinates into ECEF (Earth-Centered Earth-Fixed) Cartesian coordinates

    :param Longitude (decimal degrees), 
    :param Latitude (decimal degrees), 
    :param Height (m)
    :return: X, Y, Z in the ECEF

    :Changes:
        2023-08-20 by Joséphine Strübing-Tardy: First version
    """
    a = 6378137 # earth semi-major axis = equatorial radius in m
    b = 6356752 # earth sem-minor axis = polar radius in m
    lat = math.radians(lat)
    long = math.radians(long)
    
    e = 1 - (b**2) / (a**2)
    N = a / np.sqrt(1 - e**2 * math.sin(lat)**2 )
    
    x = (N+h) * math.cos(lat) * math.cos(long)
    y = (N+h) * math.cos(lat) * math.sin(long)
    z = (N * (1-e**2) +h) * math.sin(lat)
    
    return x, y, z

# # checking for Torun 
# long = dms2dd("18° 35' 53.30' E")
# lat = dms2dd("53° 00' 49.50' N")
# print("GPS coordinates in decimal degrees:",long, lat)

# h = 65 # above sea-level, in m 
# x, y, z = gps2ecef(long, lat, h)
# print("X:", x)
# print("Y:", y)
# print("Z:", z)

def sun_ecliptic_position(julian_date):
    """
    Returns the position of the Sun in ecliptic coordinates for given Julian Date

    :param julian_date
    :return: l - ecliptic longitude
    :rerutn: b - ecliptic latitude
    :return: r - sun-earth distance [au]
    :return: e - obliquity of the ecliptic
    :return: r_km - sun-earth distance [km]

    :Changes:
        2023-08-20 by Joséphine Strübing-Tardy: First version 
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
    # obliquity of the ecliptic 
    e = 23.439 - 0.0000003 * n 
    
    return l, b, r, e, r_km

def sun_ecliptic_position_astropy(julian_date):
    time = Time(julian_date, format="jd")
    #get coordinates of sun
    sun_coords = get_body('sun', time)
    # sun_ecliptic = sun_coords.transform_to(BarycentricMeanEcliptic())
    sun_ecliptic = sun_coords.transform_to(GeocentricMeanEcliptic())
    l = sun_ecliptic.lon.degree
    b = sun_ecliptic.lat.degree
    r = sun_ecliptic.distance.au
    r_km = sun_ecliptic.distance.to(u.km).value
    e = 23.439292 
    return l, b, r, e, r_km
 
def sun_equatorial_position(l_deg, b_deg, r, e_deg):
    """
    Returns position of the Sun in equatorial coordinates for a given set of ecliptic coordinates

    :Changes:
        2023-08-20 by Joséphine Strübing-Tardy: First version
    """
    # Convert input values to radians
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    e = np.deg2rad(e_deg)
    
    # Calculate the equatorial coordinates
    ra = np.arctan2(np.cos(e) * np.sin(l), np.cos(l))
    dec = np.arcsin(np.sin(e) * np.sin(l))
    
    return np.rad2deg(ra), np.rad2deg(dec)  # Convert back to degrees

# l, b, r, e, r_km = sun_ecliptic_position(12000)
# print(sun_equatorial_position(l, b, r, e))

def earth_velocity_vector(julian_date):
    """
    Computes Earth's Vx, Vy, Vz in the Earth-Centered, Earth-Fixed coordinates for given Julian Date

    :Changes:
        2023-08-20 by Joséphine Strübing-Tardy: First version
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

def earth_velocity_xyz(mjd):
    """
    Computes Earth's velocity in space in the Earth-Centered, Earth-Fixed coordinates for given Modified Julian Date

    :param mjd - Modified Julian Date
    :return: [vx, vy, vz] in km/s
    :Changes:
        2023-08-20 by Joséphine Strübing-Tardy: First version
    """
    jd = Time(mjd, format='mjd').jd
    vx, vy, vz = earth_velocity_vector(jd)
    return [vx.value, vy.value, vz.value] # km/s


from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris, SkyCoord
from astropy.time import Time
from astropy import units as u

def earth_velocity_xyz_astropy(mjd_time):
    print('using astropy -----------------------')
    # Definiujemy czas w formacie MJD
    time = Time(mjd_time, format='mjd')
    
    # Uzyskanie współrzędnych Słońca w układzie GCRS (Geocentric Celestial Reference System)
    sun_gcrs = get_body('sun', time)
    # Konwersja na geocentryczne współrzędne ekliptyczne
    sun_ecliptic = sun_gcrs.transform_to(BarycentricTrueEcliptic())
    print('sun_ecliptic lon: ', sun_ecliptic.lon)
    print('sun_ecliptic lat: ', sun_ecliptic.lat)

def icrs_sun_astropy(jd):
    time = Time(jd, format='jd')

    # Pobranie współrzędnych Słońca w układzie ICRS
    sun = get_sun(time)

    # Współrzędne Słońca
    ra = sun.ra  # Rektascensja
    dec = sun.dec  # Deklinacja
    distance = sun.distance  # Odległość
    print(sun)


def sun_speed_astropy(mjd):
    t = Time(mjd, format='mjd')

    # Lokalizacja obserwatora na powierzchni Ziemi (można dostosować do własnych potrzeb)
    location = EarthLocation(lat=0*u.deg, lon=0*u.deg, height=0*u.m)  # Na równiku, południk zerowy

    # Tworzymy układ współrzędnych AltAz dla tej lokalizacji i czasu
    altaz = AltAz(obstime=t, location=location)

    # Współrzędne Słońca w układzie Galaktyki
    sun_galactic = SkyCoord(l=90*u.deg, b=0*u.deg, frame='galactic')

    # Przekształcamy do układu AltAz (to jest nasze przybliżenie dla układu ECEF)
    sun_ecef = sun_galactic.transform_to(altaz)

    # Wektor prędkości Słońca w Galaktyce (w przybliżeniu)
    v_sun_galactic = 230 * u.km/u.s  # 230 km/s w kierunku l=90°

    # Teraz musimy przekształcić wektor prędkości do układu ECEF (AltAz)
    # Zrobimy to w przybliżony sposób, zakładając, że składowa azymutalna to składowa X, a składowa wysokości to składowa Z.
    v_sun_ecef = [
        v_sun_galactic * np.cos(sun_ecef.az.to(u.rad).value) * np.cos(sun_ecef.alt.to(u.rad).value),
        v_sun_galactic * np.sin(sun_ecef.az.to(u.rad).value) * np.cos(sun_ecef.alt.to(u.rad).value),
        v_sun_galactic * np.sin(sun_ecef.alt.to(u.rad).value)
    ]

    return [v.value*1e3 for v in v_sun_ecef]
