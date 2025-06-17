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

import earth_movement as em

print("==============================================")

current_jd = Time.now()
print(current_jd)
print()

#Returns the position of the Sun in ecliptic coordinates for given Julian Date
    # julian_date
    # l - ecliptic longitude
    # b - ecliptic latitude
    # r - sun-earth distance [au]
    # e - obliquity of the ecliptic
    # r_km - sun-earth distance [km]

#julian_date = 2460780

time = Time(current_jd, format="jd")
print(time)
#get coordinates of sun
sun_coords = get_body('sun', time)
# print()
# print(sun_coords)
# sun_ecliptic = sun_coords.transform_to(BarycentricMeanEcliptic())
sun_ecliptic = sun_coords.transform_to(GeocentricMeanEcliptic())
# print()
# print(sun_ecliptic)
l = sun_ecliptic.lon.degree
b = sun_ecliptic.lat.degree
r = sun_ecliptic.distance.au
r_km = sun_ecliptic.distance.to(u.km).value
e = 23.439292

print()
print('Ecliptic longitude of Sun [deg] =',l)
print('Ecliptic latitude of Sun [deg] =', b)
print('sun-earth distance [au] =', r)
#print(e, r_km)
print()


#print(em.earth_velocity_vector(current_jd))

#mjd = current_jd - 2400001
#print(em.earth_velocity_xyz(mjd))

