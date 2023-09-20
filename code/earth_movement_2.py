from astropy.coordinates import solar_system_ephemeris, EarthLocation, get_sun
from astropy.time import Time

t = Time("2023-09-20T12:00:00")

with solar_system_ephemeris.set('builtin'):
    sun = get_sun(t)

print(sun.cartesian.xyz)