# import numpy as np
# import re

# from astropy.time import Time
# import astropy.units as u
# import astropy.coordinates as coord

# from astropy.coordinates import ( 
#     AltAz,
#     BarycentricMeanEcliptic,
#     BarycentricTrueEcliptic,
#     EarthLocation,
#     Galactic,
#     GCRS,
#     GeocentricMeanEcliptic,
#     get_body,
#     get_sun,
#     SkyCoord, 
# )

# import earth_movement as em

# print("==============================================")

# current_jd = Time.now()
# print('jd =', current_jd)
# print()

# #Returns the position of the Sun in ecliptic coordinates for given Julian Date
#     # julian_date
#     # l - ecliptic longitude
#     # b - ecliptic latitude
#     # r - sun-earth distance [au]
#     # e - obliquity of the ecliptic
#     # r_km - sun-earth distance [km]

# #julian_date = 2460780

# time = Time(current_jd, format="jd")
# print('time =', time)
# #get coordinates of sun
# sun_coords = get_body('sun', time)
# # print()
# # print(sun_coords)
# # sun_ecliptic = sun_coords.transform_to(BarycentricMeanEcliptic())
# sun_ecliptic = sun_coords.transform_to(GeocentricMeanEcliptic())
# # print()
# # print(sun_ecliptic)
# l = sun_ecliptic.lon.degree
# b = sun_ecliptic.lat.degree
# r = sun_ecliptic.distance.au
# r_km = sun_ecliptic.distance.to(u.km).value
# e = 23.439292

# print()
# print('Ecliptic longitude of Sun [deg] =',l)
# print('Ecliptic latitude of Sun [deg] =', b)
# print('sun-earth distance [au] =', r)
# #print(e, r_km)
# print()


# #print(em.earth_velocity_vector(current_jd))

# #mjd = current_jd - 2400001
# #print(em.earth_velocity_xyz(mjd))

# =============================================

from astropy.coordinates import get_body_barycentric_posvel, ICRS, get_body, SkyCoord, ITRS, Galactic,Galactocentric, CartesianRepresentation, CartesianDifferential
from astropy.time import Time
import astropy.units as u
import numpy as np



current_jd = Time.now()
time = Time(current_jd, format="jd")

un_coords = get_body('earth', time)


#=== Koordynaty i prędkości Ziemi w układzie barycentrycznym ICRS ===
earth_pv = get_body_barycentric_posvel("earth", time)
sun_pv = get_body_barycentric_posvel("sun", time)

print()
pozycja_earth_icrs = earth_pv[0].xyz.to(u.m)
print("Ziemia - pozycja w ICRS (m):", pozycja_earth_icrs)

v_earth_icrs = earth_pv[1].xyz.to(u.m/u.s)
print("Ziemia - prędkość w ICRS (m/s):", v_earth_icrs)
v_sun_icrs = sun_pv[1].xyz.to(u.m/u.s)

pos_icrs = CartesianRepresentation(pozycja_earth_icrs[0], pozycja_earth_icrs[1], pozycja_earth_icrs[2] )
vel_icrs = CartesianDifferential(v_earth_icrs[0], v_earth_icrs[1], v_earth_icrs[2])

#=== Transformacja do układu ITRS ===
pos_with_vel = pos_icrs.with_differentials(vel_icrs)
coordinates_icrs = SkyCoord(pos_with_vel, frame=ICRS(), obstime=time)
coordinates_itrs = coordinates_icrs.transform_to(ITRS(obstime=time))

pos_itrs = coordinates_itrs.cartesian.xyz
vel_itrs = coordinates_itrs.velocity.d_xyz
print()
print("ITRS position (m):", pos_itrs.to(u.m))
print("ITRS velocity (m/s):", vel_itrs.to(u.m/u.s))


#=== Transformacja do układu galaktycznego ===
coordinates_galactic = coordinates_itrs.transform_to(Galactic())
pos_galactic = coordinates_galactic.cartesian.xyz
vel_galactic = coordinates_galactic.velocity.d_xyz

print()
print("Galactic position (m):", pos_galactic.to(u.m))
print("Galactic velocity (m/s):", vel_galactic.to(u.m/u.s))
print()
vel_galactic_norm = np.sqrt(vel_galactic[0]**2 + vel_galactic[1]**2 + vel_galactic[2]**2)
print("Galactic velocity norm (m/s):", vel_galactic_norm.to(u.km/u.s))


#=== Transformacja do układu galaktycznego centrycznego ===
coordinates_galactocentric = coordinates_itrs.transform_to(Galactocentric())
pos_galactocentric = coordinates_galactocentric.cartesian.xyz
vel_galactocentric = coordinates_galactocentric.velocity.d_xyz

print()
print("Galactocentric position (m):", pos_galactocentric.to(u.m))
print("Galactocentric velocity (m/s):", vel_galactocentric.to(u.m/u.s))
print()
vel_galactocentric_norm = np.sqrt(vel_galactocentric[0]**2 + vel_galactocentric[1]**2 + vel_galactocentric[2]**2)
print("Galactocentric velocity norm (m/s):", vel_galactocentric_norm.to(u.km/u.s))