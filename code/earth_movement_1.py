from astropy.coordinates import get_body_barycentric_posvel, ICRS, get_body, SkyCoord, ITRS, Galactic,Galactocentric, CartesianRepresentation, CartesianDifferential
from astropy.time import Time
import astropy.units as u
import numpy as np

current_jd = Time.now()
#print(current_jd)
time = Time(current_jd, format="jd")
#time = Time(59000.0, format="mjd")
#print(time)
mjd = time-2400000.5
#print("MJD:", mjd.value)
un_coords = get_body('earth', time)



#=== Koordynaty i prędkości Ziemi i Słońca w układzie barycentrycznym ICRS ===
earth_pv = get_body_barycentric_posvel("earth", time)
sun_pv = get_body_barycentric_posvel("sun", time)

#print()
pozycja_earth_icrs = earth_pv[0].xyz.to(u.m)
#print("Ziemia - pozycja w ICRS (m):", pozycja_earth_icrs)
pozycja_sun_icrs = sun_pv[0].xyz.to(u.m)

v_earth_icrs = earth_pv[1].xyz.to(u.m/u.s)
#print("Ziemia - prędkość w ICRS (m/s):", v_earth_icrs)
v_sun_icrs = sun_pv[1].xyz.to(u.m/u.s)

pos_icrs_earth = CartesianRepresentation(pozycja_earth_icrs[0], pozycja_earth_icrs[1], pozycja_earth_icrs[2] )
pos_icrs_sun = CartesianRepresentation(pozycja_sun_icrs[0], pozycja_sun_icrs[1], pozycja_sun_icrs[2] )

vel_icrs_earth = CartesianDifferential(v_earth_icrs[0], v_earth_icrs[1], v_earth_icrs[2])
vel_icrs_sun = CartesianDifferential(v_sun_icrs[0], v_sun_icrs[1], v_sun_icrs[2])



#=== Transformacja do układu ITRS ===
pos_with_vel_earth = pos_icrs_earth.with_differentials(vel_icrs_earth)
coordinates_icrs_earth = SkyCoord(pos_with_vel_earth, frame=ICRS(), obstime=time)
coordinates_itrs_earth = coordinates_icrs_earth.transform_to(ITRS(obstime=time))

pos_itrs = coordinates_itrs_earth.cartesian.xyz
vel_itrs_earth = coordinates_itrs_earth.velocity.d_xyz

coordinates_icrs_sun = SkyCoord(pos_icrs_sun.with_differentials(vel_icrs_sun), frame=ICRS(), obstime=time)
coordinates_itrs_sun = coordinates_icrs_sun.transform_to(ITRS(obstime=time))
vel_itrs_sun = coordinates_itrs_sun.velocity.d_xyz


'''print('=============================================')
print("Earth velocity in ITRS (Earth' center) at MJD " + str(mjd.value))
#print()
#print("ITRS position (m):", pos_itrs.to(u.m))
print("ITRS velocity (m/s):", vel_itrs_earth.to(u.m/u.s))
print()
print("ITRS velocity norm (km/s):", np.sqrt(vel_itrs_earth[0]**2 + vel_itrs_earth[1]**2 + vel_itrs_earth[2]**2).to(u.km/u.s))

print()
sum_vel_ITRS = vel_itrs_earth + vel_itrs_sun
print("sum_vel_ITRS =", sum_vel_ITRS)
print("sum_vel_ITRS_norm = ", np.sqrt(sum_vel_ITRS[0]**2 + sum_vel_ITRS[1]**2 + sum_vel_ITRS[2]**2).to(u.km/u.s))

print()
print()
#=== Transformacja do układu galaktycznego ===
coordinates_galactic_earth = coordinates_itrs_earth.transform_to(Galactic())
pos_galactic_earth = coordinates_galactic_earth.cartesian.xyz
vel_galactic_earth = coordinates_galactic_earth.velocity.d_xyz

coordinates_galactic_sun = coordinates_itrs_sun.transform_to(Galactic())
vel_galactic_sun = coordinates_galactic_sun.velocity.d_xyz


print('=============================================')
print("Earth velocity in Galactic (Sun's barycenter) at MJD " + str(mjd.value))
#print("Galactic position (m):", pos_galactic.to(u.m))
print("Galactic velocity (m/s):", vel_galactic_earth.to(u.m/u.s))
#print()
vel_galactic_norm = np.sqrt(vel_galactic_earth[0]**2 + vel_galactic_earth[1]**2 + vel_galactic_earth[2]**2)
print()
print("Galactic velocity norm (m/s):", vel_galactic_norm.to(u.km/u.s))

print()
sum_vel = vel_galactic_earth + vel_galactic_sun
print("sum_vel =", sum_vel)
print("sum_vel_norm = ", np.sqrt(sum_vel[0]**2 + sum_vel[1]**2 + sum_vel[2]**2).to(u.km/u.s))
'''

print()
print()
#=== Transformacja do układu galaktycznego centrycznego ===
coordinates_galactocentric_earth = coordinates_itrs_earth.transform_to(Galactocentric())
pos_galactocentric_earth = coordinates_galactocentric_earth.cartesian.xyz
vel_galactocentric_earth = coordinates_galactocentric_earth.velocity.d_xyz

vel_galactocentric_sun = coordinates_itrs_sun.transform_to(Galactocentric()).velocity.d_xyz


print('=============================================')
print("Earth velocity in Galactocentric (galaxy center) at MJD " + str(mjd.value))
#print("Galactocentric position (m):", pos_galactocentric.to(u.m))
print("Galactocentric velocity (m/s):", vel_galactocentric_earth.to(u.m/u.s))
print()
vel_galactocentric_norm = np.sqrt(vel_galactocentric_earth[0]**2 + vel_galactocentric_earth[1]**2 + vel_galactocentric_earth[2]**2)
print("Galactocentric velocity norm (km/s):", vel_galactocentric_norm.to(u.km/u.s))

'''print()
sum_vel = vel_galactocentric_sun + vel_galactocentric_earth
print("sum_vel =", sum_vel)
print("sum_vel_norm = ", np.sqrt(sum_vel[0]**2 + sum_vel[1]**2 + sum_vel[2]**2).to(u.km/u.s))
'''

print()
print()
print("=============================================")
print("ITRS sum transform to Galactocentric at MJD " + str(mjd.value))

sum_pos_icrs = pos_icrs_earth + pos_icrs_sun
sum_vel_icrs = CartesianDifferential(v_earth_icrs[0], v_earth_icrs[1], v_earth_icrs[2]) + CartesianDifferential(v_sun_icrs[0], v_sun_icrs[1], v_sun_icrs[2])

#print("sum_pos_icrs =", sum_pos_icrs)
#print("sum_vel_icrs =", sum_vel_icrs)

sum_coordinates = SkyCoord(sum_pos_icrs.with_differentials(vel_icrs_earth), frame=ICRS(), obstime=time)

sum_coordinates_itrs = sum_coordinates.transform_to(ITRS(obstime=time))
sum_coordinates_galactocentric = sum_coordinates_itrs.transform_to(Galactocentric())
vel_sum_galactocentric = sum_coordinates_galactocentric.velocity.d_xyz

print("Galactocentric velocity sum (m/s):", vel_sum_galactocentric.to(u.m/u.s))
print()
print("Galactocentric velocity sum norm (km/s):", np.sqrt(vel_sum_galactocentric[0]**2 + vel_sum_galactocentric[1]**2 + vel_sum_galactocentric[2]**2).to(u.km/u.s))