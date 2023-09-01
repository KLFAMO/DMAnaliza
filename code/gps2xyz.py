"""
Author: Joséphine Strübing-Tardy
"""

import math
import numpy as np
import re

def dms2dd(lat):
    # Converts GPS coordinates from hours, minutes, second to decimal degrees format    
    deg, minutes, seconds, direction =  re.split('[°\'"]', lat)
    newlat = (float(deg) + float(minutes)/60 + float(seconds)/(60*60)) * (-1 if direction in ['W', 'S'] else 1)
    return newlat

def dd2dms(dd):
    # Converts coordinates from decimal degrees to hours, minute, seconds
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
    Input: Longitude (decimal degrees), 
           Latitude (decimal degrees), 
           Height (m)
    Output: X, Y, Z in the ECEF
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

# checking for Torun 

long = dms2dd("18° 35' 53.30' E")
lat = dms2dd("53° 00' 49.50' N")
print("GPS coordinates in decimal degrees:",long, lat)

h = 65 # above sea-level, in m 

x, y, z = gps2ecef(long, lat, h)
print("X:", x)
print("Y:", y)
print("Z:", z)

print('\nPTB')
print(gps2ecef(10.459046, 52.295024, 50))

print('\nKRISS')
print(gps2ecef(127.372131, 36.387554, 50))

print('\nNMIJ')
print(gps2ecef(140.099934, 36.112252, 50))

print('\nSYRTE')
print(gps2ecef(2.334838, 48.835866, 50))
