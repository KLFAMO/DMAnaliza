from astropy.time import Time
from astropy.coordinates import SkyCoord, Galactocentric
import astropy.units as u
from empm import gmst_from_jd, mjd_to_jd
import numpy as np
import math


def sun_velocity_in_icrs_from_galactic():
    """
    Prędkość Słońca względem środka Galaktyki wyrażona w ICRS.

    Kroki:
    1) W układzie galaktycznym używamy dekompozycji (U, V, W):
       - U: dodatnie w stronę środka Galaktyki,
       - V: dodatnie w kierunku rotacji galaktycznej (wzdłuż +l=90°),
       - W: dodatnie w stronę północnego bieguna galaktycznego (NGP).
       Przyjmujemy: prędkość LSR (Local Standard of Rest) ~ V0 oraz
       ruch własny Słońca (tzw. peculiar motion) wg Schönrich+2010.

       Typowe wartości:
         V0 ≈ 220–240 km/s  (tu: 238.0 km/s)
         (U, V, W)_☉ ≈ (11.1, 12.24, 7.25) km/s

       Prędkość Słońca względem Galaktyki w układzie galaktycznym:
         v_gal = [ U_☉,
                   V0 + V_☉,
                   W_☉ ]

    2) Obrót do ICRS macierzą gal->ICRS.
       Używamy stałej macierzy (definicja osi galaktycznych względem J2000/ICRS),
       zgodnej z implementacją w Astropy (wywodzi się z IAU 1958/2000 z korektami).

    Wynik: v_icrs w km/s.

    Źródła:
      - Schönrich, Binney, Dehnen (2010), MNRAS 403, 1829 (peculiar motion Słońca)
      - Definicja osi galaktycznych i macierze w dokumentacji Astropy / IAU
    """
    # (1) Parametry modelu prędkości
    V0 = 238.0  # km/s, przyjęta prędkość rotacji LSR (można podmienić na 220–240)
    U_sun = 11.1   # km/s (do środka Galaktyki)
    V_sun = 12.24  # km/s (zgodnie z kierunkiem rotacji)
    W_sun = 7.25   # km/s (do NGP)

    # (1) Wektor w układzie galaktycznym (U, V, W)
    v_gal = np.array([
        U_sun,
        V0 + V_sun,
        W_sun
    ])  # km/s

    # (2) Macierz obrotu: Galactic -> ICRS (J2000).
    # Ta macierz jest standardowa; numery znajdziesz w astropy/erfa/iau.
    R = np.array([
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [ 0.4941094279, -0.4448296300,  0.7469822445],
        [-0.8676661490, -0.1980763734,  0.4559837762]
    ])

    # Mnożenie macierzy przez wektor: v_ICRS = R * v_gal
    v_icrs = R @ v_gal
    return v_icrs  # km/s



def earth_velocity_itrf_from_mjd_astropy(mjd):
    """
    Wersja z wykorzystaniem astropy + jej efemeryd.
    Kroki:
      1) Z astropy bierzemy prędkość Ziemi (geocentrum) względem barycentrum Układu Słonecznego:
         get_body_barycentric_posvel('earth', t) zwraca (pozycja, prędkość) w GCRS/ICRS.
      2) Dodajemy naszą składową galaktyczną Słońca (w ICRS) — jak w funkcji powyżej.
      3) Obracamy R3(GMST) do ITRF (jak wyżej).

    Uwaga:
      - Ta wersja lepiej oddaje *zmienność* prędkości orbitalnej (prędkość i kierunek),
        bo korzysta z efemeryd, a nie stałej wartości 29.78 km/s.
      - Składowa galaktyczna nadal jest modelem parametrycznym (V0, U,V,W).
    """
    t = Time(mjd, format='mjd', scale='utc')

    # (1) prędkość geocentrum względem barycentrum (ICRS) z efemeryd
    from astropy.coordinates import get_body_barycentric_posvel
    pos_earth, v_earth = get_body_barycentric_posvel('earth', t)
    pos_earth = pos_earth.xyz.to(u.km).value  # km, wektor 3D
    v_earth_icrs = v_earth.xyz.to(u.km/u.s).value  # km/s, wektor 3D

    # (2) składowa Słońca w ICRS (jak wyżej)
    v_sun_icrs = sun_velocity_in_icrs_from_galactic()  # km/s

    # (3) suma (ICRS)
    v_tot_icrs = v_sun_icrs + v_earth_icrs  # km/s

    # (4) ICRS -> ITRF przez GMST (obrót R3)
    jd = mjd_to_jd(mjd)
    gmst = gmst_from_jd(jd)
    cosg = math.cos(gmst); sing = math.sin(gmst)
    R3 = np.array([[ cosg,  sing, 0],
                   [-sing,  cosg, 0],
                   [ 0,     0,    1]])
    v_itrf = R3 @ v_tot_icrs
    return pos_earth, v_itrf # km/s

