import math
import numpy as np

def mjd_to_jd(mjd):
    """Zamiana Modified Julian Date -> Julian Date."""
    return mjd + 2400000.5

def gmst_from_jd(jd):
    """
    Przybliżony GMST (Greenwich Mean Sidereal Time) w radianach.
    Wystarczy do obrotu z układu inercjalnego do Ziemi-obrotowego (ECEF/ITRF) bez polary motion.
    Wzór IAU uproszczony.
    """
    T = (jd - 2451545.0) / 36525.0  # wieki juliańskie od J2000.0
    # GMST w sekundach czasu
    gmst_sec = (67310.54841
                + (876600.0 * 3600 + 8640184.812866) * T
                + 0.093104 * T**2
                - 6.2e-6 * T**3)
    # na stopnie
    gmst_deg = (gmst_sec / 240.0) % 360.0
    return math.radians(gmst_deg)

def earth_orbital_velocity_icrs(mjd):
    """
    Bardzo uproszczona prędkość orbitalna Ziemi wokół Słońca w układzie równikowym (ICRS ~ J2000).
    Używamy klasycznych przybliżeń Meeusa: pozycja Słońca na ekliptyce i przesunięcie o 90° do kierunku ruchu.
    Dokładność: rząd 1 km/s, czyli OK do zgrubnej analizy.
    Zwraca wektor [vx, vy, vz] w km/s.
    """
    jd = mjd_to_jd(mjd)
    T = (jd - 2451545.0) / 36525.0

    # Średnia anomalia Słońca (w stopniach)
    M = 357.52911 + 35999.05029 * T - 0.0001537 * T**2
    M = math.radians(M % 360.0)

    # Długość ekliptyczna Słońca (pozorna) – Meeus
    L0 = 280.46646 + 36000.76983 * T + 0.0003032 * T**2
    L0 = L0 % 360.0
    # Równanie środka
    C = (1.914602 - 0.004817 * T - 0.000014 * T**2) * math.sin(M) \
        + (0.019993 - 0.000101 * T) * math.sin(2 * M) \
        + 0.000289 * math.sin(3 * M)
    # ekliptyczna długość Słońca
    lam_sun = math.radians((L0 + C) % 360.0)

    # Dla Ziemi wektor pozycji jest przeciwny, a prędkość orbitalna jest styczna,
    # czyli przesunięta o +90° w kierunku rosnącej długości
    lam_earth = (lam_sun + math.pi) % (2 * math.pi)
    tangential_dir = lam_earth + math.pi / 2.0

    # Średnia prędkość orbitalna ~29.78 km/s
    v_orb = 29.78  # km/s

    # W układzie ekliptycznym (x_ecl, y_ecl, z_ecl=0)
    vx_ecl = v_orb * math.cos(tangential_dir)
    vy_ecl = v_orb * math.sin(tangential_dir)
    vz_ecl = 0.0

    # Obrót ekliptyka -> równik (ICRS) przez nachylenie osi
    eps = math.radians(23.439291 - 0.0130042 * T)  # nachylenie ekliptyki
    vx_eq = vx_ecl
    vy_eq = vy_ecl * math.cos(eps)
    vz_eq = vy_ecl * math.sin(eps)

    return np.array([vx_eq, vy_eq, vz_eq])  # km/s

def sun_velocity_in_icrs_from_galactic():
    """
    Prędkość Słońca względem środka Galaktyki w układzie ICRS.
    Robimy to tak:
    1. w galaktycznym (U,V,W): U do środka Galaktyki, V w kierunku rotacji, W do NGP.
       bierzemy LSR + ruch własny Słońca.
    2. Mnożymy przez macierz obrotu gal->ICRS (jak w astropy).
    Wynik w km/s.
    """
    # Parametry (można dostroić):
    V0 = 238.0  # km/s, prędkość LSR (możesz zmienić na 220–240)
    # ruch własny Słońca wg Schönrich+2010
    U_sun = 11.1   # km/s (do środka Galaktyki)
    V_sun = 12.24  # km/s (w kierunku rotacji)
    W_sun = 7.25   # km/s (do NGP)

    # razem w galaktycznym
    v_gal = np.array([
        U_sun,
        V0 + V_sun,
        W_sun
    ])  # km/s

    # macierz z układu galaktycznego do ICRS (J2000)
    R = np.array([
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [ 0.4941094279, -0.4448296300,  0.7469822445],
        [-0.8676661490, -0.1980763734,  0.4559837762]
    ])

    v_icrs = R @ v_gal
    return v_icrs  # km/s

def earth_velocity_itrf_from_mjd(mjd):
    """
    Główna funkcja o którą chodziło:
    - bierze MJD
    - liczy prędkość Słońca w Galaktyce (w ICRS)
    - dodaje prędkość orbitalną Ziemi wokół Słońca (w ICRS)
    - obraca do układu związanego z Ziemią (ITRF/ECEF) przez GMST
      (pomijamy nutację, precesję biegunów, x_p, y_p)
    - zwraca wektor [vx, vy, vz] w m/s w osi ITRF:
        X: Greenwich, równik
        Y: 90°E, równik
        Z: biegun północny
    """
    jd = mjd_to_jd(mjd)

    # 1) prędkość "galaktyczna" Słońca
    v_sun_icrs = sun_velocity_in_icrs_from_galactic()  # km/s

    # 2) prędkość orbitalna Ziemi
    v_earth_orb_icrs = earth_orbital_velocity_icrs(mjd)  # km/s

    # 3) suma: prędkość środka Ziemi względem Galaktyki
    v_tot_icrs = v_sun_icrs + v_earth_orb_icrs  # km/s

    # 4) obrót ICRS -> ITRF przez GMST (tylko Z)
    gmst = gmst_from_jd(jd)
    cosg = math.cos(gmst)
    sing = math.sin(gmst)
    R3 = np.array([
        [ cosg,  sing, 0.0],
        [-sing,  cosg, 0.0],
        [ 0.0,   0.0,  1.0]
    ])
    v_itrf = R3 @ v_tot_icrs  # km/s

    # na m/s
    return v_itrf * 1000.0


from astropy.time import Time
from astropy.coordinates import SkyCoord, Galactocentric
import astropy.units as u
import numpy as np

def earth_velocity_itrf_from_mjd_astropy(mjd):
    """
    Wersja 'luxury': wymaga astropy i jej efemeryd.
    1. bierzemy prędkość Ziemi (geocentrum) względem barycentrum z astropy
    2. dodajemy prędkość Słońca względem Galaktyki (jak wcześniej)
    3. obracamy do ITRF przez GMST
    """
    t = Time(mjd, format='mjd', scale='utc')
    # prędkość geocentrum względem Słońca/barycentrum:
    # Astropy: get_body_barycentric_posvel('earth', t)[1]
    from astropy.coordinates import get_body_barycentric_posvel
    _, v_earth = get_body_barycentric_posvel('earth', t)
    v_earth_icrs = v_earth.xyz.to(u.km/u.s).value  # km/s

    # nasza prędkość Słońca w ICRS
    v_sun_icrs = sun_velocity_in_icrs_from_galactic()  # km/s

    v_tot_icrs = v_sun_icrs + v_earth_icrs  # km/s

    jd = mjd_to_jd(mjd)
    gmst = gmst_from_jd(jd)
    cosg = math.cos(gmst); sing = math.sin(gmst)
    R3 = np.array([[ cosg,  sing, 0],
                   [-sing,  cosg, 0],
                   [ 0,     0,    1]])
    v_itrf = R3 @ v_tot_icrs
    return v_itrf * 1000.0  # m/s

if __name__ == "__main__":
    mjd_test = 59000.5
    v_itrf = earth_velocity_itrf_from_mjd(mjd_test)
    v_itrf_astropy = earth_velocity_itrf_from_mjd_astropy(mjd_test)
    print("Earth velocity in ITRF at MJD {:.1f}:".format(mjd_test))
    print("  Custom code:   vx={:.2f} m/s, vy={:.2f} m/s, vz={:.2f} m/s".format(v_itrf[0], v_itrf[1], v_itrf[2]))
    print("  Astropy code:  vx={:.2f} m/s, vy={:.2f} m/s, vz={:.2f} m/s".format(v_itrf_astropy[0], v_itrf_astropy[1], v_itrf_astropy[2]))
    print("Difference (custom - astropy):")
    print("  dvx={:.2f} m/s, dvy={:.2f} m/s, dvz={:.2f} m/s".format(v_itrf[0]-v_itrf_astropy[0], v_itrf[1]-v_itrf_astropy[1], v_itrf[2]-v_itrf_astropy[2]))