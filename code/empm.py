import math
import numpy as np

def mjd_to_jd(mjd):
    """
    Zamiana Modified Julian Date (MJD) -> Julian Date (JD).

    Skąd 2400000.5?
    - Definicje:
        JD  = liczba dni (z częścią ułamkową) od południa UT 1 stycznia 4713 p.n.e.
        MJD = JD - 2400000.5  (przesunięty tak, by numery były krótsze i zaczynały się o północy, a nie w południe)
    - Zatem JD = MJD + 2400000.5

    Źródła: standard astronomiczny, np. Meeus „Astronomical Algorithms”.
    """
    return mjd + 2400000.5


def gmst_from_jd(jd):
    """
    Przybliżony GMST (Greenwich Mean Sidereal Time) w radianach.

    Co liczymy?
    - GMST to "średni czas gwiazdowy" w Greenwich. Umożliwia obrót wektora z układu
      inercjalnego (GCRS/ICRS) do obracającego się z Ziemią (ECEF/ITRF), zakładając
      pominięcie kilku subtelności (nutacji, ruchu bieguna, UT1-UTC).

    Skąd wzór i stałe?
    - Używany tu wielomian (w sekundach czasu) jest klasycznym przybliżeniem
      spotykanym m.in. u Vallado i w materiałach IAU (blisko IAU 1982/2000).
      Zapis:
        T = (JD - 2451545.0)/36525 (wieki juliańskie od epoki J2000.0)
        GMST_sec = 67310.54841
                   + (876600 * 3600 + 8640184.812866) * T
                   + 0.093104 * T^2
                   - 6.2e-6 * T^3
      gdzie:
        - 67310.54841 s  to GMST w J2000.0,
        - 876600*3600 + 8640184.812866 = 1 dzień gwiazdowy * liczba sekund/obr
          w formalnym rozwinięciu na T (tu bez dokładnej separacji),
        - współczynniki T^2 i T^3 to drobne poprawki sekularne.

    Konwersje:
    - sekundy czasu -> stopnie: 1° = 240 s   (360° = 24h -> 1h = 15°, 1s czasu = 15 arcsec = 1/240°)
    - na końcu zwracamy wartość w radianach, zredukowaną mod 2π (poprzez % 360°).

    Uwaga:
    - To „mean” ST (bez nutacji, polarmotion, bez ERA/UT1), wystarczające do obrotu o R3(GMST)
      dla zastosowań, gdzie błędy rzędu ~0.1° nie są krytyczne.

    Źródła: Vallado „Fundamentals of Astrodynamics and Applications”, IAU 1982/2000,
            IERS Conventions 2010 (dokładniejsze procedury).
    """
    T = (jd - 2451545.0) / 36525.0  # wieki juliańskie od J2000.0

    # GMST w sekundach czasu (patrz opis wyżej)
    gmst_sec = (67310.54841
                + (876600.0 * 3600 + 8640184.812866) * T
                + 0.093104 * T**2
                - 6.2e-6 * T**3)

    # sekundy -> stopnie; redukcja do [0,360)
    gmst_deg = (gmst_sec / 240.0) % 360.0

    # stopnie -> radiany
    return math.radians(gmst_deg)


def earth_orbital_velocity_icrs(mjd):
    """
    Uproszczona prędkość orbitalna Ziemi wokół Słońca w układzie równikowym (ICRS ~ J2000).

    Idea:
    - Używamy klasycznej, *analitycznej* (nie-efemerydalnej) aproksymacji wg Meeusa:
      najpierw pozycja Słońca na ekliptyce (długość ekliptyczna λ_☉), potem:
        * pozycja Ziemi jest przeciwnie skierowana (λ_⊕ = λ_☉ + 180°),
        * prędkość orbitalna jest styczna do orbity, czyli kierunek = λ_⊕ + 90°,
      przyjmujemy stałą wartość modułu prędkości ~29.78 km/s (średnia prędkość orbitalna).

    Dokładność:
    - Rząd 1 km/s (wystarczające do sanity-check i zgrubnej analizy).
      Dokładniej robi to wersja z astropy.

    Kroki wzorów Meeusa:
    1) Epoka:
       T = (JD - 2451545.0)/36525
    2) Średnia anomalia Słońca M (w stopniach):
       M = 357.52911 + 35999.05029*T - 0.0001537*T^2
       (określa położenie na orbicie eliptycznej względem peryhelium)
    3) Średnia długość L0 (pozorna) Słońca (stopnie):
       L0 = 280.46646 + 36000.76983*T + 0.0003032*T^2
    4) Równanie środka C (korekta eliptyczności):
       C = (1.914602 - 0.004817*T - 0.000014*T^2)*sin M
           + (0.019993 - 0.000101*T)*sin(2M)
           + 0.000289*sin(3M)
    5) Długość ekliptyczna Słońca: λ_☉ = L0 + C
       Z tego wyprowadzamy kierunek styczny do orbity Ziemi.

    Transformacje:
    - Ekliptyczny -> równikowy: obrót o nachylenie ekliptyki ε(T) wokół osi x:
      [x_eq, y_eq, z_eq]^T = R1(+ε) * [x_ecl, y_ecl, z_ecl]^T
      Wektor prędkości w płaszczyźnie ekliptyki ma z_ecl = 0.

    Źródło: J. Meeus, „Astronomical Algorithms”.
    """
    jd = mjd_to_jd(mjd)
    T = (jd - 2451545.0) / 36525.0

    # 2) Średnia anomalia Słońca (stopnie -> rad)
    M = 357.52911 + 35999.05029 * T - 0.0001537 * T**2
    M = math.radians(M % 360.0)

    # 3) Średnia długość Słońca L0 (stopnie), redukcja mod 360°
    L0 = 280.46646 + 36000.76983 * T + 0.0003032 * T**2
    L0 = L0 % 360.0

    # 4) Równanie środka C (stopnie -> bezpośrednio na sin M, M w radianach)
    C = (1.914602 - 0.004817 * T - 0.000014 * T**2) * math.sin(M) \
        + (0.019993 - 0.000101 * T) * math.sin(2 * M) \
        + 0.000289 * math.sin(3 * M)

    # 5) Długość ekliptyczna Słońca λ_☉ (rad)
    lam_sun = math.radians((L0 + C) % 360.0)

    # Pozycja Ziemi jest przeciwnie skierowana: λ_⊕ = λ_☉ + π
    lam_earth = (lam_sun + math.pi) % (2 * math.pi)

    # Kierunek wektora prędkości orbitalnej to styczna: +90° do kierunku promienia
    tangential_dir = lam_earth + math.pi / 2.0

    # Moduł prędkości orbitalnej ~ 29.78 km/s (średnia po orbicie)
    v_orb = 29.78  # km/s

    # W układzie ekliptycznym (z_ecl = 0, bo prędkość leży w płaszczyźnie ekliptyki)
    vx_ecl = v_orb * math.cos(tangential_dir)
    vy_ecl = v_orb * math.sin(tangential_dir)
    vz_ecl = 0.0

    # Obrót ekliptyka -> równik (ICRS) o nachylenie ekliptyki ε(T).
    # Uwaga: dokładniejszy model ε znajdziesz u Meeusa; tu prosty liniowy przybliżony:
    eps = math.radians(23.439291 - 0.0130042 * T)  # w stopniach -> rad

    # Macierz R1(+ε) zastosowana "ręcznie" (bo z_ecl = 0, upraszcza się):
    # x_eq = x_ecl
    # y_eq = y_ecl * cos ε
    # z_eq = y_ecl * sin ε
    vx_eq = vx_ecl
    vy_eq = vy_ecl * math.cos(eps)
    vz_eq = vy_ecl * math.sin(eps)

    return np.array([vx_eq, vy_eq, vz_eq])  # km/s

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

def earth_velocity_itrf_from_mjd(mjd):
    """
    Główna funkcja:
      - wejście: MJD,
      - wyjście: wektor prędkości Ziemi (geocentrum) względem Galaktyki
                 wyrażony w osi ITRF/ECEF (X do Greenwich na równiku,
                 Y do 90°E na równiku, Z do bieguna północnego), w m/s.

    Składanie wektora prędkości:
      v_total(ICRS) = v_sun(ICRS) + v_earth_orbit(ICRS)

    Transformacja układów:
      ICRS (inercjalny) -> ITRF (obracający się z Ziemią)
      - tu stosujemy jedynie obrót wokół osi Z o kąt GMST:
            v_ITRF = R3(GMST) * v_ICRS
        gdzie:
            R3(G) = [[ cosG,  sinG, 0],
                     [-sinG,  cosG, 0],
                     [   0,     0,  1]]
        To odpowiada rotacji „nieba” o GMST, czyli rzutuje wektor inercjalny
        na chwilowe osie związane z Ziemią (pomijamy precesję-nutację, polarmotion, i UT1-UTC).

    Uwaga:
      - Nie dodajemy (ω × r) punktu na powierzchni – zwracamy prędkość „geocentrum” w ITRF.
        Jeśli chcesz prędkość laboratorium, wtedy do wyniku w ITRF dodaj (ω × r) w ITRF.
    """
    jd = mjd_to_jd(mjd)

    # (1) Składowa galaktyczna Słońca (w ICRS)
    v_sun_icrs = sun_velocity_in_icrs_from_galactic()  # km/s

    # (2) Składowa orbitalna Ziemi wokół Słońca (w ICRS)
    v_earth_orb_icrs = earth_orbital_velocity_icrs(mjd)  # km/s

    # (3) Suma w ICRS: prędkość geocentrum względem Galaktyki
    v_tot_icrs = v_sun_icrs + v_earth_orb_icrs  # km/s

    # (4) Obrót ICRS -> ITRF przez GMST (obrót o kąt GMST wokół Z)
    gmst = gmst_from_jd(jd)
    cosg = math.cos(gmst)
    sing = math.sin(gmst)
    R3 = np.array([
        [ cosg,  sing, 0.0],
        [-sing,  cosg, 0.0],
        [ 0.0,   0.0,  1.0]
    ])
    v_itrf = R3 @ v_tot_icrs  # km/s

    # km/s -> m/s
    return v_itrf * 1000.0


from astropy.time import Time
from astropy.coordinates import SkyCoord, Galactocentric
import astropy.units as u
import numpy as np

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
    _, v_earth = get_body_barycentric_posvel('earth', t)
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
    return v_itrf * 1000.0  # m/s

#mjd_test = Time.now()
#time = Time(mjd_test, format="jd")
print("=============================================")
if __name__ == "__main__":
    mjd_test = 59000.0
    #mjd_test = 61615.37540498899
    v_itrf = earth_velocity_itrf_from_mjd(mjd_test)
    v_itrf_astropy = earth_velocity_itrf_from_mjd_astropy(mjd_test)
    print("Earth velocity in ITRF at MJD {:.1f}:".format(mjd_test))
    print("  Custom code:   vx={:.2f} m/s, vy={:.2f} m/s, vz={:.2f} m/s".format(v_itrf[0], v_itrf[1], v_itrf[2]))
    print("  Astropy code:  vx={:.2f} m/s, vy={:.2f} m/s, vz={:.2f} m/s".format(v_itrf_astropy[0], v_itrf_astropy[1], v_itrf_astropy[2]))
    print("Custom Norm: v={:.2f} km/s".format((np.sqrt(v_itrf[0]**2 + v_itrf[1]**2 + v_itrf[2]**2))/1000.0))
    print("Astropy Norm: v={:.2f} km/s".format((np.sqrt(v_itrf_astropy[0]**2 + v_itrf_astropy[1]**2 + v_itrf_astropy[2]**2))/1000.0))
    #print("Difference (custom - astropy):")
    #print("  dvx={:.2f} m/s, dvy={:.2f} m/s, dvz={:.2f} m/s".format(v_itrf[0]-v_itrf_astropy[0], v_itrf[1]-v_itrf_astropy[1], v_itrf[2]-v_itrf_astropy[2]))