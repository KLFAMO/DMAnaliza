from astropy.coordinates import ICRS, SkyCoord, ITRS, Galactocentric, CartesianRepresentation, CartesianDifferential, get_body_barycentric_posvel, galactocentric_frame_defaults
from astropy.time import Time
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

OMEGA_EARTH = 7.2921150e-5 / u.s

velocity_cache = {} # zmienna globalna która pamięta obliczone prędkości dla całkowitych części MJD, aby przyspieszyć kolejne wywołania funkcji earth_velocity_fast() dla tych samych całkowitych części MJD.
t = Time(Time.now(), format="mjd")
omega_vec = [0, 0, 7.2921150e-5] /u.s # Earth's angular velocity

def ITRS_to_ICRS(t, point_on_earth): 
    #=== Koordynaty miejsca na Ziemi w układzie ITRS ===
    t = Time(t, format="mjd")
    on_earth_itrs = SkyCoord(point_on_earth, frame=ITRS(obstime=t))
    
    #=== Sprawdzanie czy predkość obrotowa Ziemi jest dobrze liczona ===
    #on_earth_velocity = CartesianDifferential(np.cross(omega_vec.to(1/u.s).value, on_earth_itrs.cartesian.xyz.to(u.m).value)*(u.m/u.s))
    #print()

    #=== Transformacja położenia do układu ICRS ===
    on_earth_icrs = on_earth_itrs.transform_to(ICRS())

    #=== Koordynaty i prędkości Ziemi w układzie ICRS ===
    pos_e, vel_e = get_body_barycentric_posvel("earth", t)
    earth_icrs = SkyCoord(CartesianRepresentation(pos_e.xyz).with_differentials(CartesianDifferential(vel_e.xyz)),frame=ICRS(),obstime=t)
    #cr = CartesianRepresentation(pos_e.xyz).with_differentials(CartesianDifferential(vel_e.xyz))
    #earth_icrs = SkyCoord(cr, frame=ICRS(), obstime=t)

    #=== Obliczanie prędkości obrotowej w układzie ICRS ===
    r_icrs = (on_earth_icrs.cartesian.xyz - earth_icrs.cartesian.xyz).to(u.m)
    v_rot_icrs = np.cross(omega_vec, r_icrs).to(u.m/u.s)

    #=== Całkowita prędkość miejsca na Ziemi w układzie ICRS ===
    v_earth_icrs = vel_e.xyz.to(u.m/u.s)
    v_total_icrs = v_earth_icrs + v_rot_icrs


    return on_earth_icrs, r_icrs, v_total_icrs


def ITRS_to_Galactocentric(t, point_on_earth):
    # --- bierzemy wynik z pierwszej funkcji ---
    t = Time(t, format="mjd")
    on_earth_icrs, r_icrs, v_total_icrs = ITRS_to_ICRS(t, point_on_earth)

    # --- tworzymy SkyCoord w ICRS z prędkością ---
    coord_icrs = SkyCoord(CartesianRepresentation(on_earth_icrs.cartesian.xyz).with_differentials(CartesianDifferential(v_total_icrs)),frame=ICRS(),obstime=t)

    # --- transformacja do Galactocentric ---
    coord_gal = coord_icrs.transform_to(Galactocentric())

    # --- wyciągamy wektory ---
    pos_gal = coord_gal.cartesian.xyz
    vel_gal = coord_gal.velocity.d_xyz
    

    return pos_gal, vel_gal


def rotation_z(angle_rad):      # macierz rotacji
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)

    return np.array([
        [c, -s, 0],
        [s,  c, 0],
        [0,  0, 1],
    ])


#=== Prędkość Ziemii ===

def earth_velocity(mjd, xyz_list=[6378137, 0, 0]): # xyz_list=[0, 0, 0]

    # Prędkość Ziemi względem halo, wyrażona w ITRS, wynik w m/s.
    t = Time(mjd, format="mjd")
    point_on_earth = CartesianRepresentation(*xyz_list) * u.m

    pos_gal, v_gal = ITRS_to_Galactocentric(t, point_on_earth)

    # Duży krok pomocniczy tylko do wyznaczenia kierunku wektora. To NIE jest ewolucja w czasie, tylko sztuczne przesunięcie przestrzenne.
    eps_time = 1000.0 * u.s

    coord0_gal = SkyCoord(CartesianRepresentation(pos_gal), frame=Galactocentric(), obstime=t)
    coord1_gal = SkyCoord(CartesianRepresentation(pos_gal + v_gal * eps_time), frame=Galactocentric(), obstime=t)

    coord0_itrs = coord0_gal.transform_to(ITRS(obstime=t))
    coord1_itrs = coord1_gal.transform_to(ITRS(obstime=t))

    v_itrs = (coord1_itrs.cartesian.xyz - coord0_itrs.cartesian.xyz) / eps_time

    return v_itrs.to(u.m / u.s)


def earth_velocity_fast(mjd):
    """
    Szybsza wersja earth_velocity().

    Dla MJD = N + f:
      1. liczy dokładnie earth_velocity(N) tylko raz,
      2. zapamiętuje wynik,
      3. dla kolejnych wywołań obraca wektor o czas f dni.
    """

    mjd0 = int(np.floor(mjd))
    frac = mjd - mjd0

    # Jeśli nie mamy jeszcze prędkości dla całkowitej części MJD, to ją liczymy i zapamiętujemy.
    # if mjd0 not in _velocity_cache: sprawdza czy mamy już prędkość dla całkowitej części MJD
    # w słowniku _velocity_cache. 
    # Jeśli nie, to obliczamy prędkość dla tej całkowitej części MJD i zapisujemy ją w słowniku pod kluczem mjd0.
    if mjd0 not in velocity_cache:
        velocity_cache[mjd0] = earth_velocity(mjd0, [0, 0, 0])

    # Pobieramy prędkość dla całkowitej części MJD z cache'u z tablicy _velocity_cache.
    v0 = velocity_cache[mjd0]

    # obliczenie kąta obrotu Ziemi w ciągu czasu frac dni 
    dt = frac * u.day
    # mnozymy u.day przez frac, aby uzyskać czas w dniach, a następnie mnożymy przez OMEGA_EARTH, aby uzyskać kąt obrotu Ziemi w ciągu tego czasu.
    # Wynik jest negowany, ponieważ chcemy obrócić wektor w kierunku przeciwnym do ruchu obrotowego Ziemi.
    angle = -(OMEGA_EARTH * dt).decompose().value

    # Tworzymy macierz obrotu wokół osi Z o obliczony kąt. 
    # Funkcja rotation_z(angle) zwraca macierz obrotu, która obraca wektor o kąt angle wokół osi Z.
    R = rotation_z(angle)

    # Obracamy wektor v0, mnożąc macierz obrotu R przez wektor v0. 
    # Wynik jest przekształcany z jednostek m/s na m/s, aby zachować spójność jednostek.
    v_rot = R @ v0.to_value(u.m / u.s)

    return v_rot * (u.m / u.s)


def plot_velocity_compare(
    velocity_functions,
    mjd_start=60000,
    mjd_stop=60000.5,
    step=0.01,
):
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")

    mjd_values = np.arange(mjd_start, mjd_stop, step)

    linestyles = ["-", "--", "-.", ":"]
    markers = ["o", "s", "^", "D", "x", "+"]

    for i, velocity_function in enumerate(velocity_functions):
        xs = []
        ys = []
        zs = []

        for mjd in mjd_values:
            v = velocity_function(mjd).to(u.km / u.s)

            xs.append(v[0].value)
            ys.append(v[1].value)
            zs.append(v[2].value)

        function_name = velocity_function.__name__

        ax.plot(
            xs,
            ys,
            zs,
            linestyle=linestyles[i % len(linestyles)],
            marker=markers[i % len(markers)],
            markersize=4,
            linewidth=2,
            alpha=0.85,
            label=function_name,
        )

        # start
        ax.scatter(
            xs[0],
            ys[0],
            zs[0],
            marker="*",
            s=120,
            color="black",
        )

        # stop
        ax.scatter(
            xs[-1],
            ys[-1],
            zs[-1],
            marker="X",
            s=80,
            color="black",
        )


    # początek układu
    ax.scatter(0, 0, 0, marker="+", s=120, color="black")

    ax.set_xlabel("vx [km/s]")
    ax.set_ylabel("vy [km/s]")
    ax.set_zlabel("vz [km/s]")
    ax.set_title("Comparison of velocity vectors in 3D")

    ax.legend(title="Function")
    plt.tight_layout()
    plt.show()

#if __name__ == "__main__":
#    plot_velocity_compare(velocity_functions=[earth_velocity_fast, earth_velocity])


ev_x, evf_x = [[], []]
ev_y, evf_y = [[], []]
ev_z, evf_z = [[], []]
mjd_values = []

EV, EVF = [[], []]

mjd_start = 60000
mjd_stop  = 60001
step = 0.05

while mjd_start < mjd_stop:
    evx = earth_velocity(mjd_start, xyz_list=[6378137, 0, 0])[0].value
    evfx = earth_velocity_fast(mjd_start)[0].value
    evy = earth_velocity(mjd_start, xyz_list=[6378137, 0, 0])[1].value
    evfy = earth_velocity_fast(mjd_start)[1].value
    evz = earth_velocity(mjd_start, xyz_list=[6378137, 0, 0])[2].value
    evfz = earth_velocity_fast(mjd_start)[2].value

    ev = np.sqrt( (earth_velocity(mjd_start, xyz_list=[6378137, 0, 0])[0].value)**2 
                 + (earth_velocity(mjd_start, xyz_list=[6378137, 0, 0])[1].value)**2 
                 + (earth_velocity(mjd_start, xyz_list=[6378137, 0, 0])[2].value)**2 )
    evf = np.sqrt( (earth_velocity_fast(mjd_start)[0].value)**2
                  + (earth_velocity_fast(mjd_start)[1].value)**2
                  + (earth_velocity_fast(mjd_start)[2].value)**2 )

    ev_x.append(evx)
    evf_x.append(evfx)
    ev_y.append(evy)
    evf_y.append(evfy)
    ev_z.append(evz)
    evf_z.append(evfz)

    mjd_values.append(mjd_start)
    EV.append(ev)
    EVF.append(evf)

    mjd_start += step

#Function compare earth_velocity vs earth_velocity_fast
print("========================================")
with open("figures/tabela_v/EV_EVF.txt", "w") as f:
    f.write("mjd\t\tev_x\tevf_x\tblad_wzg[%]\n")

    for mjd, ev, evf in zip(mjd_values, EV, EVF):
        f.write(f"{mjd:.1f}\t{ev:.3f}\t{evf:.3f}\t{100*abs((ev-evf)/ev):.3f}\n")


with open("figures/tabela_v/EV_xyz.txt", "w") as f:
    f.write("mjd\t\tev_x\tevf_x\t\tblad_wzg_x[%]\t\tev_y\tevf_y\t\tblad_wzg_y[%]\t\tev_z\tevf_z\t\tblad_wzg_z[%]\n")

    for mjd, evx, evfx, evy, evfy, evz, evfz in zip(mjd_values, ev_x, evf_x, ev_y, evf_y, ev_z, evf_z):
        f.write(f"{mjd:.1f}\t{evx:.3f}\t{evfx:.3f}\t{100*abs((evx-evfx)/evx):.3f}\t{evy:.3f}\t{evfy:.3f}\t{100*abs((evy-evfy)/evy):.3f}\t{evz:.3f}\t{evfz:.3f}\t{100*abs((evz-evfz)/evz):.3f}\n")