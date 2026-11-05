
"""
ВАЖНО:
Это самостоятельная расчётная модель первого приближения, а не точная
копия исторического полёта. Программа интегрирует уравнения движения,
гравитацию, вращение Земли, сопротивление атмосферы, расход топлива и
разделение боковых блоков. Угол тангажа задаётся встроенной программой управления. Двигатели
выключаются, когда апоцентр и перицентр одновременно оказываются рядом
с целевой высотой и эксцентриситет достаточно мал.
"""

import csv
import math
import os
try:
    import simplekml
except ImportError:
    print("Не найдена библиотека simplekml.")
    print("Установи её один раз командой: python -m pip install simplekml")
    input("Нажми Enter для закрытия...")
    raise SystemExit(1)

# =========================== НАСТРОЙКИ ==============================

DT = 0.01                          # шаг интегрирования, с
OUTPUT_INTERVAL = 0.01             # интервал записи CSV, с

MU = 3.98600435436096e14           # гравитационный параметр Земли, м^3/с^2
EARTH_RADIUS = 6_371_000.0         # средний радиус в этой модели, м
EARTH_OMEGA = 7.292115146706924e-5 # вращение Земли, рад/с
G0 = 9.80665

LAUNCH_LAT_DEG = 45.92             # 45°55′12″ с. ш.
LAUNCH_LON_DEG = 63.34777777777778 # 63°20′52″ в. д.
LAUNCH_AZIMUTH_DEG = 90.0          # 90° = на восток

TARGET_ALTITUDE_M = 300_000.0      # целевая высота орбиты
CD_AREA_M2 = 8.366                 # эффективное Cd*A, м^2
ATMOSPHERE_LIMIT_M = 140_000.0

# Четыре боковых блока (РД-107)
BOOSTER_COUNT = 4
BOOSTER_DRY_EACH_KG = 2_500.0
BOOSTER_PROPELLANT_EACH_KG = 45_090.0
BOOSTER_THRUST_SL_EACH_N = 83_000.0 * 9.80665
BOOSTER_THRUST_VAC_EACH_N = 102_000.0 * 9.80665
BOOSTER_ISP_SL_S = 243.0
BOOSTER_ISP_VAC_S = 299.0

# Центральный блок (РД-108)
CORE_DRY_KG = 5_371.47
CORE_PROPELLANT_KG = 96_694.70
CORE_THRUST_SL_N = 76_000.0 * 9.80665
CORE_THRUST_VAC_N = 96_000.0 * 9.80665
CORE_ISP_SL_S = 239.0
CORE_ISP_VAC_S = 302.0

# Таблица желаемого угла траектории к местному горизонту.
# Формат: (высота, градусы). Между точками линейная интерполяция.
# 90° — вертикально вверх; 0° — горизонтально.
"""
Эти углы подобраны численным расчётом так, чтобы получить
почти круговую орбиту около TARGET_ALTITUDE_M в данной модели.
Угол отсчитывается от местного горизонта: отрицательный угол
означает небольшой наклон тяги ниже горизонта для гашения
избыточной радиальной скорости перед выключением двигателя.
"""
PITCH_PROGRAM = [
    
    (0.0,        90.0),
    (10_000.0,   86.0),
    (25_000.0,   76.0),
    (50_000.0,   59.14),
    (80_000.0,   40.12),
    (110_000.0,  11.58),
    (150_000.0,  -6.87),
    (200_000.0, -13.61),
    (250_000.0, -10.30),
    (300_000.0, -10.30),
]

TARGET_ALTITUDE_TOLERANCE_M = 10_000.0  # допустимое отклонение апо-/перицентра
TARGET_ECCENTRICITY_MAX = 0.003         # критерий почти круговой орбиты

MAX_FLIGHT_TIME_S = 600.0

# ====================================================================

def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def mul(a, k):
    return (a[0] * k, a[1] * k, a[2] * k)


def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    n = norm(a)
    if n < 1e-12:
        return (0.0, 0.0, 0.0)
    return mul(a, 1.0/n)


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])


def atmosphere_density(h):
    """Упрощённая экспоненциальная атмосфера."""
    if h < 0.0:
        h = 0.0
    if h >= ATMOSPHERE_LIMIT_M:
        return 0.0
    return 1.225 * math.exp(-h / 8500.0)


def pressure_fraction(h):
    """Приближённое отношение давления к уровню моря."""
    return math.exp(-max(0.0, h) / 8500.0) if h < ATMOSPHERE_LIMIT_M else 0.0


def interp_pitch(h):
    if h <= PITCH_PROGRAM[0][0]:
        return PITCH_PROGRAM[0][1]
    for i in range(len(PITCH_PROGRAM)-1):
        h0, p0 = PITCH_PROGRAM[i]
        h1, p1 = PITCH_PROGRAM[i+1]
        if h <= h1:
            f = (h-h0)/(h1-h0)
            return p0 + f*(p1-p0)
    return PITCH_PROGRAM[-1][1]


def orbital_parameters(r, v):
    rn = norm(r)
    vn = norm(v)
    hvec = cross(r, v)
    hmag = norm(hvec)
    energy = 0.5*vn*vn - MU/rn
    evec = sub(mul(cross(v, hvec), 1.0/MU), mul(r, 1.0/rn))
    ecc = norm(evec)
    if energy < 0.0 and ecc < 1.0:
        a = -MU/(2.0*energy)
        rp = a*(1.0-ecc)
        ra = a*(1.0+ecc)
        return max(0.0, rp-EARTH_RADIUS), max(0.0, ra-EARTH_RADIUS), ecc
    return 0.0, float("inf"), ecc


def inertial_to_geographic(r, time_s):
    """Переводит инерциальные координаты в долготу/широту Земли для KML."""
    angle = EARTH_OMEGA * time_s
    c = math.cos(angle)
    s = math.sin(angle)
    # Поворот ECI -> ECEF на отрицательный угол вращения Земли.
    x = c * r[0] + s * r[1]
    y = -s * r[0] + c * r[1]
    z = r[2]
    radius = math.sqrt(x*x + y*y + z*z)
    longitude = math.degrees(math.atan2(y, x))
    latitude = math.degrees(math.atan2(z, math.hypot(x, y)))
    altitude = radius - EARTH_RADIUS
    return longitude, latitude, altitude


def write_kml(path, trajectory_points):
    """Создаёт KML-трек средствами simplekml."""
    if not trajectory_points:
        raise ValueError("Нет точек траектории для записи KML")

    kml = simplekml.Kml()
    kml.document.name = "Траектория полёта Р-7 — модель"

    # Стартовая метка в точных координатах космодрома Байконур.
    start = kml.newpoint(
        name="Старт — космодром Байконур",
        coords=[(LAUNCH_LON_DEG, LAUNCH_LAT_DEG, 0.0)],
    )
    start.altitudemode = simplekml.AltitudeMode.absolute

    # Точки имеют формат (долгота, широта, высота над уровнем моря).
    flight = kml.newlinestring(
        name="Траектория полёта Р-7",
        coords=trajectory_points,
    )
    flight.altitudemode = simplekml.AltitudeMode.absolute
    flight.extrude = 0
    flight.tessellate = 1
    flight.style.linestyle.color = "ff0095ff"
    flight.style.linestyle.width = 3

    end_lon, end_lat, end_alt = trajectory_points[-1]
    finish = kml.newpoint(
        name="Конец расчёта / выключение двигателей",
        coords=[(end_lon, end_lat, end_alt)],
    )
    finish.altitudemode = simplekml.AltitudeMode.absolute


    kml.save(path, format=False)


def acceleration(state, mass, thrust_n, thrust_dir, h):
    r = state[0:3]
    v = state[3:6]
    rn = norm(r)
    gravity = mul(r, -MU/(rn**3))
    # Атмосфера вращается вместе с Землёй вокруг оси Z.
    v_atm = cross((0.0, 0.0, EARTH_OMEGA), r)
    v_rel = sub(v, v_atm)
    speed_rel = norm(v_rel)
    rho = atmosphere_density(h)
    drag = (0.0, 0.0, 0.0)
    if speed_rel > 0.0 and rho > 0.0:
        drag = mul(v_rel, -0.5*rho*CD_AREA_M2*speed_rel/mass)
    thrust_acc = mul(thrust_dir, thrust_n/mass) if mass > 0.0 else (0,0,0)
    return add(add(gravity, drag), thrust_acc)


def derivative(state, mass, thrust_n, thrust_dir):
    r = state[0:3]
    v = state[3:6]
    h = norm(r)-EARTH_RADIUS
    a = acceleration(state, mass, thrust_n, thrust_dir, h)
    return (v[0], v[1], v[2], a[0], a[1], a[2])


def rk4_step(state, dt, mass, thrust_n, thrust_dir):
    k1 = derivative(state, mass, thrust_n, thrust_dir)
    s2 = [state[i] + 0.5*dt*k1[i] for i in range(6)]
    k2 = derivative(s2, mass, thrust_n, thrust_dir)
    s3 = [state[i] + 0.5*dt*k2[i] for i in range(6)]
    k3 = derivative(s3, mass, thrust_n, thrust_dir)
    s4 = [state[i] + dt*k3[i] for i in range(6)]
    k4 = derivative(s4, mass, thrust_n, thrust_dir)
    return [state[i] + dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0
            for i in range(6)]


def main():
    lat = math.radians(LAUNCH_LAT_DEG)
    lon = math.radians(LAUNCH_LON_DEG)
    az = math.radians(LAUNCH_AZIMUTH_DEG)

    # Стартовая точка в инерциальной системе координат.
    r0 = (EARTH_RADIUS*math.cos(lat)*math.cos(lon),
          EARTH_RADIUS*math.cos(lat)*math.sin(lon),
          EARTH_RADIUS*math.sin(lat))
    radial = unit(r0)
    east = (-math.sin(lon), math.cos(lon), 0.0)
    north = (-math.sin(lat)*math.cos(lon),
             -math.sin(lat)*math.sin(lon), math.cos(lat))
    horizontal = unit(add(mul(north, math.cos(az)), mul(east, math.sin(az))))
    # Начальная скорость от вращения Земли.
    v0 = cross((0.0, 0.0, EARTH_OMEGA), r0)
    state = list(r0 + v0)

    booster_prop = BOOSTER_COUNT*BOOSTER_PROPELLANT_EACH_KG
    booster_dry = BOOSTER_COUNT*BOOSTER_DRY_EACH_KG
    core_prop = CORE_PROPELLANT_KG
    core_dry = CORE_DRY_KG
    boosters_attached = True
    t = 0.0
    rows = []
    kml_points = []
    last_cutoff_reason = "достигнут предел времени"

    def current_mass():
        return core_dry + core_prop + (booster_dry + booster_prop if boosters_attached else 0.0)

    while t <= MAX_FLIGHT_TIME_S:
        r = tuple(state[0:3])
        v = tuple(state[3:6])
        h = norm(r)-EARTH_RADIUS
        mass = current_mass()
        rho = atmosphere_density(h)
        pfrac = pressure_fraction(h)

        # Направление тяги: в плоскости радиус-вектора и горизонтальной скорости.
        radial = unit(r)
        v_horizontal = sub(v, mul(radial, dot(v, radial)))
        if norm(v_horizontal) < 1e-6:
            tangent = horizontal
        else:
            tangent = unit(v_horizontal)
        gamma = math.radians(interp_pitch(h))
        thrust_dir = unit(add(mul(radial, math.sin(gamma)),
                              mul(tangent, math.cos(gamma))))

        thrust = 0.0
        mdot = 0.0
        if boosters_attached and booster_prop > 0.0:
            tb = BOOSTER_COUNT * (
                BOOSTER_THRUST_SL_EACH_N*pfrac +
                BOOSTER_THRUST_VAC_EACH_N*(1.0-pfrac))
            isp_b = BOOSTER_ISP_SL_S*pfrac + BOOSTER_ISP_VAC_S*(1.0-pfrac)
            thrust += tb
            mdot += tb/(isp_b*G0)
        if core_prop > 0.0:
            tc = CORE_THRUST_SL_N*pfrac + CORE_THRUST_VAC_N*(1.0-pfrac)
            isp_c = CORE_ISP_SL_S*pfrac + CORE_ISP_VAC_S*(1.0-pfrac)
            thrust += tc
            mdot += tc/(isp_c*G0)

        # Орбитальные параметры до текущего шага.
        peri, apo, ecc = orbital_parameters(r, v)
        v_inertial = norm(v)
        v_surface = norm(sub(v, cross((0.0, 0.0, EARTH_OMEGA), r)))

        if int(round(t/OUTPUT_INTERVAL)) == len(rows):
            rows.append({
                "время (с)": t,
                "высота (м)": h,
                "поверхностная скорость (м/с)": v_surface,
                "орбитальная скорость (м/с)": v_inertial,
                "масса (кг)": mass,
                "апоцентр (м)": apo,
                "перицентр (м)": peri,
            })
            kml_points.append(inertial_to_geographic(r, t))

        # Главный критерий выключения: не только апоцентр, но и перицентр
        # должны находиться рядом с заданной высотой. Это не позволяет
        # выключить двигатель на вытянутой эллиптической траектории.
        radial_speed = dot(r, v)/norm(r)
        tangential_speed = math.sqrt(max(0.0, v_inertial*v_inertial-radial_speed*radial_speed))
        orbit_reached = (
            h > 150_000.0
            and math.isfinite(apo)
            and abs(apo - TARGET_ALTITUDE_M) <= TARGET_ALTITUDE_TOLERANCE_M
            and abs(peri - TARGET_ALTITUDE_M) <= TARGET_ALTITUDE_TOLERANCE_M
            and ecc <= TARGET_ECCENTRICITY_MAX
        )
        if orbit_reached:
            last_cutoff_reason = "достигнута почти круговая орбита на целевой высоте"
            break

        if thrust <= 0.0 or mass <= 0.0:
            last_cutoff_reason = "топливо закончилось"
            break

        # Половинный шаг расхода массы для согласования массы с интегрированием.
        dm = min(mdot*DT, booster_prop+core_prop)
        # Сначала расходуются боковые блоки и центральный двигатель одновременно.
        if boosters_attached and booster_prop > 0.0:
            tb_now = BOOSTER_COUNT * (
                BOOSTER_THRUST_SL_EACH_N*pfrac + BOOSTER_THRUST_VAC_EACH_N*(1.0-pfrac))
            isp_b_now = BOOSTER_ISP_SL_S*pfrac + BOOSTER_ISP_VAC_S*(1.0-pfrac)
            mdot_b = tb_now/(isp_b_now*G0)
        else:
            mdot_b = 0.0
        if core_prop > 0.0:
            tc_now = CORE_THRUST_SL_N*pfrac + CORE_THRUST_VAC_N*(1.0-pfrac)
            isp_c_now = CORE_ISP_SL_S*pfrac + CORE_ISP_VAC_S*(1.0-pfrac)
            mdot_c = tc_now/(isp_c_now*G0)
        else:
            mdot_c = 0.0

        old_booster_prop = booster_prop
        old_core_prop = core_prop
        if boosters_attached and booster_prop > 0.0:
            booster_prop = max(0.0, booster_prop-mdot_b*DT)
        if core_prop > 0.0:
            core_prop = max(0.0, core_prop-mdot_c*DT)

        # При исчерпании топлива боковых блоков сбрасываем их сухую массу.
        if boosters_attached and booster_prop <= 1e-6:
            boosters_attached = False
            booster_dry = 0.0

        # RK4 с массой на середине шага.
        mass_mid = max(core_dry + core_prop +
                       (booster_dry + booster_prop if boosters_attached else 0.0), 1.0)
        state = rk4_step(state, DT, mass_mid, thrust, thrust_dir)
        t += DT

        if core_prop <= 1e-6:
            last_cutoff_reason = "израсходовано топливо центрального блока"
            break

    # Пересчитываем финальные параметры.
    r = tuple(state[0:3])
    v = tuple(state[3:6])
    final_h = norm(r)-EARTH_RADIUS
    final_peri, final_apo, final_ecc = orbital_parameters(r, v)
    final_mass = current_mass()

    output_dir = os.path.dirname(os.path.abspath(__file__))
    output_name = os.path.join(output_dir, "vostok_simulation.csv")
    kml_name = os.path.join(output_dir, "vostok_trajectory.kml")
    csv_columns = [
        "время (с)",
        "высота (м)",
        "поверхностная скорость (м/с)",
        "орбитальная скорость (м/с)",
        "масса (кг)",
        "апоцентр (м)",
        "перицентр (м)",
    ]
    with open(output_name, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=csv_columns, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)

    write_kml(kml_name, kml_points)

    print("\nРасчёт завершён.")
    print(f"Время полёта: {t:.2f} с")
    print(f"Высота в конце расчёта: {final_h/1000:.2f} км")
    print(f"Инерциальная скорость: {norm(v):.2f} м/с")
    print(f"Масса в конце: {final_mass:.2f} кг")
    print(f"Апоцентр: {final_apo/1000:.2f} км")
    print(f"Перицентр: {final_peri/1000:.2f} км")
    print(f"Эксцентриситет: {final_ecc:.6f}")
    print(f"Причина остановки: {last_cutoff_reason}")
    reached = (
        math.isfinite(final_apo)
        and abs(final_apo - TARGET_ALTITUDE_M) <= TARGET_ALTITUDE_TOLERANCE_M
        and abs(final_peri - TARGET_ALTITUDE_M) <= TARGET_ALTITUDE_TOLERANCE_M
        and final_ecc <= TARGET_ECCENTRICITY_MAX
    )
    print(f"Целевая почти круговая орбита достигнута: {'ДА' if reached else 'НЕТ'}")
    if not reached:
        print("ВНИМАНИЕ: расчёт не достиг заданной орбиты до остановки двигателей/лимита времени.")
        if math.isfinite(final_apo):
            print(f"Отклонение апоцентра: {(final_apo-TARGET_ALTITUDE_M)/1000:+.2f} км")
        print(f"Отклонение перицентра: {(final_peri-TARGET_ALTITUDE_M)/1000:+.2f} км")
    else:
        print(f"Остаток топлива: {booster_prop + core_prop:.2f} кг")
    print(f"CSV сохранён: {os.path.abspath(output_name)}")
    print(f"KML сохранён: {os.path.abspath(kml_name)}")
    input("\nНажми Enter для закрытия...")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"\nОшибка: {exc}")
        input("Нажми Enter для закрытия...")
