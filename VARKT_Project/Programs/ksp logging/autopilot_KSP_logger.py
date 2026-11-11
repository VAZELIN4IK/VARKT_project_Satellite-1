"""
Автопилот kRPC для KSP (RP-1).
Выход на орбиту за одно включение (учет 1 зажигания).
Отсечка по достижении апоцентра и малого эксцентриситета.
Исправлено разделение орбитальной и поверхностной скорости.
Частота записи данных в CSV увеличена до 20 Гц (20 раз в секунду).
"""

import csv
import math
import os
import sys
import time

try:
    import krpc
except ImportError:
    print("Не найден пакет krpc. Установить: py -m pip install krpc")
    raise SystemExit(1)

try:
    import simplekml
except ImportError:
    print("Внимание: пакет 'simplekml' не установлен. KML-файл не будет создан.")
    simplekml = None

# =========================== НАСТРОЙКИ ==============================

RPC_PORT = 50000
STREAM_PORT = 50001

TARGET_ALTITUDE_M = 300_000.0

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

LAUNCH_AZIMUTH_DEG = 90.0
OUTPUT_INTERVAL = 0.05
MAX_FLIGHT_TIME_S = 600.0

# ====================================================================

def interp_pitch(altitude_m: float) -> float:
    if altitude_m <= PITCH_PROGRAM[0][0]:
        return PITCH_PROGRAM[0][1]
    for i in range(len(PITCH_PROGRAM) - 1):
        h0, p0 = PITCH_PROGRAM[i]
        h1, p1 = PITCH_PROGRAM[i + 1]
        if altitude_m <= h1:
            f = (altitude_m - h0) / (h1 - h0)
            return p0 + f * (p1 - p0)
    return PITCH_PROGRAM[-1][1]

def write_kml(path: str, points: list[tuple[float, float, float]]):
    if not points or simplekml is None:
        return
    kml = simplekml.Kml()
    kml.document.name = "Траектория полёта Р-7 (kRPC / KSP)"
    start_lon, start_lat, start_alt = points[0]
    kml.newpoint(name="Старт", coords=[(start_lon, start_lat, start_alt)], altitudemode=simplekml.AltitudeMode.absolute)
    
    line = kml.newlinestring(name="Траектория", coords=points)
    line.altitudemode = simplekml.AltitudeMode.absolute
    line.extrude = 0
    line.tessellate = 1
    line.style.linestyle.color = "ff0095ff"
    line.style.linestyle.width = 3

    end_lon, end_lat, end_alt = points[-1]
    kml.newpoint(name="Выключение двигателей", coords=[(end_lon, end_lat, end_alt)], altitudemode=simplekml.AltitudeMode.absolute)
    kml.save(path, format=False)

def save_data(rows, kml_points, csv_path, kml_path, final_met=None, final_alt=None, final_mass=None, final_apo=None, final_peri=None, final_ecc=None, cutoff_reason=""):
    if not rows:
        return
    csv_columns = ["время (с)", "высота (м)", "поверхностная скорость (м/с)", "орбитальная скорость (м/с)", "масса (кг)", "апоцентр (м)", "перицентр (м)"]
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=csv_columns, delimiter=";")
        writer.writeheader()
        writer.writerows(rows)
    print(f"CSV сохранён: {csv_path} (Всего строк: {len(rows)})")
    
    if simplekml is not None and kml_points:
        try:
            write_kml(kml_path, kml_points)
            print(f"KML сохранён: {kml_path}")
        except Exception as e:
            print(f"Ошибка записи KML: {e}")
        
    print("\n--- Отчёт о полёте ---")
    if final_met is not None:
        print(f"Время полёта: {final_met:.2f} с")
        print(f"Высота в конце: {final_alt / 1000:.2f} км")
        print(f"Масса в конце: {final_mass:.2f} кг")
        if final_apo is not None and not math.isnan(final_apo):
            print(f"Апоцентр: {final_apo / 1000:.2f} км")
            print(f"Перицентр: {final_peri / 1000:.2f} км")
            print(f"Эксцентриситет: {final_ecc:.6f}")
    print(f"Причина остановки: {cutoff_reason}")
    print("----------------------")

def main():
    conn = None
    vessel = None
    ap = None
    start_time = None
    rows = []
    kml_points = []
    cutoff_reason = "Неизвестно"
    data_saved = False
    
    output_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(output_dir, "flight_data.csv")
    kml_path = os.path.join(output_dir, "trajectory.kml")

    try:
        conn = krpc.connect(name="R-7 Autopilot", rpc_port=RPC_PORT, stream_port=STREAM_PORT)
        space_center = conn.space_center
        vessel = space_center.active_vessel
        if vessel is None:
            raise RuntimeError("В KSP не выбран активный корабль.")

        vessel.control.sas = False
        vessel.control.rcs = False

        ap = vessel.auto_pilot
        ap.engaged = True
        ap.target_pitch_and_heading(90.0, LAUNCH_AZIMUTH_DEG)

        vessel.control.throttle = 1.0

        print("=====================================================")
        print("Автопилот готов.")
        print("1. Включите двигатели вручную (Ступень 4 / Пробел).")
        print("2. Автопилот сам отпустит площадку.")
        print("3. Отсечка произойдет при достижении орбиты.")
        print(f"4. Запись данных: 20 раз в секунду.")
        print("=====================================================")
        print("Для АВАРИЙНОЙ ОСТАНОВКИ нажмите Ctrl+C.")

        orbital_frame = vessel.orbit.body.non_rotating_reference_frame
        surface_frame = vessel.orbit.body.reference_frame

        flight_orbital = vessel.flight(orbital_frame)
        flight_surface = vessel.flight(surface_frame)

        # --- ОЖИДАНИЕ ЗАЖИГАНИЯ ДВИГАТЕЛЕЙ ---
        while vessel.thrust < 100.0:
            time.sleep(0.1)

        print("Двигатели запущены! Ожидание выхода на тягу...")

        # --- ОЖИДАНИЕ TWR >= 1 ---
        while vessel.thrust < (vessel.mass * vessel.orbit.body.surface_gravity):
            time.sleep(0.1)

        print("Тяга достаточна. Отсоединяем стартовые зажимы.")
        vessel.control.activate_next_stage()

        # Время t=0 начинается СТРОГО после отрыва от площадки
        start_time = vessel.met
        last_output_time = -OUTPUT_INTERVAL
        last_orbit_print = 0
        cutoff_reason = "достигнут предел времени (600 с)"

        # --- ОСНОВНОЙ ЦИКЛ ПОЛЁТА ---
        while vessel.met - start_time <= MAX_FLIGHT_TIME_S:
            met = vessel.met
            altitude = flight_surface.mean_altitude

            # Управление тангажом
            pitch_deg = interp_pitch(altitude)
            ap.target_pitch_and_heading(pitch_deg, LAUNCH_AZIMUTH_DEG)

            # Сбор и запись телеметрии (20 раз в секунду)
            if met - last_output_time >= OUTPUT_INTERVAL:
                last_output_time = met
                
                orbital_speed = flight_orbital.speed
                surface_speed = flight_surface.speed
                
                apo = vessel.orbit.apoapsis_altitude
                peri = vessel.orbit.periapsis_altitude
                ecc = vessel.orbit.eccentricity
                mass = vessel.mass

                rows.append({
                    "время (с)": met - start_time,
                    "высота (м)": altitude,
                    "поверхностная скорость (м/с)": surface_speed,
                    "орбитальная скорость (м/с)": orbital_speed,
                    "масса (кг)": mass,
                    "апоцентр (м)": apo,
                    "перицентр (м)": peri,
                })

                lat = flight_surface.latitude
                lon = flight_surface.longitude
                kml_points.append((lon, lat, altitude))

            # --- ПРОВЕРКА УСЛОВИЯ ВЫКЛЮЧЕНИЯ ДВИГАТЕЛЕЙ ---
            try:
                apo = vessel.orbit.apoapsis_altitude
                peri = vessel.orbit.periapsis_altitude
                ecc = vessel.orbit.eccentricity

                if altitude > 150_000.0:
                    if (TARGET_ALTITUDE_M * 0.925 <= apo <= TARGET_ALTITUDE_M * 1.075)\
                        and (TARGET_ALTITUDE_M * 0.925 <= peri <= TARGET_ALTITUDE_M * 1.075) and ecc <= 0.01:
                        cutoff_reason = f"Достигнута орбита (Апо: {apo/1000:.1f} км, Экс: {ecc:.3f})"
                        print(f"\n[!] ОТСЕЧКА: Орбита достигнута!")
                        break
                        
                    if apo > TARGET_ALTITUDE_M + 100_000:
                        cutoff_reason = f"Перелет! Апоцентр {apo/1000:.1f} км (страховка)"
                        print(f"\n[!] ОТСЕЧКА: Апоцентр улетел слишком высоко!")
                        break
                    
            except Exception:
                pass

            time.sleep(0.015)

        # Штатное завершение полёта
        print(f"\n[!] Выключение двигателя. Причина: {cutoff_reason}")
        
        vessel.control.throttle = 0.0
        ap.engaged = False

        final_met = vessel.met - start_time
        final_alt = flight_surface.mean_altitude
        final_mass = vessel.mass
        try:
            final_apo = vessel.orbit.apoapsis_altitude
            final_peri = vessel.orbit.periapsis_altitude
            final_ecc = vessel.orbit.eccentricity
        except Exception:
            final_apo = final_peri = final_ecc = float("nan")

        save_data(rows, kml_points, csv_path, kml_path, final_met, final_alt, final_mass, final_apo, final_peri, final_ecc, cutoff_reason)
        data_saved = True

    except KeyboardInterrupt:
        cutoff_reason = "Аварийная остановка пользователем (Ctrl+C)"
        print("\n\n[!] Получен сигнал Ctrl+C. Аварийное завершение...")
    except Exception as exc:
        cutoff_reason = f"Ошибка выполнения: {exc}"
        print(f"\n[!] Ошибка автопилота: {exc}", file=sys.stderr)
    
    finally:
        if vessel is not None:
            try:
                vessel.control.throttle = 0.0
            except Exception:
                pass
        if ap is not None:
            try:
                ap.engaged = False
            except Exception:
                pass
                
        if not data_saved and rows:
            print("Экстренное сохранение собранных данных...")
            final_met = final_alt = final_mass = final_apo = final_peri = final_ecc = None
            if vessel is not None and start_time is not None:
                try:
                    final_met = vessel.met - start_time
                    final_alt = vessel.flight().mean_altitude
                    final_mass = vessel.mass
                    final_apo = vessel.orbit.apoapsis_altitude
                    final_peri = vessel.orbit.periapsis_altitude
                    final_ecc = vessel.orbit.eccentricity
                except Exception:
                    pass
            save_data(rows, kml_points, csv_path, kml_path, final_met, final_alt, final_mass, final_apo, final_peri, final_ecc, cutoff_reason)
            
        if conn is not None:
            try:
                conn.close()
            except Exception:
                pass

if __name__ == "__main__":
    main()
