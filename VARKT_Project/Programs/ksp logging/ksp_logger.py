
"""
Оптимизированный низконагрузочный логгер KSP/RSS через kRPC.
Сохраняет в CSV: время, высоту, поверхностную/орбитальную скорости, массу, апоцентр и перицентр.
Гео-координаты для KML рассчитываются в реальном времени во время полета.
Файлы пишутся строго в директорию скрипта.
"""

import csv
import json
import math
import os
import sys
import time
from datetime import datetime

try:
    import krpc
except ImportError:
    print("Не найден пакет krpc. Установи: py -m pip install krpc")
    raise SystemExit(1)

try:
    import simplekml
except ImportError:
    print("Внимание: Пакет 'simplekml' не установлен. Автоэкспорт в KML будет пропущен.")
    print("Чтобы включить KML, выполни: pip install simplekml")
    simplekml = None

# ---------------------------- НАСТРОЙКИ ----------------------------
RPC_PORT = 50000
STREAM_PORT = 50001
SAMPLE_INTERVAL_S = 0.05  # 20 Гц
FLUSH_EVERY_N_ROWS = 20
OUTPUT_PREFIX = "ksp_logger"
# ------------------------------------------------------------------

def finite_or_blank(value):
    """Защита от NaN/inf."""
    try:
        number = float(value)
        return number if math.isfinite(number) else ""
    except (TypeError, ValueError):
        return ""

def main():
    conn = None
    streams = []
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    csv_name = os.path.join(script_dir, f"{OUTPUT_PREFIX}_{stamp}.csv")
    meta_name = os.path.join(script_dir, f"{OUTPUT_PREFIX}_{stamp}_body.json")
    kml_name = os.path.join(script_dir, f"{OUTPUT_PREFIX}_{stamp}_trajectory.kml")
    
    # Список для накопления готовых KML-координат во время полета
    kml_coords = []

    try:
        conn = krpc.connect(
            name="Vostok Short Trajectory Logger",
            rpc_port=RPC_PORT,
            stream_port=STREAM_PORT,
        )
        vessel = conn.space_center.active_vessel
        if vessel is None:
            raise RuntimeError("В KSP не выбран активный корабль.")

        body = vessel.orbit.body
        body_fixed_frame = body.reference_frame
        orbit_frame = body.non_rotating_reference_frame
        surface_frame = vessel.surface_reference_frame
        
        space_center = conn.space_center
        flight = vessel.flight(body_fixed_frame)
        orbit = vessel.orbit

        metadata = {
            "vessel_name": vessel.name,
            "body_name": body.name,
            "target_sample_interval_s": SAMPLE_INTERVAL_S,
            "target_sample_rate_hz": 1.0 / SAMPLE_INTERVAL_S,
            "logger_version": "short-csv-1.1",
        }

        with open(meta_name, "w", encoding="utf-8") as meta_file:
            json.dump(metadata, meta_file, ensure_ascii=False, indent=2)

        def add_stream(method, *args):
            stream = conn.add_stream(method, *args)
            streams.append(stream)
            return stream

        # Создание потоков данных
        ut_stream = add_stream(getattr, space_center, "ut")
        alt_stream = add_stream(getattr, flight, "mean_altitude")
        v_surf_stream = add_stream(vessel.velocity, surface_frame)
        v_orbit_stream = add_stream(vessel.velocity, orbit_frame)
        mass_stream = add_stream(getattr, vessel, "mass")
        apo_stream = add_stream(getattr, orbit, "apoapsis_altitude")
        peri_stream = add_stream(getattr, orbit, "periapsis_altitude")
        
        # Потоки для KML-координат в реальном времени
        lat_stream = add_stream(getattr, flight, "latitude")
        lon_stream = add_stream(getattr, flight, "longitude")

        start_ut = space_center.ut

        headers = [
            "Время (с)",
            "Высота (м)",
            "поверхностная скорость (м/с)",
            "орбитальная скорость (м/с)",
            "масса (кг)",
            "апоцентр (м)",
            "перицентр (м)"
        ]

        print(f"Подключено: {vessel.name}; тело: {body.name}")
        print(f"CSV (сокращенный): {csv_name}")
        print(f"Целевая частота: {1.0 / SAMPLE_INTERVAL_S:.0f} Гц")
        print("Запусти ракету.\nОстановка: Ctrl+C.")

        row_count = 0
        next_tick = time.monotonic()

        with open(csv_name, "w", newline="", encoding="utf-8-sig") as csv_file:
            writer = csv.writer(csv_file, delimiter=";")
            writer.writerow(headers)

            try:
                while True:
                    current_ut = ut_stream()
                    t = finite_or_blank(current_ut - start_ut)
                    
                    v_surf = v_surf_stream()
                    v_orbit = v_orbit_stream()
                    
                    speed_surf = math.sqrt(v_surf[0]**2 + v_surf[1]**2 + v_surf[2]**2)
                    speed_orbit = math.sqrt(v_orbit[0]**2 + v_orbit[1]**2 + v_orbit[2]**2)
                    
                    current_alt = alt_stream()

                    row = [
                        t,
                        finite_or_blank(current_alt),
                        finite_or_blank(speed_surf),
                        finite_or_blank(speed_orbit),
                        finite_or_blank(mass_stream()),
                        finite_or_blank(apo_stream()),
                        finite_or_blank(peri_stream()),
                    ]
                    
                    writer.writerow(row)
                    row_count += 1

                    # Собираем KML каждую 5-ю итерацию (~4 Гц), чтобы разгрузить лог
                    if row_count % 5 == 0:
                        try:
                            kml_alt = current_alt if current_alt >= 0 else 0
                            kml_coords.append((lon_stream(), lat_stream(), kml_alt))
                        except Exception:
                            pass

                    if row_count % FLUSH_EVERY_N_ROWS == 0:
                        csv_file.flush()

                    next_tick += SAMPLE_INTERVAL_S
                    delay = next_tick - time.monotonic()
                    if delay > 0:
                        time.sleep(delay)
                    else:
                        next_tick = time.monotonic()

            except KeyboardInterrupt:
                print("\nЗапись остановлена пользователем.")
            finally:
                csv_file.flush()
                print(f"Записано строк в CSV: {row_count}")

        # === ВЫДЕЛЕННЫЙ БЛОК: ЗАПИСЬ ГОТОВОГО ТРЕКА KML ===
        if simplekml and kml_coords:
            print("Сохраняю собранную траекторию в KML...")
            try:
                kml = simplekml.Kml()
                ls = kml.newlinestring(name=f"Траектория {vessel.name}")
                ls.coords = kml_coords
                ls.altitudemode = simplekml.AltitudeMode.absolute
                ls.extrude = 1
                ls.style.linestyle.color = simplekml.Color.red
                ls.style.linestyle.width = 3
                ls.style.polystyle.color = simplekml.Color.changealphaint(40, simplekml.Color.red)
                
                kml.save(kml_name)
                print(f"KML-траектория успешно сохранена: {kml_name}")
            except Exception as kml_err:
                print(f"Не удалось записать KML файл: {kml_err}")
        else:
            print("Не удалось рассчитать координаты для KML или пакет simplekml отсутствует.")

    except Exception as exc:
        print(f"Ошибка логгера: {exc}", file=sys.stderr)
        raise
    finally:
        for stream in streams:
            try:
                stream.remove()
            except Exception:
                pass
        if conn is not None:
            try:
                conn.close()
            except Exception:
                pass

if __name__ == "__main__":
    main()
