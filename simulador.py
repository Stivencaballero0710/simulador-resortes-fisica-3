from __future__ import annotations

import json
import math
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit

BASE = Path(__file__).resolve().parent
VERSION = "2.0.0"


def Rigidez(G_Pa, d_m, D_m, n_espiras):
    return G_Pa * (d_m ** 4) / (8.0 * (D_m ** 3) * n_espiras)


def Cargar_Configuracion(Ruta_Archivo=BASE / "parametros.json", Semilla_Aleatoria=None):
    if Path(Ruta_Archivo).exists():
        Configuracion = json.loads(Path(Ruta_Archivo).read_text(encoding="utf-8"))
    else:
        Configuracion = {
            "fisica": {
                "g": 9.81,
                "d": 0.001,
                "D": 0.025,
                "n": 20,
                "m_ref": 0.25,
                "m_val": 0.3,
                "m_accesorio": 0.08,
                "A0": 0.01,
                "limite": 0.22,
                "duracion": 15,
                "dt": 0.005,
                "sigma": 0.00005,
                "paso_cuantizacion": 0.0001,
                "zeta_ref": 0.02,
                "semilla": 2250191,
                "rtol": 1e-10,
                "atol": 1e-12,
                "max_step": 0.03
            },
            "masas_estaticas": [0.15, 0.2, 0.25, 0.35],
            "repeticiones": 3,
            "barrido_zeta": [0, 0.02, 0.1, 0.3, 1, 1.5],
            "sensibilidad": {"G": 0.05, "d": 0.01, "D": 0.01},
            "materiales": [
                {"codigo": "AC", "nombre": "Acero Oteva 90", "G_Pa": 79600000000},
                {"codigo": "IN", "nombre": "Inoxidable AISI 302", "E_Pa": 185000000000, "nu": 0.3},
                {"codigo": "TI", "nombre": "Titanio B120VCA", "G_Pa": 41400000000}
            ]
        }

    Parametros_Fisicos = Configuracion["fisica"]
    if Semilla_Aleatoria is not None:
        Parametros_Fisicos["semilla"] = int(Semilla_Aleatoria)

    for Mat in Configuracion["materiales"]:
        if "E_Pa" in Mat:
            Mat["G_Pa"] = Mat["E_Pa"] / (2.0 * (1.0 + Mat["nu"]))
        Mat["k_N_m"] = Rigidez(
            Mat["G_Pa"],
            Parametros_Fisicos["d"],
            Parametros_Fisicos["D"],
            Parametros_Fisicos["n"]
        )
        Mat["b_ref_kg_s"] = 2.0 * Parametros_Fisicos["zeta_ref"] * math.sqrt(Mat["k_N_m"] * Parametros_Fisicos["m_ref"])

    return Configuracion


def Validar_Sistema(Masa, Rigidez_K, Amortiguamiento_B, Amplitud_A):
    if not all(np.isfinite(V) for V in (Masa, Rigidez_K, Amortiguamiento_B, Amplitud_A)):
        raise ValueError("Parametros deben ser finitos.")
    if Masa <= 0 or Rigidez_K <= 0 or Amortiguamiento_B < 0 or Amplitud_A < 0:
        raise ValueError("Masa y rigidez positivas; amortiguamiento y amplitud no negativos.")


def Integrar(Masa, Rigidez_K, Amortiguamiento_B, Amplitud_A=0.010, Vector_Tiempo=None,
             Tolerancia_Rel=1e-10, Tolerancia_Abs=1e-12, Paso_Maximo=0.03):
    Validar_Sistema(Masa, Rigidez_K, Amortiguamiento_B, Amplitud_A)
    if Vector_Tiempo is None:
        Vector_Tiempo = np.arange(3001, dtype=float) * 0.005
    Vector_Tiempo = np.asarray(Vector_Tiempo, dtype=float)

    if (Vector_Tiempo.ndim != 1 or len(Vector_Tiempo) < 2 or not np.all(np.isfinite(Vector_Tiempo))
            or Vector_Tiempo[0] != 0 or np.any(np.diff(Vector_Tiempo) <= 0)):
        raise ValueError("Vector de tiempo invalido.")

    Solucion = solve_ivp(
        lambda _, Y: [Y[1], -(Amortiguamiento_B * Y[1] + Rigidez_K * Y[0]) / Masa],
        (0, float(Vector_Tiempo[-1])), [Amplitud_A, 0.0],
        t_eval=Vector_Tiempo, method="DOP853",
        rtol=Tolerancia_Rel, atol=Tolerancia_Abs, max_step=Paso_Maximo
    )
    if not Solucion.success:
        raise RuntimeError(Solucion.message)

    Posicion, Velocidad = Solucion.y
    Aceleracion = -(Amortiguamiento_B * Velocidad + Rigidez_K * Posicion) / Masa
    Energia = 0.5 * Masa * (Velocidad ** 2) + 0.5 * Rigidez_K * (Posicion ** 2)
    return Vector_Tiempo, Posicion, Velocidad, Aceleracion, Energia


def Analitica(Vector_Tiempo, Masa, Rigidez_K, Amortiguamiento_B, Amplitud_A=0.010):
    Validar_Sistema(Masa, Rigidez_K, Amortiguamiento_B, Amplitud_A)
    Vector_Tiempo = np.asarray(Vector_Tiempo, dtype=float)
    W0 = math.sqrt(Rigidez_K / Masa)
    Gamma = Amortiguamiento_B / (2.0 * Masa)
    Zeta = Gamma / W0

    if abs(Zeta - 1.0) <= 1e-12:
        Exp_Term = np.exp(-W0 * Vector_Tiempo)
        Posicion = Amplitud_A * (1.0 + W0 * Vector_Tiempo) * Exp_Term
        Velocidad = -Amplitud_A * (W0 ** 2) * Vector_Tiempo * Exp_Term
    elif Zeta < 1.0:
        Wd = math.sqrt(W0 * W0 - Gamma * Gamma)
        Exp_Term = np.exp(-Gamma * Vector_Tiempo)
        Posicion = Amplitud_A * Exp_Term * (np.cos(Wd * Vector_Tiempo) + (Gamma / Wd) * np.sin(Wd * Vector_Tiempo))
        Velocidad = -Amplitud_A * (W0 ** 2 / Wd) * Exp_Term * np.sin(Wd * Vector_Tiempo)
    else:
        Q = math.sqrt(Zeta * Zeta - 1.0)
        R1 = -W0 / (Zeta + Q)
        R2 = -W0 * (Zeta + Q)
        C1 = -Amplitud_A * R2 / (R1 - R2)
        C2 = Amplitud_A * R1 / (R1 - R2)
        Posicion = C1 * np.exp(R1 * Vector_Tiempo) + C2 * np.exp(R2 * Vector_Tiempo)
        Velocidad = C1 * R1 * np.exp(R1 * Vector_Tiempo) + C2 * R2 * np.exp(R2 * Vector_Tiempo)
    return Posicion, Velocidad


def Periodo(Masa, Rigidez_K, Amortiguamiento_B):
    Wd_Cuadrado = Rigidez_K / Masa - (Amortiguamiento_B / (2.0 * Masa)) ** 2
    return 2.0 * math.pi / math.sqrt(Wd_Cuadrado) if Wd_Cuadrado > 0 else np.nan


def Regimen(Zeta):
    if Zeta == 0:
        return "Sin amortiguar"
    if abs(Zeta - 1.0) <= 1e-12:
        return "Critico"
    return "Subamortiguado" if Zeta < 1.0 else "Sobreamortiguado"


def Tiempo_Al_5(Masa, Rigidez_K, Amortiguamiento_B):
    W0 = math.sqrt(Rigidez_K / Masa)
    Zeta = Amortiguamiento_B / (2.0 * math.sqrt(Rigidez_K * Masa))
    if Zeta == 0:
        return math.inf
    if Zeta < 1.0 - 1e-12:
        return math.log(1.0 / (0.05 * math.sqrt(1.0 - Zeta * Zeta))) / (Zeta * W0)

    def Cruce(T_Val):
        return float(Analitica(np.array([T_Val]), Masa, Rigidez_K, Amortiguamiento_B, 1.0)[0][0]) - 0.05

    Limite_Sup = 10.0 / W0
    while Cruce(Limite_Sup) > 0:
        Limite_Sup *= 2.0
        if Limite_Sup > 1e12:
            raise RuntimeError("No se pudo acotar el tiempo al 5 por ciento.")
    return float(brentq(Cruce, 0, Limite_Sup, xtol=1e-12))


def Senal_Sintetica(Limpia, Generador_Rng, Parametros):
    Paso_Q = Parametros["paso_cuantizacion"]
    Ruido = Generador_Rng.normal(0, Parametros["sigma"], size=np.shape(Limpia))
    return np.round((np.asarray(Limpia) + Ruido) / Paso_Q) * Paso_Q


def Modelo_Ajuste(T_Val, Amplitud_A, Gamma_Val, Wd_Val, Offset_Val):
    return Amplitud_A * np.exp(-Gamma_Val * T_Val) * (np.cos(Wd_Val * T_Val) + Gamma_Val / Wd_Val * np.sin(Wd_Val * T_Val)) + Offset_Val


def Ajustar_Dinamica(Vector_Tiempo, Vector_Posicion, Masa):
    Dt = float(Vector_Tiempo[1] - Vector_Tiempo[0])
    Espectro = np.abs(np.fft.rfft(Vector_Posicion - np.mean(Vector_Posicion)))
    Frecuencias = np.fft.rfftfreq(len(Vector_Tiempo), Dt)
    Wd0 = 2.0 * math.pi * Frecuencias[1 + int(np.argmax(Espectro[1:]))]
    A0 = max(float(np.max(np.abs(Vector_Posicion))), 1e-8)
    Parametros_Opt, Covarianza = curve_fit(
        Modelo_Ajuste, Vector_Tiempo, Vector_Posicion,
        p0=[A0, 0.02 * Wd0, Wd0, 0],
        bounds=([0, 0, 0.5 * Wd0, -2 * A0], [5 * A0, 2 * Wd0, 1.5 * Wd0, 2 * A0]),
        method="trf", x_scale="jac", max_nfev=20000,
        ftol=1e-12, xtol=1e-12, gtol=1e-12
    )
    A_Est, Gamma_Est, Wd_Est, Offset_Est = map(float, Parametros_Opt)
    Errores = np.sqrt(np.maximum(np.diag(Covarianza), 0))
    Residuos = Vector_Posicion - Modelo_Ajuste(Vector_Tiempo, *Parametros_Opt)
    return dict(
        T_s=2.0 * math.pi / Wd_Est,
        b_kg_s=2.0 * Masa * Gamma_Est,
        A_m=A_Est,
        gamma_s_1=Gamma_Est,
        wd_rad_s=Wd_Est,
        offset_m=Offset_Est,
        u_T_ajuste_s=float(2.0 * math.pi * Errores[2] / (Wd_Est ** 2)),
        u_b_ajuste_kg_s=float(2.0 * Masa * Errores[1]),
        u_A_ajuste_m=float(Errores[0]),
        RMSE_m=float(np.sqrt(np.mean(Residuos ** 2)))
    )


def Verificar(Configuracion):
    Parametros = Configuracion["fisica"]
    Filas = []
    for Mat in Configuracion["materiales"]:
        for Zeta_Val in Configuracion["barrido_zeta"]:
            Masa = Parametros["m_val"]
            Rigidez_K = Mat["k_N_m"]
            Amort_B = 2.0 * Zeta_Val * math.sqrt(Rigidez_K * Masa)
            T_Vec = np.arange(round(Parametros["duracion"] / Parametros["dt"]) + 1, dtype=float) * Parametros["dt"]
            T_Vec, Pos_Num, Vel_Num, Acel_Num, Energ_Num = Integrar(
                Masa, Rigidez_K, Amort_B, Parametros["A0"], T_Vec,
                Parametros["rtol"], Parametros["atol"], Parametros["max_step"]
            )
            Pos_Ana, Vel_Ana = Analitica(T_Vec, Masa, Rigidez_K, Amort_B, Parametros["A0"])
            Err_X = float(np.max(np.abs(Pos_Num - Pos_Ana)))
            Err_V = float(np.max(np.abs(Vel_Num - Vel_Ana)))
            Escala = max(float(Energ_Num[0]), 1e-30)
            Deriva = float(np.max(np.abs(Energ_Num - Energ_Num[0])) / Escala)
            Aumento = float(max(0, np.max(np.diff(Energ_Num))) / Escala)
            Condiciones = (
                abs(Pos_Num[0] - Parametros["A0"]) < 1e-14 and
                abs(Vel_Num[0]) < 1e-14 and
                abs(Acel_Num[0] + Rigidez_K * Parametros["A0"] / Masa) < 1e-12
            )
            Energia_Correcta = Deriva < 1e-7 if Amort_B == 0 else Aumento < 1e-8
            Filas.append(dict(
                Caso=Mat["codigo"],
                zeta=Zeta_Val,
                error_max_x_m=Err_X,
                error_max_v_m_s=Err_V,
                deriva_relativa_E=Deriva,
                aumento_relativo_E=Aumento,
                condiciones_iniciales=bool(Condiciones),
                energia_correcta=bool(Energia_Correcta),
                aprobado=bool(Condiciones and Energia_Correcta and Err_X < 1e-8 and Err_V < 1e-7)
            ))
    return pd.DataFrame(Filas)


def Barra_Visual(Posicion_M, Amplitud_M, Ancho=16):
    if Amplitud_M <= 0:
        return "[      |      ]"
    Fraccion = max(-1.0, min(1.0, Posicion_M / Amplitud_M))
    Centro = Ancho // 2
    Pos_Indice = int(round(Centro + Fraccion * (Centro - 1)))
    Pos_Indice = max(0, min(Ancho - 1, Pos_Indice))
    Caracteres = [" "] * Ancho
    Caracteres[Centro] = "|"
    if Pos_Indice > Centro:
        for I in range(Centro + 1, Pos_Indice + 1):
            Caracteres[I] = ">"
    elif Pos_Indice < Centro:
        for I in range(Pos_Indice, Centro):
            Caracteres[I] = "<"
    else:
        Caracteres[Centro] = "O"
    return "[" + "".join(Caracteres) + "]"


def Mostrar_Tabla_Simulacion(T_Vec, X_Vec, V_Vec, A_Vec, E_Vec, Amplitud_M, Total_Muestras=16):
    Longitud = len(T_Vec)
    Salto = max(1, Longitud // Total_Muestras)
    Indices = list(range(0, Longitud, Salto))
    if Indices[-1] != Longitud - 1:
        Indices.append(Longitud - 1)

    print("\n" + "=" * 76)
    print("EVOLUCION TEMPORAL DEL MOVIMIENTO (Muestreo representativo)")
    print("=" * 76)
    print(f"{'Tiempo (s)':>10} | {'Posicion (cm)':>14} | {'Velocidad (m/s)':>15} | {'Acel (m/s2)':>12} | {'Visual':^18}")
    print("-" * 76)
    for Idx in Indices:
        T_Val = T_Vec[Idx]
        X_Cm = X_Vec[Idx] * 100.0
        V_Val = V_Vec[Idx]
        A_Val = A_Vec[Idx]
        Grafica_Texto = Barra_Visual(X_Vec[Idx], Amplitud_M)
        print(f"{T_Val:>10.3f} | {X_Cm:>+14.4f} | {V_Val:>+15.4f} | {A_Val:>+12.4f} | {Grafica_Texto:^18}")
    print("=" * 76)
    print(f"Total de pasos calculados: {Longitud} puntos temporales.")


def Preguntar_Exportacion(Dataframe_Datos, Nombre_Sugerido="simulacion_resorte"):
    print("\nDesea exportar los resultados?")
    print("1. Exportar a Excel (.xlsx)")
    print("2. Exportar a CSV (.csv)")
    print("3. Exportar ambos (Excel y CSV)")
    print("4. No exportar")
    try:
        Opcion_Exp = int(input("Seleccione una opcion (1-4): "))
    except ValueError:
        print("Opcion no numerica. No se exportara.")
        return

    if Opcion_Exp not in (1, 2, 3):
        return

    Carpeta_Salida = BASE / "resultados"
    Carpeta_Salida.mkdir(parents=True, exist_ok=True)
    Marca_Tiempo = datetime.now().strftime("%Y%m%d_%H%M%S")
    Ruta_Base = Carpeta_Salida / f"{Nombre_Sugerido}_{Marca_Tiempo}"

    if Opcion_Exp in (1, 3):
        Ruta_Excel = Ruta_Base.with_suffix(".xlsx")
        Dataframe_Datos.to_excel(Ruta_Excel, index=False, engine="openpyxl")
        print(f"Archivo Excel guardado en: {Ruta_Excel}")

    if Opcion_Exp in (2, 3):
        Ruta_Csv = Ruta_Base.with_suffix(".csv")
        Dataframe_Datos.to_csv(Ruta_Csv, index=False, encoding="utf-8-sig")
        print(f"Archivo CSV guardado en: {Ruta_Csv}")


def Simular_Resorte_Individual(Configuracion):
    Materiales = Configuracion["materiales"]
    Fisica = Configuracion["fisica"]

    print("\n" + "=" * 60)
    print("SELECCION DE MATERIAL PARA EL RESORTE")
    print("=" * 60)
    for I, Mat in enumerate(Materiales, 1):
        print(f"{I}. {Mat['nombre']} (k = {Mat['k_N_m']:.2f} N/m)")
    print("4. Material personalizado (Definir modulo de corte G)")

    try:
        Opcion_Mat = int(input("\nSeleccione el material (1-4): "))
    except ValueError:
        print("Entrada invalida. Se utilizara Titanio por defecto.")
        Opcion_Mat = 3

    if Opcion_Mat in (1, 2, 3):
        Mat_Elegido = Materiales[Opcion_Mat - 1]
        Nombre_Mat = Mat_Elegido["nombre"]
        Rigidez_K = Mat_Elegido["k_N_m"]
    else:
        Nombre_Mat = "Personalizado"
        print("\nIngrese el modulo de corte G en GPa (Ejemplo: 41.4 para Titanio, 79.6 para Acero):")
        try:
            G_Gpa = float(input("G en GPa: "))
            Rigidez_K = Rigidez(G_Gpa * 1e9, Fisica["d"], Fisica["D"], Fisica["n"])
        except ValueError:
            print("Valor invalido. Se usa k = 20.0 N/m.")
            Rigidez_K = 20.0

    print("\n" + "-" * 60)
    print("PARAMETROS DE ENTRADA")
    print("-" * 60)
    print("Seleccione la masa total suspendidida:")
    print("1. 150 gramos (0.15 kg)")
    print("2. 200 gramos (0.20 kg)")
    print("3. 250 gramos (0.25 kg)")
    print("4. 300 gramos (0.30 kg - Masa de validacion)")
    print("5. 350 gramos (0.35 kg)")
    print("6. Ingresar otro valor manual en gramos")

    try:
        Opcion_Masa = int(input("Seleccione masa (1-6): "))
    except ValueError:
        Opcion_Masa = 4

    Mapeo_Masas = {1: 0.15, 2: 0.20, 3: 0.25, 4: 0.30, 5: 0.35}
    if Opcion_Masa in Mapeo_Masas:
        Masa_Total = Mapeo_Masas[Opcion_Masa]
    else:
        try:
            Masa_Gramo = float(input("Ingrese masa en gramos (ejemplo 280): "))
            Masa_Total = max(0.05, Masa_Gramo / 1000.0)
        except ValueError:
            Masa_Total = 0.30

    print("\nSeleccione el amortiguamiento:")
    print("1. Sin amortiguamiento (zeta = 0.0)")
    print("2. Muy leve / Aire (zeta = 0.02)")
    print("3. Amortiguamiento medio (zeta = 0.10)")
    print("4. Amortiguamiento alto (zeta = 0.30)")
    print("5. Amortiguamiento critico (zeta = 1.00)")
    print("6. Sobreamortiguado (zeta = 1.50)")
    print("7. Ingresar valor personalizado de zeta")

    try:
        Opcion_Zeta = int(input("Seleccione amortiguamiento (1-7): "))
    except ValueError:
        Opcion_Zeta = 2

    Mapeo_Zeta = {1: 0.0, 2: 0.02, 3: 0.10, 4: 0.30, 5: 1.0, 6: 1.5}
    if Opcion_Zeta in Mapeo_Zeta:
        Factor_Zeta = Mapeo_Zeta[Opcion_Zeta]
    else:
        try:
            Factor_Zeta = float(input("Ingrese zeta (ejemplo 0.05): "))
            Factor_Zeta = max(0.0, Factor_Zeta)
        except ValueError:
            Factor_Zeta = 0.02

    print("\nIngrese la amplitud inicial en centimetros (Ejemplo: 1 para 1 cm = 0.01 m):")
    try:
        Amp_Cm = float(input("Amplitud inicial en cm [Enter para 1.0]: ") or "1.0")
        Amplitud_M = max(0.001, Amp_Cm / 100.0)
    except ValueError:
        Amplitud_M = 0.010

    print("\nIngrese la duracion de la simulacion en segundos [Enter para 10]:")
    try:
        Duracion_S = float(input("Duracion en segundos: ") or "10.0")
        Duracion_S = max(1.0, min(100.0, Duracion_S))
    except ValueError:
        Duracion_S = 10.0

    Amortiguamiento_B = 2.0 * Factor_Zeta * math.sqrt(Rigidez_K * Masa_Total)
    Dt = 0.005
    Vector_Tiempo = np.arange(round(Duracion_S / Dt) + 1, dtype=float) * Dt

    T_Vec, X_Vec, V_Vec, A_Vec, E_Vec = Integrar(
        Masa_Total, Rigidez_K, Amortiguamiento_B, Amplitud_M, Vector_Tiempo
    )

    W0 = math.sqrt(Rigidez_K / Masa_Total)
    T_Periodo = Periodo(Masa_Total, Rigidez_K, Amortiguamiento_B)
    Elongacion_Eq_Cm = 100.0 * (Masa_Total * Fisica["g"] / Rigidez_K)
    Elongacion_Max_Cm = Elongacion_Eq_Cm + (Amplitud_M * 100.0)
    Limite_Espacio_Cm = Fisica["limite"] * 100.0
    Holgura_Cm = Limite_Espacio_Cm - Elongacion_Max_Cm
    Aceleracion_Max = float(np.max(np.abs(A_Vec)))
    Tiempo_5 = Tiempo_Al_5(Masa_Total, Rigidez_K, Amortiguamiento_B)

    print("\n" + "=" * 65)
    print(f"RESULTADOS DE LA SIMULACION - {Nombre_Mat.upper()}")
    print("=" * 65)
    print(f"Rigidez calculada (k)           : {Rigidez_K:10.3f} N/m")
    print(f"Masa total suspendida (m)       : {Masa_Total:10.3f} kg ({Masa_Total * 1000.0:.0f} g)")
    print(f"Amortiguamiento viscoso (b)     : {Amortiguamiento_B:10.4f} kg/s")
    print(f"Factor de amortiguamiento (zeta): {Factor_Zeta:10.4f}")
    print(f"Regimen de movimiento           : {Regimen(Factor_Zeta)}")
    print(f"Frecuencia angular propia (w0)  : {W0:10.4f} rad/s")
    if np.isfinite(T_Periodo):
        print(f"Periodo de oscilacion (T)       : {T_Periodo:10.4f} s")
        print(f"Frecuencia natural (f)          : {1.0 / T_Periodo:10.4f} Hz")
    else:
        print("Periodo de oscilacion (T)       : No oscilatorio (Sin periodo)")
    print(f"Elongacion en equilibrio        : {Elongacion_Eq_Cm:10.2f} cm")
    print(f"Elongacion maxima alcanzada     : {Elongacion_Max_Cm:10.2f} cm")
    print(f"Limite maximo disponible        : {Limite_Espacio_Cm:10.2f} cm")
    if Holgura_Cm >= 0:
        print(f"Holgura de seguridad            : {Holgura_Cm:10.2f} cm [CUMPLE: Cabe]")
    else:
        print(f"Holgura de seguridad            : {Holgura_Cm:10.2f} cm [ALERTA: Excede limite]")
    print(f"Aceleracion maxima              : {Aceleracion_Max:10.4f} m/s2")
    if np.isfinite(Tiempo_5):
        print(f"Tiempo de atenuacion al 5%      : {Tiempo_5:10.2f} s")
    else:
        print("Tiempo de atenuacion al 5%      : Infinito (Sin amortiguamiento)")
    print("=" * 65)

    Mostrar_Tabla_Simulacion(T_Vec, X_Vec, V_Vec, A_Vec, E_Vec, Amplitud_M)

    Df_Export = pd.DataFrame({
        "Tiempo_s": T_Vec,
        "Posicion_m": X_Vec,
        "Posicion_cm": X_Vec * 100.0,
        "Velocidad_m_s": V_Vec,
        "Aceleracion_m_s2": A_Vec,
        "Energia_J": E_Vec
    })

    Preguntar_Exportacion(Df_Export, f"simulacion_{Nombre_Mat.replace(' ', '_')}")


def Comparar_Todos_Los_Materiales(Configuracion):
    Materiales = Configuracion["materiales"]
    Fisica = Configuracion["fisica"]
    Masa_Val = Fisica["m_val"]
    Amplitud_M = Fisica["A0"]

    print("\n" + "=" * 90)
    print(f"COMPARACION DE MATERIALES (Masa = {Masa_Val:.2f} kg, Amplitud inicial = {Amplitud_M * 100:.1f} cm)")
    print("=" * 90)
    print(f"{'Material':<20} | {'k (N/m)':>8} | {'T (s)':>8} | {'Acel Max(m/s2)':>14} | {'Eq.(cm)':>8} | {'Holgura':>9} | {'Estado':<10}")
    print("-" * 90)

    Filas_Comparacion = []
    for Mat in Materiales:
        K_Val = Mat["k_N_m"]
        B_Val = Mat["b_ref_kg_s"]
        T_Per = Periodo(Masa_Val, K_Val, B_Val)
        Eq_Cm = 100.0 * (Masa_Val * Fisica["g"] / K_Val)
        Holgura_Cm = 100.0 * (Fisica["limite"] - (Masa_Val * Fisica["g"] / K_Val + Amplitud_M))
        T_Vec = np.arange(round(5.0 / Fisica["dt"]) + 1, dtype=float) * Fisica["dt"]
        _, _, _, A_Vec, _ = Integrar(Masa_Val, K_Val, B_Val, Amplitud_M, T_Vec)
        Acel_Max = float(np.max(np.abs(A_Vec)))
        Estado_Texto = "Cumple" if Holgura_Cm > 0 else "Excede"

        print(f"{Mat['nombre']:<20} | {K_Val:>8.2f} | {T_Per:>8.4f} | {Acel_Max:>14.4f} | {Eq_Cm:>8.2f} | {Holgura_Cm:>8.2f}cm | {Estado_Texto:<10}")

        Filas_Comparacion.append({
            "Material": Mat["nombre"],
            "Codigo": Mat["codigo"],
            "Rigidez_N_m": K_Val,
            "Periodo_s": T_Per,
            "Aceleracion_Max_m_s2": Acel_Max,
            "Elongacion_Eq_cm": Eq_Cm,
            "Holgura_cm": Holgura_Cm,
            "Cabe_En_Recorrido": Holgura_Cm > 0
        })

    print("=" * 90)
    print("CONCLUSION:")
    print("  El Titanio B120VCA presenta la menor aceleracion transmitida (0.552 m/s2),")
    print("  reduciendo las sacudidas en casi un 48% respecto al Acero, y mantiene")
    print("  una holgura positiva de seguridad dentro de los 22 cm asignados.")
    print("=" * 90)

    Df_Comp = pd.DataFrame(Filas_Comparacion)
    Preguntar_Exportacion(Df_Comp, "comparacion_materiales")


def Ejecutar_Menu_Verificacion(Configuracion):
    print("\nEjecutando verificacion numerica contra soluciones analiticas...")
    Df_Verif = Verificar(Configuracion)
    Total_Casos = len(Df_Verif)
    Aprobados = int(Df_Verif["aprobado"].sum())
    Max_Error_X = float(Df_Verif["error_max_x_m"].max())
    Max_Error_V = float(Df_Verif["error_max_v_m_s"].max())

    print("\n" + "=" * 70)
    print("REPORTE DE VERIFICACION NUMERICA")
    print("=" * 70)
    print(f"Casos verificados                : {Total_Casos}")
    print(f"Casos aprobados                  : {Aprobados} / {Total_Casos}")
    print(f"Error maximo en posicion (x)     : {Max_Error_X:.3e} m")
    print(f"Error maximo en velocidad (v)    : {Max_Error_V:.3e} m/s")
    print(f"Resultado final                  : {'APROBADO' if Aprobados == Total_Casos else 'REVISAR'}")
    print("=" * 70)

    Preguntar_Exportacion(Df_Verif, "verificacion_numerica")


def Menu_Principal():
    Configuracion = Cargar_Configuracion()

    while True:
        print("\n" + "=" * 50)
        print("     SIMULADOR DE RESORTES - FISICA 3")
        print("=" * 50)
        print("1. Simular resorte individual (Personalizado)")
        print("2. Comparar los tres materiales (AC, IN, TI)")
        print("3. Ejecutar verificacion numerica del modelo")
        print("4. Cambiar semilla aleatoria")
        print("5. Salir")
        print("=" * 50)

        try:
            Opcion = int(input("Seleccione una opcion (1-5): "))
        except ValueError:
            print("Por favor ingrese un numero entero valido.")
            continue

        if Opcion == 1:
            Simular_Resorte_Individual(Configuracion)
        elif Opcion == 2:
            Comparar_Todos_Los_Materiales(Configuracion)
        elif Opcion == 3:
            Ejecutar_Menu_Verificacion(Configuracion)
        elif Opcion == 4:
            try:
                Nueva_Semilla = int(input("Ingrese nueva semilla (entero no negativo): "))
                if Nueva_Semilla >= 0:
                    Configuracion["fisica"]["semilla"] = Nueva_Semilla
                    print(f"Semilla actualizada a: {Nueva_Semilla}")
                else:
                    print("La semilla debe ser positiva.")
            except ValueError:
                print("Valor invalido.")
        elif Opcion == 5:
            print("\nGracias por usar el simulador. Hasta pronto.")
            break
        else:
            print("Opcion no valida. Seleccione un numero entre 1 y 5.")


def Main():
    Menu_Principal()


if __name__ == "__main__":
    Main()
