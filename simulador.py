#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Suspensión vertical idealizada: simulación, ajuste, gráficas e interfaz.

Datos computacionales. No son mediciones con CASSY-Lab.
Uso: python simulador.py
     python simulador.py --interactivo
     python simulador.py --verificar
Guía de uso: docs/guia.md.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit

BASE = Path(__file__).resolve().parent
VERSION = "1.0.0"


def cargar_configuracion(path=BASE / "parametros.json", semilla=None):
    cfg = json.loads(Path(path).read_text(encoding="utf-8"))
    p = cfg["fisica"]
    positivos = ("g", "d", "D", "n", "m_ref", "m_val", "limite",
                 "duracion", "dt", "paso_cuantizacion", "rtol", "atol", "max_step")
    for nombre in positivos:
        if not np.isfinite(p[nombre]) or p[nombre] <= 0:
            raise ValueError(f"{nombre} debe ser finito y positivo.")
    for nombre in ("A0", "m_accesorio", "sigma", "zeta_ref"):
        if not np.isfinite(p[nombre]) or p[nombre] < 0:
            raise ValueError(f"{nombre} debe ser finito y no negativo.")
    if p["A0"] <= 0 or not 0 <= p["zeta_ref"] < 1:
        raise ValueError("El protocolo de ajuste requiere A0 > 0 y zeta_ref < 1.")
    if p["m_accesorio"] >= min(p["m_ref"], p["m_val"]):
        raise ValueError("Las masas totales deben superar la masa del accesorio.")
    if p["dt"] > p["duracion"] or not np.isclose(
            p["duracion"] / p["dt"], round(p["duracion"] / p["dt"]),
            rtol=0, atol=1e-8):
        raise ValueError("duracion/dt debe ser un entero positivo.")
    if semilla is not None:
        p["semilla"] = int(semilla)
    if not isinstance(p["semilla"], int) or p["semilla"] < 0:
        raise ValueError("La semilla debe ser un entero no negativo.")
    if not isinstance(cfg["repeticiones"], int) or cfg["repeticiones"] < 3:
        raise ValueError("Se requieren al menos tres repeticiones.")
    if len(cfg["masas_estaticas"]) < 3 or any(
            not np.isfinite(m) or m <= 0 for m in cfg["masas_estaticas"]):
        raise ValueError("Se requieren al menos tres masas positivas.")
    if len(set(cfg["masas_estaticas"])) != len(cfg["masas_estaticas"]):
        raise ValueError("Las masas de caracterización deben ser distintas.")
    if any(not np.isfinite(z) or z < 0 for z in cfg["barrido_zeta"]):
        raise ValueError("Los valores de zeta deben ser finitos y no negativos.")
    if any(not 0 <= cfg["sensibilidad"][s] < 1 for s in ("G", "d", "D")):
        raise ValueError("Las variaciones relativas deben pertenecer a [0,1).")
    codigos = [m["codigo"] for m in cfg["materiales"]]
    if sorted(codigos) != ["AC", "IN", "TI"]:
        raise ValueError("El protocolo requiere los códigos AC, IN y TI.")
    for mat in cfg["materiales"]:
        if "E_Pa" in mat:
            if not -1 < mat["nu"] < 0.5:
                raise ValueError("El coeficiente de Poisson debe estar entre -1 y 0.5.")
            mat["G_Pa"] = mat["E_Pa"] / (2 * (1 + mat["nu"]))
        if not np.isfinite(mat["G_Pa"]) or mat["G_Pa"] <= 0:
            raise ValueError("El módulo de corte debe ser positivo.")
        mat["k_N_m"] = rigidez(mat["G_Pa"], p["d"], p["D"], p["n"])
        mat["b_ref_kg_s"] = 2 * p["zeta_ref"] * math.sqrt(mat["k_N_m"] * p["m_ref"])
    return cfg


def rigidez(G, d, D, n):
    """Rigidez de una hélice circular idealizada, en N/m."""
    return G * d**4 / (8 * D**3 * n)


def malla(p):
    return np.arange(round(p["duracion"] / p["dt"]) + 1, dtype=float) * p["dt"]


def validar_sistema(m, k, b, amplitud):
    if not all(np.isfinite(v) for v in (m, k, b, amplitud)):
        raise ValueError("Los parámetros físicos deben ser finitos.")
    if m <= 0 or k <= 0 or b < 0 or amplitud < 0:
        raise ValueError("m y k positivos; b y amplitud no negativos.")


def integrar(m, k, b, amplitud=0.010, t=None,
             rtol=1e-10, atol=1e-12, max_step=0.03):
    """DOP853; x se mide desde equilibrio y v(0)=0. Devuelve t,x,v,a,E."""
    validar_sistema(m, k, b, amplitud)
    if t is None:
        t = np.arange(3001, dtype=float) * 0.005
    t = np.asarray(t, dtype=float)
    if (t.ndim != 1 or len(t) < 2 or not np.all(np.isfinite(t))
            or t[0] != 0 or np.any(np.diff(t) <= 0)):
        raise ValueError("t debe empezar en 0 y ser finito, unidimensional y creciente.")
    if not all(np.isfinite(v) and v > 0 for v in (rtol, atol, max_step)):
        raise ValueError("Las tolerancias y max_step deben ser positivos.")
    sol = solve_ivp(
        lambda _, y: [y[1], -(b*y[1] + k*y[0]) / m],
        (0, float(t[-1])), [amplitud, 0.0], t_eval=t,
        method="DOP853", rtol=rtol, atol=atol, max_step=max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    x, v = sol.y
    a = -(b*v + k*x) / m
    energia = 0.5*m*v*v + 0.5*k*x*x
    return t, x, v, a, energia


def resolver(p, m, k, b, amplitud=None, t=None):
    return integrar(m, k, b, p["A0"] if amplitud is None else amplitud,
                    malla(p) if t is None else t,
                    p["rtol"], p["atol"], p["max_step"])


def analitica(t, m, k, b, amplitud=0.010):
    """Soluciones independientes de solve_ivp para los tres regímenes."""
    validar_sistema(m, k, b, amplitud)
    t = np.asarray(t, dtype=float)
    if np.any(~np.isfinite(t)) or np.any(t < 0):
        raise ValueError("Los tiempos analíticos deben ser finitos y no negativos.")
    w0 = math.sqrt(k / m)
    gamma = b / (2*m)
    zeta = gamma / w0
    if abs(zeta - 1) <= 1e-12:
        e = np.exp(-w0*t)
        x = amplitud * (1 + w0*t) * e
        v = -amplitud * w0*w0 * t * e
    elif zeta < 1:
        wd = math.sqrt(w0*w0 - gamma*gamma)
        e = np.exp(-gamma*t)
        x = amplitud * e * (np.cos(wd*t) + gamma/wd*np.sin(wd*t))
        v = -amplitud * w0*w0/wd * e * np.sin(wd*t)
    else:
        q = math.sqrt(zeta*zeta - 1)
        r1 = -w0 / (zeta + q)  # evita cancelar dos números próximos
        r2 = -w0 * (zeta + q)
        c1 = -amplitud*r2 / (r1-r2)
        c2 = amplitud*r1 / (r1-r2)
        x = c1*np.exp(r1*t) + c2*np.exp(r2*t)
        v = c1*r1*np.exp(r1*t) + c2*r2*np.exp(r2*t)
    return x, v


def periodo(m, k, b):
    wd2 = k/m - (b/(2*m))**2
    return 2*math.pi/math.sqrt(wd2) if wd2 > 0 else np.nan


def regimen(zeta):
    if zeta == 0:
        return "Sin amortiguar"
    if abs(zeta - 1) <= 1e-12:
        return "Crítico"
    return "Subamortiguado" if zeta < 1 else "Sobreamortiguado"


def tiempo_al_5(m, k, b):
    """Cota de la amplitud exponencial (zeta<1) o cruce monótono (zeta>=1)."""
    w0 = math.sqrt(k/m)
    zeta = b/(2*math.sqrt(k*m))
    if zeta == 0:
        return math.inf
    if zeta < 1 - 1e-12:
        return math.log(1/(0.05*math.sqrt(1-zeta*zeta))) / (zeta*w0)
    def cruce(t):
        return float(analitica(np.array([t]), m, k, b, 1.0)[0][0]) - 0.05
    superior = 10/w0
    while cruce(superior) > 0:
        superior *= 2
        if superior > 1e12:
            raise RuntimeError("No se pudo acotar el cruce del 5 %.")
    return float(brentq(cruce, 0, superior, xtol=1e-12))


def senal_sintetica(limpia, rng, p):
    """Ruido gaussiano seguido de cuantización; no son lecturas instrumentales."""
    q = p["paso_cuantizacion"]
    return np.round((np.asarray(limpia) + rng.normal(
        0, p["sigma"], size=np.shape(limpia))) / q) * q


def ajustar_estatica(x, fuerza):
    """OLS de F=k*x+c; errores estándar descriptivos de la regresión."""
    coef, cov = np.polyfit(np.asarray(x), np.asarray(fuerza), 1, cov=True)
    k, c = map(float, coef)
    if k <= 0:
        raise RuntimeError("La caracterización produjo una rigidez no positiva.")
    resid = np.asarray(fuerza) - (k*np.asarray(x) + c)
    sst = np.sum((fuerza - np.mean(fuerza))**2)
    return dict(k_ajustado_N_m=k, c_N=c,
                u_k_regresion_N_m=float(math.sqrt(max(cov[0, 0], 0))),
                u_c_regresion_N=float(math.sqrt(max(cov[1, 1], 0))),
                R2=float(1-np.sum(resid**2)/sst))


def modelo_ajuste(t, A, gamma, wd, offset):
    return A*np.exp(-gamma*t)*(np.cos(wd*t) + gamma/wd*np.sin(wd*t)) + offset


def ajustar_dinamica(t, x, masa):
    """Estima A, gamma, wd y offset desde la señal sintética subamortiguada."""
    dt = float(t[1]-t[0])
    espectro = np.abs(np.fft.rfft(x-np.mean(x)))
    frecuencias = np.fft.rfftfreq(len(t), dt)
    wd0 = 2*math.pi*frecuencias[1 + int(np.argmax(espectro[1:]))]
    A0 = max(float(np.max(np.abs(x))), 1e-8)
    popt, pcov = curve_fit(
        modelo_ajuste, t, x,
        p0=[A0, 0.02*wd0, wd0, 0],
        bounds=([0, 0, 0.5*wd0, -2*A0],
                [5*A0, 2*wd0, 1.5*wd0, 2*A0]),
        method="trf", x_scale="jac", max_nfev=20000,
        ftol=1e-12, xtol=1e-12, gtol=1e-12)
    if not np.all(np.isfinite(popt)) or not np.all(np.isfinite(pcov)):
        raise RuntimeError("El ajuste dinámico no produjo parámetros identificables.")
    A, gamma, wd, offset = map(float, popt)
    errores = np.sqrt(np.maximum(np.diag(pcov), 0))
    residuos = x-modelo_ajuste(t, *popt)
    return dict(T_s=2*math.pi/wd, b_kg_s=2*masa*gamma,
                A_m=A, gamma_s_1=gamma, wd_rad_s=wd, offset_m=offset,
                u_T_ajuste_s=float(2*math.pi*errores[2]/wd**2),
                u_b_ajuste_kg_s=float(2*masa*errores[1]),
                u_A_ajuste_m=float(errores[0]),
                RMSE_m=float(np.sqrt(np.mean(residuos**2))))


def exportar_csv(df, ruta):
    """Punto decimal para programación y copia con coma decimal para hojas de cálculo."""
    ruta = Path(ruta)
    df.to_csv(ruta, index=False, encoding="utf-8-sig", float_format="%.12g")
    df.to_csv(ruta.with_name(ruta.stem+"_es.csv"), index=False, sep=";",
              decimal=",", encoding="utf-8-sig", float_format="%.12g")


def dataframe_serie(sol, sintetica=None):
    t, x, v, a, E = sol
    df = pd.DataFrame({"t_s":t, "x_limpia_m":x, "v_m_s":v,
                       "a_m_s2":a, "energia_J":E})
    if sintetica is not None:
        df.insert(2, "x_sintetica_m", sintetica)
    return df


def verificar(cfg):
    """18 casos: solución analítica, condiciones iniciales y balance de energía."""
    p = cfg["fisica"]
    filas = []
    for mat in cfg["materiales"]:
        for z in cfg["barrido_zeta"]:
            m, k = p["m_val"], mat["k_N_m"]
            b = 2*z*math.sqrt(k*m)
            t, x, v, a, E = resolver(p, m, k, b)
            xa, va = analitica(t, m, k, b, p["A0"])
            ex = float(np.max(np.abs(x-xa)))
            ev = float(np.max(np.abs(v-va)))
            escala = max(float(E[0]), 1e-30)
            deriva = float(np.max(np.abs(E-E[0])) / escala)
            aumento = float(max(0, np.max(np.diff(E))) / escala)
            condiciones = (abs(x[0]-p["A0"]) < 1e-14 and
                           abs(v[0]) < 1e-14 and
                           abs(a[0] + k*p["A0"]/m) < 1e-12)
            energia_ok = deriva < 1e-7 if b == 0 else aumento < 1e-8
            filas.append(dict(Caso=mat["codigo"], zeta=z,
                error_max_x_m=ex, error_max_v_m_s=ev,
                deriva_relativa_E=deriva, aumento_relativo_E=aumento,
                condiciones_iniciales=bool(condiciones),
                energia_correcta=bool(energia_ok),
                aprobado=bool(condiciones and energia_ok and ex < 1e-8 and ev < 1e-7)))
    return pd.DataFrame(filas)


def comparar_informe(tablas, destino):
    """Solo compara con cifras históricas redondeadas; nunca altera los resultados."""
    ref = json.loads((BASE/"referencia_informe.json").read_text(encoding="utf-8"))
    esquemas = {
        2: ["Caso","m_kg","dl_1_m","dl_2_m","dl_3_m"],
        3: ["Registro","m_kg","T_s","b_kg_s","A_m"],
        4: ["Caso","b_referencia_kg_s","T_previsto_s","T_simulado_medio_s","E_T_porcentaje"],
        5: ["Caso","dl_equilibrio_cm","T_nominal_s","a_max_m_s2","holgura_cm"],
        6: ["zeta","Regimen","b_kg_s","T_s","tiempo_5_s"],
        7: ["Caso","k_min_N_m","k_max_N_m","holgura_k_min_cm"]}
    filas = []
    for n, cols in esquemas.items():
        actual = tablas[n]
        for valores in ref[f"tabla{n}"]:
            referencia = dict(zip(cols, valores))
            llave = cols[0]
            if n == 2:
                mask = (actual["Caso"] == referencia["Caso"]) & np.isclose(
                    actual["m_kg"], referencia["m_kg"])
            elif llave == "zeta":
                mask = np.isclose(actual[llave], referencia[llave])
            else:
                mask = actual[llave] == referencia[llave]
            if not np.any(mask):
                continue
            nueva = actual.loc[mask].iloc[0]
            for col, valor in referencia.items():
                if col == llave or valor is None or isinstance(valor, str):
                    continue
                obtenido = float(nueva[col])
                if not np.isfinite(obtenido):
                    continue
                filas.append(dict(Tabla=n, Registro=str(referencia[llave]),
                    magnitud=col, informe=float(valor), recalculado=obtenido,
                    diferencia=obtenido-float(valor)))
    exportar_csv(pd.DataFrame(filas), destino/"comparacion_con_informe.csv")


def generar(cfg, destino):
    p = cfg["fisica"]
    # Un solo flujo PCG64. Primero toda la estática AC, IN, TI;
    # después dinámica AC, IN, TI, REF antes de VAL y repeticiones ascendentes.
    rng = np.random.Generator(np.random.PCG64(p["semilla"]))
    tablas, ajustes, series = {}, [], {}
    estatica = []
    nr = cfg["repeticiones"]
    for mat in cfg["materiales"]:
        filas_mat = []
        for masa in cfg["masas_estaticas"]:
            dl = senal_sintetica(np.repeat(masa*p["g"]/mat["k_N_m"], nr), rng, p)
            fila = dict(Caso=mat["codigo"], m_kg=masa)
            fila.update({f"dl_{i+1}_m":float(v) for i, v in enumerate(dl)})
            fila.update(Media_m=float(np.mean(dl)), mg_N=masa*p["g"],
                        s_dl_m=float(np.std(dl, ddof=1)),
                        u_media_dl_m=float(np.std(dl, ddof=1)/math.sqrt(nr)))
            estatica.append(fila)
            filas_mat.append(fila)
        d = pd.DataFrame(filas_mat)
        ajustes.append(dict(Caso=mat["codigo"],
            **ajustar_estatica(d["Media_m"].to_numpy(), d["mg_N"].to_numpy())))
    tablas[2] = pd.DataFrame(estatica)
    ajuste_estatico = pd.DataFrame(ajustes)

    dinamica = []
    for mat in cfg["materiales"]:
        for clase, masa in (("REF", p["m_ref"]), ("VAL", p["m_val"])):
            sol = resolver(p, masa, mat["k_N_m"], mat["b_ref_kg_s"])
            for r in range(1, nr+1):
                nombre = f'{mat["codigo"]}-{clase}-{r}'
                ruido = senal_sintetica(sol[1], rng, p)
                df = dataframe_serie(sol, ruido)
                series[nombre] = df
                exportar_csv(df, destino/"series"/f"{nombre}.csv")
                dinamica.append(dict(Registro=nombre, Caso=mat["codigo"],
                    Condicion=clase, Repeticion=r, m_kg=masa,
                    **ajustar_dinamica(sol[0], ruido, masa)))
    tablas[3] = pd.DataFrame(dinamica)

    contraste, nominales, predicciones = [], [], {}
    for mat in cfg["materiales"]:
        cod, k, m = mat["codigo"], mat["k_N_m"], p["m_val"]
        kfit = float(ajuste_estatico.loc[ajuste_estatico.Caso == cod, "k_ajustado_N_m"].iloc[0])
        ref = tablas[3].loc[(tablas[3].Caso == cod) & (tablas[3].Condicion == "REF")]
        val = tablas[3].loc[(tablas[3].Caso == cod) & (tablas[3].Condicion == "VAL")]
        bfit = float(ref.b_kg_s.mean())
        Tpred = periodo(m, kfit, bfit)
        Tval = float(val.T_s.mean())
        contraste.append(dict(Caso=cod, b_referencia_kg_s=bfit,
            T_previsto_s=Tpred, T_simulado_medio_s=Tval,
            E_T_porcentaje=abs(Tpred-Tval)/Tval*100,
            s_T_s=float(val.T_s.std(ddof=1)),
            u_T_media_s=float(val.T_s.std(ddof=1)/math.sqrt(nr)),
            s_b_ref_kg_s=float(ref.b_kg_s.std(ddof=1))))
        pred = dataframe_serie(resolver(p, m, kfit, bfit))
        predicciones[cod] = pred
        exportar_csv(pred, destino/"series"/f"{cod}-PRED.csv")
        sol = resolver(p, m, k, mat["b_ref_kg_s"])
        holgura = p["limite"]-(m*p["g"]/k+p["A0"])
        nominales.append(dict(Caso=cod, k_N_m=k,
            dl_equilibrio_cm=100*m*p["g"]/k,
            T_nominal_s=periodo(m, k, mat["b_ref_kg_s"]),
            a_max_m_s2=float(np.max(np.abs(sol[3]))), holgura_cm=100*holgura,
            cabe=holgura > 0,
            zeta_validacion=mat["b_ref_kg_s"]/(2*math.sqrt(k*m))))
    tablas[4], tablas[5] = pd.DataFrame(contraste), pd.DataFrame(nominales)

    ti = next(mat for mat in cfg["materiales"] if mat["codigo"] == "TI")
    barrido, curvas = [], {}
    for zeta in cfg["barrido_zeta"]:
        m, k = p["m_val"], ti["k_N_m"]
        b = 2*zeta*math.sqrt(k*m)
        sol = resolver(p, m, k, b)
        curvas[zeta] = sol
        exportar_csv(dataframe_serie(sol), destino/"series"/f"TI-zeta-{zeta:g}.csv")
        barrido.append(dict(zeta=zeta, Regimen=regimen(zeta), b_kg_s=b,
            T_s=periodo(m, k, b) if zeta < 1 else np.nan,
            tiempo_5_s=tiempo_al_5(m, k, b),
            tipo_tiempo="Cota exponencial" if 0 < zeta < 1 else
                ("No finito" if zeta == 0 else "Cruce monótono"),
            a_max_m_s2=float(np.max(np.abs(sol[3])))))
    tablas[6] = pd.DataFrame(barrido)

    sensibilidad = []
    s = cfg["sensibilidad"]
    for mat in cfg["materiales"]:
        kmin = rigidez(mat["G_Pa"]*(1-s["G"]), p["d"]*(1-s["d"]),
                       p["D"]*(1+s["D"]), p["n"])
        kmax = rigidez(mat["G_Pa"]*(1+s["G"]), p["d"]*(1+s["d"]),
                       p["D"]*(1-s["D"]), p["n"])
        sensibilidad.append(dict(Caso=mat["codigo"], k_min_N_m=kmin, k_max_N_m=kmax,
            holgura_k_min_cm=100*(p["limite"]-(p["m_val"]*p["g"]/kmin+p["A0"]))))
    tablas[7] = pd.DataFrame(sensibilidad)
    tablas[8] = pd.DataFrame([
        ["Recipiente y accesorio","Recuperar un recipiente apropiado para la carga.",
         "Propuesto; sin montaje ni pesaje."],
        ["Sujeciones","Usar uniones desmontables y revisar su estado.",
         "Propuestas; verificación pendiente."],
        ["Resorte y soporte","Caracterizar componentes disponibles antes de seleccionarlos.",
         "Alternativas idealizadas en Python."]],
        columns=["Elemento_previsto","Criterio_de_reutilizacion","Estado"])

    parametros = [(k, str(v), "Parámetro del escenario") for k, v in p.items()]
    parametros += [(f'{m["codigo"]}: G_Pa', str(m["G_Pa"]), m["procedencia"])
                   for m in cfg["materiales"]]
    tablas[1] = pd.DataFrame(parametros, columns=["Parametro","Valor_SI","Procedencia"])
    for n, df in sorted(tablas.items()):
        exportar_csv(df, destino/f"tabla_{n}.csv")
    exportar_csv(ajuste_estatico, destino/"ajustes_estaticos.csv")
    comprobaciones = verificar(cfg)
    exportar_csv(comprobaciones, destino/"verificacion_numerica.csv")
    comparar_informe(tablas, destino)
    graficar(cfg, tablas, series, predicciones, curvas, destino)
    return tablas, comprobaciones


def graficar(cfg, tablas, series, predicciones, curvas, destino):
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10,
                         "axes.grid":True, "grid.alpha":0.20})
    def guardar(fig, nombre):
        fig.tight_layout()
        fig.savefig(destino/"graficas"/f"{nombre}.png", dpi=240)
        fig.savefig(destino/"graficas"/f"{nombre}.pdf")
        plt.close(fig)
    def limitar(ax):
        ax.set_xlim(0, min(5, cfg["fisica"]["duracion"]))

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for mat in cfg["materiales"]:
        df = tablas[2].loc[tablas[2].Caso == mat["codigo"]]
        ax.plot(df.Media_m*100, df.mg_N, "o-", color=mat["color"], label=mat["nombre"])
    ax.set(xlabel="Elongación simulada (cm)", ylabel="Fuerza (N)")
    ax.legend()
    guardar(fig, "figura_1_estatica")

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for mat in cfg["materiales"]:
        df = series[f'{mat["codigo"]}-VAL-1']
        axes[0].plot(df.t_s, df.x_limpia_m*100, color=mat["color"], label=mat["nombre"])
        axes[1].plot(df.t_s, df.a_m_s2, color=mat["color"])
    axes[0].set_ylabel("Posición (cm)")
    axes[0].legend()
    axes[1].set(xlabel="Tiempo (s)", ylabel="Aceleración (m/s²)")
    limitar(axes[1])
    guardar(fig, "figura_2_posicion_aceleracion")

    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    for ax, mat in zip(axes, cfg["materiales"]):
        df = series[f'{mat["codigo"]}-VAL-1']
        pred = predicciones[mat["codigo"]]
        ax.plot(df.t_s, df.x_sintetica_m*100, ".", ms=1.5, alpha=.6,
                color=mat["color"], label="Serie sintética reservada 1")
        ax.plot(pred.t_s, pred.x_limpia_m*100, color="#333333", lw=1,
                label="Predicción desde referencia")
        ax.set(title=mat["nombre"], ylabel="x (cm)")
        limitar(ax)
    axes[0].legend(fontsize=8, ncol=2)
    axes[-1].set_xlabel("Tiempo (s)")
    guardar(fig, "figura_3_contraste")

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
    df = tablas[5].set_index("Caso").loc[[m["codigo"] for m in cfg["materiales"]]]
    labels = ["Acero", "Inoxidable", "Titanio"]
    colores = [m["color"] for m in cfg["materiales"]]
    axes[0].bar(labels, df.a_max_m_s2, color=colores)
    axes[1].bar(labels, df.holgura_cm, color=colores)
    axes[0].set_ylabel("Aceleración máxima (m/s²)")
    axes[1].set_ylabel("Holgura (cm)")
    axes[1].axhline(0, color="#333333", lw=1)
    guardar(fig, "figura_4_comparacion")

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    p = cfg["fisica"]
    for zeta, sol in curvas.items():
        t, x, _, _, E = sol
        axes[0].plot(t, x/p["A0"], label=f"ζ = {zeta:g}")
        axes[1].plot(t, E/E[0])
    axes[0].set_ylabel("x / A₀")
    axes[0].legend(ncol=3, fontsize=8)
    axes[1].set(xlabel="Tiempo (s)", ylabel="E / E₀")
    limitar(axes[1])
    guardar(fig, "figura_5_amortiguamiento")


def crear_interfaz(cfg, destino=None, animar=True):
    """Interfaz local Matplotlib: material, masa total, zeta y separación inicial."""
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from matplotlib.widgets import Button, RadioButtons, Slider

    p = cfg["fisica"]
    mats = {m["codigo"]:m for m in cfg["materiales"]}
    estado = {"codigo":"TI", "indice":0, "pausa":False}
    fig = plt.figure(figsize=(11, 7.3))
    fig.suptitle("Suspensión vertical | Exploración numérica", fontsize=14)
    ax_resorte = fig.add_axes([.08, .36, .20, .54])
    ax_tiempo = fig.add_axes([.37, .42, .59, .48])
    ax_radio = fig.add_axes([.06, .11, .21, .18])
    radio = RadioButtons(ax_radio, ("AC", "IN", "TI"), active=2)
    radio.set_active(2)
    sliders = [
        Slider(fig.add_axes([.40, .29, .52, .03]), "Masa total (kg)",
               max(.09, p["m_accesorio"]+.01), .50, valinit=p["m_val"], valstep=.005),
        Slider(fig.add_axes([.40, .23, .52, .03]), "ζ", 0, 2, valinit=.30, valstep=.01),
        Slider(fig.add_axes([.40, .17, .52, .03]), "A₀ (mm)", 0, 20,
               valinit=p["A0"]*1000, valstep=.25)]
    boton = Button(fig.add_axes([.73, .095, .20, .05]), "Pausa / continuar")
    reset = Button(fig.add_axes([.49, .095, .20, .05]), "Reiniciar tiempo")
    info = fig.text(.36, .015, "", fontsize=9)
    ax_resorte.set(xlim=(-.8, .8), xlabel="Hélice esquemática", ylabel="Elongación (m)")
    ax_resorte.set_xticks([])
    resorte, = ax_resorte.plot([], [], lw=2)
    masa_dibujada, = ax_resorte.plot([], [], "s", ms=18)
    ax_resorte.axhline(p["limite"], color="#a33434", ls="--", label="Recorrido asignado")
    curva, = ax_tiempo.plot([], [], lw=1.5)
    punto, = ax_tiempo.plot([], [], "o", ms=6)
    cursor = ax_tiempo.axvline(0, color="#777777", lw=.8)
    ax_tiempo.set(xlabel="Tiempo (s)", ylabel="Posición respecto al equilibrio (mm)")
    ax_tiempo.grid(alpha=.2)
    fase = np.linspace(0, 1, 401)

    def dibujar_indice(i):
        mat = mats[estado["codigo"]]
        t, x, _, _, _ = estado["sol"]
        i = int(i) % len(t)
        fin = estado["equilibrio"] + x[i]
        resorte.set_data(.24*np.sin(2*np.pi*p["n"]*fase), fase*fin)
        masa_dibujada.set_data([0], [fin])
        punto.set_data([t[i]], [1000*x[i]])
        cursor.set_xdata([t[i], t[i]])
        for artist in (resorte, masa_dibujada, curva, punto):
            artist.set_color(mat["color"])
        return resorte, masa_dibujada, punto, cursor

    def actualizar(_=None):
        mat = mats[estado["codigo"]]
        m, zeta, A_mm = [s.val for s in sliders]
        k, A = mat["k_N_m"], A_mm/1000
        b = 2*zeta*math.sqrt(k*m)
        estado["sol"] = resolver(p, m, k, b, amplitud=A)
        estado["equilibrio"] = m*p["g"]/k
        estado["indice"] = 0
        t, x, _, _, _ = estado["sol"]
        h = p["limite"]-(estado["equilibrio"]+A)
        curva.set_data(t, 1000*x)
        ax_tiempo.set_xlim(0, p["duracion"])
        ax_tiempo.set_ylim(-max(A_mm*1.15, 1), max(A_mm*1.15, 1))
        extremo = max(p["limite"], estado["equilibrio"]+A)*1.18
        ax_resorte.set_ylim(extremo, min(0, estado["equilibrio"]-A)-.01)
        condicion = "Cabe en el recorrido supuesto" if h > 0 else "Excede el recorrido supuesto"
        info.set_text(f'{mat["nombre"]} | k = {k:.3f} N/m | b = {b:.4f} kg/s | {regimen(zeta)}\n'
                      f'Holgura = {h*1000:.1f} mm | {condicion}\n'
                      "Datos simulados. ζ controla el sistema; no es una propiedad medida del metal.")
        dibujar_indice(0)
        fig.canvas.draw_idle()

    def elegir(codigo):
        estado["codigo"] = codigo
        actualizar()

    def pausar(_):
        estado["pausa"] = not estado["pausa"]

    def reiniciar(_):
        estado["indice"] = 0
        dibujar_indice(0)
        fig.canvas.draw_idle()

    radio.on_clicked(elegir)
    for s in sliders:
        s.on_changed(actualizar)
    boton.on_clicked(pausar)
    reset.on_clicked(reiniciar)
    actualizar()
    ani = None
    if animar:
        paso = max(1, round(.04/p["dt"]))
        def avanzar(_):
            if not estado["pausa"]:
                estado["indice"] = (estado["indice"]+paso) % len(estado["sol"][0])
            return dibujar_indice(estado["indice"])
        ani = FuncAnimation(fig, avanzar, interval=40, blit=False,
                            cache_frame_data=False)
    # Se conservan controles y animación mientras vive la figura.
    fig._controles_simulador = (radio, sliders, boton, reset, ani, estado)
    if destino is not None:
        fig.savefig(Path(destino)/"figura_A1_interfaz.png", dpi=180)
    return fig, ani


def exportar_gif(cfg, destino):
    """Animación TI, m reservada y zeta=0.30; 15 s, 20 fotogramas/s por defecto."""
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    p = cfg["fisica"]
    ti = next(m for m in cfg["materiales"] if m["codigo"] == "TI")
    m, k, zeta = p["m_val"], ti["k_N_m"], .30
    b = 2*zeta*math.sqrt(k*m)
    t, x, _, _, _ = resolver(p, m, k, b)
    eq = m*p["g"]/k
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4.7),
                                  gridspec_kw={"width_ratios":[1, 2.5]})
    ax1.set(xlim=(-.7, .7), ylim=(max(p["limite"], eq+p["A0"])*1.15, 0),
            ylabel="Elongación (m)", title="Hélice esquemática")
    ax1.set_xticks([])
    ax1.axhline(p["limite"], color="#a33434", ls="--")
    linea, = ax1.plot([], [], color=ti["color"], lw=2)
    masa, = ax1.plot([], [], "s", color=ti["color"], ms=18)
    ax2.plot(t, 1000*x, color=ti["color"], lw=1)
    punto, = ax2.plot([], [], "o", color=ti["color"])
    ax2.set(xlim=(0, p["duracion"]), xlabel="Tiempo (s)", ylabel="x (mm)")
    leyenda = fig.suptitle("Titanio B120VCA | Datos simulados")
    fig.tight_layout(rect=(0, 0, 1, .94))
    u = np.linspace(0, 1, 401)
    fps = 20
    tiempos_video = np.linspace(0, p["duracion"], max(2, round(p["duracion"]*fps)+1))
    indices = np.minimum(np.round(tiempos_video/p["dt"]).astype(int), len(t)-1)
    def animar(j):
        i = indices[j]
        linea.set_data(.24*np.sin(2*np.pi*p["n"]*u), (eq+x[i])*u)
        masa.set_data([0], [eq+x[i]])
        punto.set_data([t[i]], [1000*x[i]])
        leyenda.set_text(f"Titanio B120VCA | ζ = {zeta:g} | t = {t[i]:.2f} s | Datos simulados")
        return linea, masa, punto, leyenda
    ani = FuncAnimation(fig, animar, frames=len(indices), interval=1000/fps, blit=False)
    ani.save(Path(destino)/"animacion_TI.gif", writer=PillowWriter(fps=fps), dpi=90)
    plt.close(fig)


def versiones():
    return {"Python":platform.python_version(), "sistema":platform.platform(),
            **{p:importlib.metadata.version(p) for p in
               ("numpy", "scipy", "pandas", "matplotlib", "Pillow")}}


def escribir_resumen(cfg, tablas, comprobaciones, destino):
    n = tablas[5]
    validos = n.loc[n.cabe]
    elegido = "Ninguno" if validos.empty else str(
        validos.sort_values("a_max_m_s2").iloc[0].Caso)
    error = float(comprobaciones.error_max_x_m.max())
    reduccion = 100*(1 - float(n.loc[n.Caso == "TI", "a_max_m_s2"].iloc[0]) /
                          float(n.loc[n.Caso == "AC", "a_max_m_s2"].iloc[0]))
    texto = (
        "RESULTADOS DE SIMULACIÓN\n"
        "No se adquirieron datos con CASSY-Lab.\n\n"
        f"Alternativa de menor aceleración con holgura positiva: {elegido}\n"
        f"Reducción nominal TI respecto a AC: {reduccion:.4f} %\n"
        f"Error máximo de posición frente a solución analítica: {error:.6e} m\n"
        f"Comprobaciones aprobadas: {int(comprobaciones.aprobado.sum())}"
        f"/{len(comprobaciones)}\n\n"
        "La verificación corresponde al modelo ideal, no al montaje físico.\n"
        "El archivo comparacion_con_informe.csv contrasta cifras nuevas con las\n"
        "cifras redondeadas del informe. No se fuerzan resultados históricos.\n"
        "Las covarianzas y dispersiones son estadísticas del registro sintético,\n"
        "no incertidumbres de sensores ni intervalos de confianza de materiales.\n"
        "Guía de archivos y unidades: docs/guia.md, en la carpeta del proyecto.\n")
    (destino/"RESUMEN.txt").write_text(texto, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=BASE/"parametros.json")
    ap.add_argument("--salida", type=Path, help="Carpeta nueva para los resultados.")
    ap.add_argument("--semilla", type=int, help="Sobrescribe la semilla del archivo de parámetros.")
    modo = ap.add_mutually_exclusive_group()
    modo.add_argument("--interactivo", action="store_true", help="Abrir controles locales.")
    modo.add_argument("--verificar", action="store_true", help="Solo verificación analítica.")
    ap.add_argument("--sin-gif", action="store_true", help="Omitir la animación en ejecución completa.")
    args = ap.parse_args()
    cfg = cargar_configuracion(args.config, args.semilla)
    import matplotlib
    if not args.interactivo:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if args.interactivo:
        if "agg" == matplotlib.get_backend().lower():
            raise RuntimeError("No hay backend gráfico. Ejecute localmente con Tk disponible.")
        fig, animacion = crear_interfaz(cfg)
        plt.show()
        return 0
    sello = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    destino = (args.salida or BASE/"resultados"/sello).resolve()
    destino.mkdir(parents=True, exist_ok=False)  # no pisa una ejecución anterior
    (destino/"series").mkdir()
    (destino/"graficas").mkdir()
    metadatos = dict(configuracion=cfg, versiones=versiones(),
        fecha_UTC=datetime.now(timezone.utc).isoformat(),
        version_programa=VERSION,
        sha256_programa=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        generador_aleatorio="NumPy Generator(PCG64)",
        orden_aleatorio="Estática AC IN TI; después dinámica por material, REF VAL, repetición.",
        datos="Simulados; sin adquisición experimental.")
    (destino/"parametros_ejecucion.json").write_text(
        json.dumps(metadatos, indent=2, ensure_ascii=False), encoding="utf-8")
    if args.verificar:
        comprobaciones = verificar(cfg)
        exportar_csv(comprobaciones, destino/"verificacion_numerica.csv")
    else:
        tablas, comprobaciones = generar(cfg, destino)
        fig, _ = crear_interfaz(cfg, destino/"graficas", animar=False)
        plt.close(fig)
        if not args.sin_gif:
            exportar_gif(cfg, destino/"graficas")
        escribir_resumen(cfg, tablas, comprobaciones, destino)
    ok = bool(comprobaciones.aprobado.all())
    (destino/"ESTADO.txt").write_text(
        ("VERIFICACION NUMERICA APROBADA" if ok else "REVISAR VERIFICACION NUMERICA")+
        "\nDatos simulados; no validación experimental.\n", encoding="utf-8")
    print(f"Archivos guardados en: {destino}")
    print(f"Verificación numérica: {'APROBADA' if ok else 'REVISAR'}")
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())

