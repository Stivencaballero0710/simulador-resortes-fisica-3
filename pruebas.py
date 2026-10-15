# -*- coding: utf-8 -*-
"""Pruebas de física y de recuperación de parámetros. Ejecutar: python pruebas.py"""
import math
import unittest

import numpy as np

from simulador import (analitica, ajustar_dinamica, cargar_configuracion, integrar,
                       periodo, senal_sintetica, tiempo_al_5, verificar)


class PruebasFisicas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = cargar_configuracion()

    def test_rigideces_nominales_del_informe(self):
        esperadas = {"AC":31.84, "IN":28.46153846153846, "TI":16.56}
        for mat in self.cfg["materiales"]:
            self.assertAlmostEqual(mat["k_N_m"], esperadas[mat["codigo"]], places=10)

    def test_oscilar_sin_amortiguamiento(self):
        m, k, A = .3, 16.56, .01
        T = 2*math.pi*math.sqrt(m/k)
        t = np.array([0., T/4, T/2, 3*T/4, T])
        _, x, v, a, E = integrar(m, k, 0, A, t)
        np.testing.assert_allclose(x, [A, 0, -A, 0, A], atol=1e-10)
        self.assertAlmostEqual(v[1], -A*math.sqrt(k/m), places=8)
        np.testing.assert_allclose(E, E[0], rtol=1e-7)
        self.assertAlmostEqual(a[0], -k*A/m, places=12)

    def test_critico_y_sobreamortiguado(self):
        t = np.linspace(0, 5, 501)
        m, k, A = .3, 16.56, .01
        w0 = math.sqrt(k/m)
        xc, vc = analitica(t, m, k, 2*math.sqrt(k*m), A)
        np.testing.assert_allclose(xc, A*(1+w0*t)*np.exp(-w0*t), atol=1e-14)
        for zeta in (1, 1.5):
            b = 2*zeta*math.sqrt(k*m)
            _, x, v, _, _ = integrar(m, k, b, A, t)
            self.assertTrue(np.all(x >= -1e-10))
            self.assertTrue(np.all(np.diff(x) <= 1e-10))
            self.assertTrue(np.all(v <= 1e-10))

    def test_tiempo_5_cruce_y_cota(self):
        m, k = .3, 16.56
        self.assertTrue(math.isinf(tiempo_al_5(m, k, 0)))
        for zeta in (.02, .1, .3, 1., 1.5):
            b = 2*zeta*math.sqrt(k*m)
            t5 = tiempo_al_5(m, k, b)
            x, _ = analitica(np.linspace(t5, t5+20, 4001), m, k, b, 1)
            self.assertLessEqual(np.max(np.abs(x)), .05+1e-10)
            if zeta >= 1:
                self.assertAlmostEqual(float(x[0]), .05, places=9)

    def test_18_comparaciones_analiticas(self):
        df = verificar(self.cfg)
        self.assertEqual(len(df), 18)
        self.assertTrue(bool(df.aprobado.all()), df.to_string())

    def test_recuperacion_desde_senal_limpia(self):
        m, k, b, A = .3, 16.56, .081388942737259, .01
        t = np.arange(3001)*.005
        x, _ = analitica(t, m, k, b, A)
        p = ajustar_dinamica(t, x, m)
        self.assertAlmostEqual(p["T_s"], periodo(m, k, b), places=6)
        self.assertAlmostEqual(p["b_kg_s"], b, places=6)
        self.assertAlmostEqual(p["A_m"], A, places=7)
        self.assertLess(p["RMSE_m"], 1e-8)

    def test_ruido_reproducible_y_cuantizado(self):
        p = self.cfg["fisica"]
        x = np.linspace(-.01, .01, 1000)
        a = senal_sintetica(x, np.random.Generator(np.random.PCG64(p["semilla"])), p)
        b = senal_sintetica(x, np.random.Generator(np.random.PCG64(p["semilla"])), p)
        np.testing.assert_array_equal(a, b)
        np.testing.assert_allclose(a/p["paso_cuantizacion"],
                                   np.round(a/p["paso_cuantizacion"]), atol=1e-10)

    def test_masa_reservada_y_holgura(self):
        p = self.cfg["fisica"]
        esperadas = {"AC":(.61000, 1.0613, 11.757),
                     "IN":(.64518, .9487, 10.660),
                     "TI":(.84583, .5520, 3.228)}
        for mat in self.cfg["materiales"]:
            k, m = mat["k_N_m"], p["m_val"]
            T, a, h = esperadas[mat["codigo"]]
            self.assertAlmostEqual(periodo(m, k, mat["b_ref_kg_s"]), T, delta=5e-6)
            self.assertAlmostEqual(k*p["A0"]/m, a, delta=5e-5)
            calculada = 100*(p["limite"]-(m*p["g"]/k+p["A0"]))
            self.assertAlmostEqual(calculada, h, delta=.0005)

    def test_entradas_invalidas(self):
        for m, k, b, A in ((0, 10, 0, .01), (.3, -1, 0, .01),
                            (.3, 10, -1, .01), (.3, 10, 0, -1),
                            (.3, math.nan, 0, .01)):
            with self.assertRaises(ValueError):
                integrar(m, k, b, A)
        with self.assertRaises(ValueError):
            integrar(.3, 10, .1, t=[0, .1, .1])
        with self.assertRaises(ValueError):
            integrar(.3, 10, .1, t=[1, 2, 3])


if __name__ == "__main__":
    unittest.main(verbosity=2)

