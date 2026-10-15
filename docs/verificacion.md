# Pruebas realizadas

[Volver al inicio](../README.md)

El 15 de octubre de 2026 se completó la ejecución incluida en [resultados_ejemplo/](../resultados_ejemplo):

| Comprobación | Resultado |
| --- | --- |
| Pruebas de `pruebas.py` | 9 de 9 aprobadas |
| Comparación con soluciones analíticas | 18 de 18 aprobadas |
| Error máximo de posición | 8,1612 × 10⁻¹² m |
| Error máximo de velocidad | 8,4891 × 10⁻¹¹ m/s |
| Exportación de tablas, gráficas y GIF | Completada |

Las pruebas cubren rigidez, oscilación, energía, amortiguamiento, ajustes, repetibilidad y entradas inválidas. Se exige un error de posición menor de 10⁻⁸ m y de velocidad menor de 10⁻⁷ m/s.

Los controles de material, masa, amortiguamiento, amplitud, pausa y reinicio se probaron con Matplotlib en modo sin pantalla. No se comprobó la apertura de una ventana nativa en Windows o macOS.

El entorno fue Python 3.12.14 en Linux. Las versiones de las dependencias están en [requirements-verificado.txt](../requirements-verificado.txt). Para repetir las pruebas, usa `verificar.bat` en Windows o `.venv/bin/python pruebas.py` en Linux y macOS.

Los errores detallados están en [verificacion_numerica.csv](../resultados_ejemplo/verificacion_numerica.csv). La configuración y la huella del código de esa ejecución se conservan en [parametros_ejecucion.json](../resultados_ejemplo/parametros_ejecucion.json); corresponden a la [versión anterior a esta reorganización](https://github.com/Stivencaballero0710/simulador-resortes-fisica-3/tree/0e64f0bf37da20c2cdd99ef0dd1e9628f94f31e3).

Estas comprobaciones verifican el cálculo del modelo. No sustituyen una validación con mediciones de laboratorio.
