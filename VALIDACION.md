# Validación de la implementación

Ejecución completada el 15 de octubre de 2026 (UTC). Todos los datos son de simulación.

## Comprobaciones ejecutadas

- `python pruebas.py`: **9 pruebas aprobadas**. Incluyen rigideces, oscilación sin amortiguamiento, regímenes crítico y sobreamortiguado, tiempo al 5 %, comparación analítica, recuperación de parámetros, reproducibilidad, resultados nominales y entradas inválidas.
- `python simulador.py --salida resultados_ejemplo`: **terminó con código 0**. Produjo ocho tablas principales, series, cinco figuras en PNG/PDF, una vista de la interfaz y una animación GIF.
- Las **18/18 comparaciones** analíticas cumplieron las condiciones iniciales, los criterios de energía y los umbrales de error.
- Error máximo de posición: **8,16120174912 × 10⁻¹² m**.
- Error máximo de velocidad: **8,4890599758 × 10⁻¹¹ m/s**.
- Se ejercitaron selección de material, masa, amortiguamiento, amplitud cero, pausa y reinicio mediante los controles de Matplotlib con el backend Agg. Esta comprobación no equivale a abrir una ventana en Windows/macOS.

Los resultados detallados están en [verificacion_numerica.csv](resultados_ejemplo/verificacion_numerica.csv). Los datos de procedencia y la huella SHA-256 del programa se conservan en [parametros_ejecucion.json](resultados_ejemplo/parametros_ejecucion.json).

## Entorno utilizado

| Componente | Versión |
| --- | --- |
| Python | 3.12.14 |
| NumPy | 2.3.5 |
| SciPy | 1.17.0 |
| pandas | 2.2.3 |
| Matplotlib | 3.10.8 |
| Pillow | 12.3.0 |
| Sistema | Linux x86_64, glibc 2.39 |

`requirements-verificado.txt` fija las versiones de las cinco dependencias directas utilizadas. `requirements.txt` mantiene intervalos compatibles para una instalación general.

## Diferencias respecto del informe

El informe conservado en `informe/` mantiene las cifras históricas. Los resultados de `resultados_ejemplo/` se calcularon de nuevo. [comparacion_con_informe.csv](resultados_ejemplo/comparacion_con_informe.csv) muestra las diferencias, sin reemplazar ni ajustar los datos para hacerlos coincidir.

No se dispuso del script original. La igualdad de la semilla no garantiza igualdad de las series si cambia el orden de generación, el algoritmo o el ajuste. Esta comprobación numérica no constituye validación experimental ni certificación de un componente físico.
