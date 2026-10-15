# Guía de uso

[Volver al inicio](../README.md)

## Usar la ventana

Elige un material y mueve los controles:

- **Masa:** masa total suspendida, incluido el accesorio.
- **ζ:** amortiguamiento; indica qué tan rápido se apaga la oscilación.
- **A₀:** separación inicial respecto al equilibrio, en milímetros.

La gráfica muestra la posición a lo largo del tiempo. La holgura indica cuánto recorrido queda disponible. Usa **Pausa / continuar** para detener el movimiento y **Reiniciar tiempo** para volver al principio.

Los cambios de la ventana se usan solo durante esa sesión. Para guardar otra configuración, edita una copia de [parametros.json](../parametros.json), cuyos valores están en unidades del SI.

## Ejecutar desde la terminal

En Linux o macOS, abre una terminal en la carpeta del proyecto:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python simulador.py --interactivo
```

En Windows puedes usar los archivos `.bat` del inicio. Si prefieres la terminal, crea el entorno con `py -3 -m venv .venv` y sustituye `.venv/bin/python` por `.venv\Scripts\python.exe` en los demás comandos.

Para generar resultados sin abrir la ventana:

```bash
.venv/bin/python simulador.py --salida resultados/ensayo_01
```

La carpeta de salida debe ser nueva. Si omites `--salida`, se crea una carpeta con fecha y hora dentro de `resultados/`.

| Opción | Qué hace |
| --- | --- |
| `--interactivo` | Abre la ventana. |
| `--sin-gif` | Genera tablas y gráficas sin la animación. |
| `--config mi_configuracion.json` | Usa tu configuración. |
| `--semilla 2250191` | Fija la secuencia del ruido simulado. |
| `--verificar` | Compara el cálculo numérico con la solución analítica. |

## Leer los resultados

Empieza por `RESUMEN.txt` y la carpeta `graficas/`. Las tablas se organizan así:

| Tabla | Contenido |
| --- | --- |
| 1 | Parámetros utilizados. |
| 2 | Elongaciones estáticas. |
| 3 | Ajustes de las oscilaciones. |
| 4 | Predicción y comparación para la masa de 0,300 kg. |
| 5 | Aceleración máxima y holgura. |
| 6 | Efecto del amortiguamiento. |
| 7 | Sensibilidad a cambios del material y la geometría. |
| 8 | Propuesta de reutilización del informe. |

En `series/` están los registros de tiempo y movimiento. `parametros_ejecucion.json` guarda la configuración, la semilla, las versiones y la huella del código. `comparacion_con_informe.csv` muestra las diferencias respecto a las cifras del informe.

Los archivos `_es.csv` usan punto y coma como separador y coma decimal; suelen abrirse mejor en Excel configurado en español. Los otros CSV usan coma y punto decimal. Las unidades aparecen en los encabezados. Un período vacío significa que no hay oscilación; `inf` indica que el movimiento no se estabiliza dentro del umbral definido.

## Si algo falla

| Problema | Qué hacer |
| --- | --- |
| Falta una biblioteca | Vuelve a ejecutar `instalar.bat`. |
| No se abre la ventana | Comprueba que Python incluya Tcl/Tk. Sin pantalla puedes generar los archivos. |
| La carpeta de salida ya existe | Elige otro nombre o quita `--salida`. |
| La exportación tarda mucho | Añade `--sin-gif`. |

Para revisar la instalación, ejecuta `verificar.bat` o `.venv/bin/python pruebas.py`.
