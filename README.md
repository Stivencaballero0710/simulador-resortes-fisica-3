# Simulador de resortes — Física 3

Proyecto del grupo F1B de la UIS. Compara resortes de acero Oteva 90 (AC), acero inoxidable AISI 302 (IN) y titanio B120VCA (TI). Puedes cambiar la masa y el amortiguamiento para ver cómo se mueve cada uno.

Los datos son simulados a partir del modelo del informe.

## Cómo abrirlo en Windows

Necesitas Python 3.11 o posterior, con Tcl/Tk para abrir la ventana.

1. Descarga el proyecto desde **Code → Download ZIP** y descomprime la carpeta.
2. Abre `instalar.bat` y espera a que termine. Este paso necesita internet.
3. Abre `abrir_simulador.bat`.

Para guardar tablas, gráficas y una animación, usa `generar_resultados.bat`. Para comprobar la instalación, usa `verificar.bat`.

![Ventana del simulador](resultados_ejemplo/graficas/figura_A1_interfaz.png)

## Dónde encontrar cada cosa

| Archivo o carpeta | Para qué sirve |
| --- | --- |
| [simulador.py](simulador.py) | Código del programa. |
| [parametros.json](parametros.json) | Masas, materiales y demás ajustes. |
| [resultados_ejemplo/](resultados_ejemplo) | Tablas y gráficas ya generadas. |
| [informe/](informe) | Informe en PDF y Word. |
| [fuentes/](fuentes) | Documentos de partida. |

[Guía de uso y comandos](docs/guia.md) · [Modelo físico](docs/modelo.md) · [Pruebas realizadas](docs/verificacion.md)
