# Simulador de Resortes - Fisica 3

Simulador interactivo de suspension vertical para resortes helicoidales:
- Acero Oteva 90 (AC)
- Inoxidable AISI 302 (IN)
- Titanio B120VCA (TI)
- Material Personalizado

Calcula rigidez, periodo de oscilacion, frecuencia propia, elongacion de equilibrio, aceleracion maxima, holgura de seguridad y evolucion temporal del sistema. Permite exportar resultados a Excel (.xlsx) y CSV (.csv).

---

## Ejecucion en Google Colab

Crea un nuevo cuaderno en [Google Colab](https://colab.research.google.com/) y ejecuta en una celda:

```python
# 1. Clonar el repositorio
!git clone https://github.com/Stivencaballero0710/simulador-resortes-fisica-3.git
%cd simulador-resortes-fisica-3

# 2. Instalar dependencias
!pip install -r requirements.txt

# 3. Ejecutar simulador interactivo
!python simulador.py
```

---

## Ejecucion en Terminal Local (Windows / Linux / Mac)

```bash
# 1. Instalar dependencias
pip install -r requirements.txt

# 2. Ejecutar
python simulador.py
```

---

## Archivos del Repositorio

- `simulador.py`: Codigo principal del simulador con interfaz por terminal.
- `parametros.json`: Parametros fisicos y materiales de referencia.
- `requirements.txt`: Librerias necesarias (numpy, scipy, pandas, openpyxl).
