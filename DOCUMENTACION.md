# Documentación del simulador de resortes — Física 3

El programa implementa la suspensión vertical idealizada del proyecto de Física 3, grupo F1B de la UIS. Compara el acero Oteva 90 (AC), el inoxidable AISI 302 (IN) y el titanio B120VCA (TI), genera registros sintéticos, ajusta sus parámetros y permite explorar el movimiento en una interfaz local.

**Versión 1.0.0.** Implementación elaborada a partir del informe. El programa se ejecutó en Python 3.12.14: pasaron las nueve pruebas y las 18 comparaciones analíticas; la exportación completa terminó correctamente. Se incluye una ejecución en resultados_ejemplo. Los controles se comprobaron sin pantalla y la vista generada se inspeccionó visualmente. Consulte VALIDACION.md para conocer el alcance y las versiones.

# Archivos incluidos

| Archivo | Función |
| --- | --- |
| simulador.py | Programa principal: modelo físico, generación, ajustes, ocho tablas, cinco gráficas, interfaz, GIF y verificación. |
| parametros.json | Condiciones del escenario, materiales, procedencia y opciones numéricas. |
| referencia_informe.json | Valores publicados en las tablas 2 a 7, usados únicamente para comparar resultados. |
| requirements.txt | Dependencias con intervalos de versiones. |
| pruebas.py | Nueve pruebas de física, cálculo y entradas. |
| LEEME.txt | Inicio rápido y comandos de ejecución. |
| DOCUMENTACION.md | Manual en texto editable y legible en GitHub. |
| DOCUMENTACION.html | Manual completo para abrir e imprimir en el navegador. |
| instalar.bat | Crea el entorno local de Python e instala las bibliotecas. |
| abrir_simulador.bat | Abre la interfaz en Windows. |
| generar_resultados.bat | Genera tablas, series, gráficas y animación. |
| verificar.bat | Ejecuta las pruebas en Windows. |

# Instalación y ejecución

En Windows, extraiga toda la carpeta del ZIP e instale Python 3.11 o posterior con Tcl/Tk. Ejecute **instalar.bat** y espere a que termine. Después ejecute **verificar.bat**. La instalación necesita internet; el simulador funciona localmente.

Para usar los controles, abra **abrir_simulador.bat**. Para generar todos los archivos, abra **generar_resultados.bat**. La ventana de consola muestra los errores y permanece abierta al terminar.

También puede usar una terminal desde la carpeta descomprimida:

```text
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe pruebas.py
.venv\Scripts\python.exe simulador.py
.venv\Scripts\python.exe simulador.py --interactivo
```

En Linux y macOS sustituya la ruta del intérprete por .venv/bin/python. La interfaz requiere un entorno gráfico y un backend de Matplotlib compatible, como Tk. El modo de exportación utiliza Agg y no necesita abrir ventanas.

| Opción | Comportamiento |
| --- | --- |
| Sin opciones | Genera el análisis completo y el GIF. |
| --interactivo | Abre solamente la interfaz. |
| --verificar | Ejecuta solamente las comparaciones con las soluciones analíticas. |
| --sin-gif | Genera el análisis completo omitiendo el GIF. |
| --config ruta.json | Carga una copia modificada de los parámetros. |
| --semilla 2250191 | Fija la semilla de esa ejecución. |
| --salida resultados/ensayo_01 | Usa una carpeta nueva con ese nombre. |

No se sobrescriben carpetas existentes. Sin --salida, se crea una carpeta con fecha, hora y microsegundos dentro de resultados.

# Modelo físico

La coordenada x se mide desde el equilibrio, con sentido positivo hacia abajo. La masa m incluye la carga y el accesorio. El peso está compensado en esa posición, por lo que no se agrega otra fuerza constante a la ecuación dinámica.

La rigidez se calcula mediante k = G d⁴/(8 D³ n). Para IN se obtiene G = E/[2(1 + ν)]. La elongación estática ideal es Δℓ = mg/k. Durante la caracterización se ajusta mg = k Δℓ + c, conservando el intercepto c.

El sistema que integra el programa es dx/dt = v y dv/dt = −(bv + kx)/m. Las condiciones iniciales son x(0) = A₀ y v(0) = 0. Se utiliza solve_ivp con DOP853, rtol = 10⁻¹⁰, atol = 10⁻¹² y paso interno máximo de 0,03 s. El paso de salida de 0,005 s es distinto del paso adaptativo interno.

Se definen ω₀ = √(k/m), γ = b/(2m) y ζ = b/(2√(km)). Cuando ζ < 1, ωd = √(ω₀² − γ²) y T = 2π/ωd. El programa no asigna un período a los regímenes crítico y sobreamortiguado.

La solución subamortiguada utilizada para comprobar la integración es x(t) = A₀ exp(−γt)[cos(ωd t) + (γ/ωd) sen(ωd t)]. El término seno permite cumplir v(0) = 0. Para ζ = 1 se utiliza x(t) = A₀(1 + ω₀t)exp(−ω₀t). Para ζ > 1 se emplea la combinación de las dos exponenciales reales, calculando la raíz lenta sin cancelación numérica.

La aceleración se obtiene de a = −(bv + kx)/m. La energía respecto al equilibrio es E = mv²/2 + kx²/2 y cumple dE/dt = −bv². No se deriva dos veces la señal ruidosa para calcular la aceleración.

La holgura es H = Límite − (mg/k + A₀). Esta comparación utiliza el recorrido geométrico asignado, no un límite elástico certificado del resorte.

# Parámetros del informe

| Magnitud | Valor |
| --- | --- |
| Gravedad | 9,81 m/s² |
| Diámetro del alambre d | 0,001 m |
| Diámetro medio D | 0,025 m |
| Espiras activas n | 20 |
| G de AC | 79,60 GPa |
| E de IN y ν supuesto | 185 GPa y 0,30 |
| G de IN calculado | 71,1538461538 GPa |
| G de TI | 41,40 GPa |
| Masas estáticas | 0,150; 0,200; 0,250 y 0,350 kg |
| Masa de referencia | 0,250 kg |
| Masa reservada | 0,300 kg |
| Masa del accesorio | 0,080 kg, incluida en la masa total |
| Separación inicial A₀ | 0,0100 m |
| Recorrido asignado | 0,220 m |
| ζ de referencia | 0,020 |
| Duración y muestreo | 15 s y 200 Hz: 3001 muestras |
| Desviación del ruido | 0,00005 m = 0,05 mm |
| Paso de cuantización | 0,0001 m = 0,1 mm |
| Semilla | 2250191 |
| Barrido independiente de ζ | 0; 0,02; 0,10; 0,30; 1; 1,5 |

Los parámetros están expresados en unidades del SI en parametros.json. El módulo del titanio procede del antecedente citado en el informe y no representa a todas las aleaciones de titanio. El amortiguamiento es una propiedad equivalente del sistema; no se trata como una constante propia de cada metal.

# Generación y tratamiento de datos

Primero se crea un generador NumPy Generator(PCG64) con la semilla indicada. Se generan todas las elongaciones estáticas en el orden AC, IN y TI; dentro de cada material se recorren las masas ascendentes y las repeticiones. Después se generan los registros dinámicos por material, con REF antes de VAL y las repeticiones ascendentes. Este orden queda registrado.

Para la estática se calculan tres elongaciones por masa y material: 36 valores distribuidos en 12 condiciones. A mg/k se añade ruido gaussiano y luego se redondea al paso de cuantización. Se calcula la media de las tres repeticiones y se ajusta fuerza frente a elongación mediante mínimos cuadrados. Se exportan k, c, R² y los errores estándar de la regresión.

Para la dinámica se fija b = 2ζ√(k mref) con mref = 0,250 kg y ζ = 0,020. Ese b se conserva cuando la masa pasa a 0,300 kg. Por tanto, la razón de amortiguamiento de la masa reservada cambia a aproximadamente 0,01826.

Se generan tres registros REF y tres VAL por material: 18 series. Cada archivo conserva t, x limpia, x sintética, v, a y E. El ruido afecta la señal de posición destinada al ajuste; la velocidad, la aceleración y la energía corresponden a la trayectoria limpia.

El ajuste dinámico estima A, γ, ωd y un desplazamiento constante mediante curve_fit. La frecuencia inicial se estima con la transformada rápida de Fourier de la señal. Después se obtienen T = 2π/ωd y b = 2mγ. Se exportan el error cuadrático medio y los errores estándar locales del ajuste.

La predicción de 0,300 kg utiliza exclusivamente la rigidez del ajuste estático y el promedio de b de los tres registros REF. Los ajustes de los registros VAL se usan para contrastar sus períodos; no se utilizan para corregir la predicción. La diferencia relativa se calcula con el promedio del período reservado como denominador.

Las dispersiones y covarianzas describen datos sintéticos y el ajuste empleado. No son incertidumbres instrumentales ni demuestran la variabilidad de un resorte físico. Los resultados conservan más cifras en los CSV que en el informe.

# Amortiguamiento y sensibilidad

En el barrido se usa TI con m = 0,300 kg. Para cada ζ se calcula nuevamente b = 2ζ√(km). Este procedimiento es independiente de la comparación REF/VAL.

Cuando 0 < ζ < 1, el tiempo al 5 % es una cota basada en la amplitud exponencial: t₅ = ln[1/(0,05√(1 − ζ²))]/(ζω₀). Esta cota garantiza que la amplitud posterior permanezca por debajo del umbral. Para ζ = 0 el tiempo no es finito. Para ζ ≥ 1 se busca con brentq el cruce monótono de x/A₀ = 0,05.

La aceleración inicial es −kA₀/m. No depende de b porque la velocidad inicial es cero. El programa calcula el máximo absoluto de aceleración de la trayectoria; no supone que aumentar el amortiguamiento reduzca ese máximo.

La sensibilidad combina G ±5 %, d ±1 % y D ±1 %, con n fijo. El mínimo k se obtiene usando G y d mínimos y D máximo; para el máximo se invierten esos extremos. Los intervalos son supuestos deterministas, no intervalos de confianza.

# Resultados nominales de referencia

Estos valores se calculan directamente a partir de los parámetros, sin ruido ni ajustes. Sirven para revisar las unidades y la configuración.

| Caso | k (N/m) | T con masa reservada (s) | a máxima (m/s²) | Holgura (cm) |
| --- | ---: | ---: | ---: | ---: |
| AC | 31,840000 | 0,609995 | 1,061333 | 11,756910 |
| IN | 28,461538 | 0,645184 | 0,948718 | 10,659730 |
| TI | 16,560000 | 0,845829 | 0,552000 | 3,228261 |

Con estas condiciones, TI conserva holgura positiva y presenta la menor aceleración. La conclusión corresponde a configuraciones idealizadas; no acredita la resistencia, el costo ni la sostenibilidad de componentes reales.

# Archivos generados

| Salida | Contenido |
| --- | --- |
| tabla_1.csv | Parámetros de la ejecución y procedencia de los módulos. |
| tabla_2.csv | Elongaciones estáticas, medias y dispersión sintética. |
| tabla_3.csv | Parámetros ajustados de los 18 registros dinámicos. |
| tabla_4.csv | Predicción reservada, promedio ajustado y diferencia relativa. |
| tabla_5.csv | Respuesta nominal, aceleración máxima y holgura. |
| tabla_6.csv | Barrido de amortiguamiento y tiempo al 5 %. |
| tabla_7.csv | Extremos de rigidez y holgura con k mínimo. |
| tabla_8.csv | Propuesta de reutilización del informe; no es un inventario medido. |
| ajustes_estaticos.csv | Pendientes, interceptos, covarianzas resumidas y R². |
| series/AC-REF-1.csv y análogos | Las 18 series con trayectoria limpia y señal sintética. |
| series/AC-PRED.csv y análogos | Las tres predicciones calculadas desde REF. |
| series/TI-zeta-*.csv | Las seis trayectorias del barrido independiente. |
| graficas/figura_1 a figura_5 | PNG de 240 ppp y PDF vectoriales. |
| graficas/figura_A1_interfaz.png | Captura generada de la interfaz de esta implementación. |
| graficas/animacion_TI.gif | Animación TI, masa reservada y ζ = 0,30. |
| verificacion_numerica.csv | Errores frente a la solución analítica y pruebas de energía. |
| comparacion_con_informe.csv | Diferencias frente a las cifras redondeadas de las tablas originales. |
| parametros_ejecucion.json | Configuración efectiva, versiones, semilla, orden de generación y huella del código. |
| RESUMEN.txt | Selección nominal y resumen de las comprobaciones de esa ejecución. |
| ESTADO.txt | Estado final de la verificación numérica. |

Cada CSV tiene una copia terminada en _es.csv, con separador punto y coma y coma decimal. El CSV estándar utiliza coma como separador y punto decimal. Los archivos se codifican en UTF-8 con marca de codificación para facilitar su lectura.

En la tabla 6, un período vacío significa que no aplica. Un tiempo inf significa que no se alcanza el umbral de manera permanente. No deben interpretarse como cero.

# Uso de la interfaz

Seleccione AC, IN o TI. La barra de masa controla la masa total suspendida, incluido el accesorio. La barra ζ cambia el amortiguamiento del sistema y recalcula b. La barra A₀ controla la separación inicial en milímetros.

La gráfica derecha muestra la posición respecto al equilibrio. La hélice de la izquierda es un esquema de la elongación; no representa la longitud real del alambre ni el esfuerzo mecánico. La línea discontinua señala el recorrido asignado. La etiqueta de holgura indica si la configuración cabe dentro de ese supuesto geométrico.

Pausa / continuar detiene o reanuda el movimiento. Reiniciar tiempo devuelve el marcador al inicio. Cambiar un parámetro recalcula la trayectoria y reinicia el tiempo. Los cambios de los controles no modifican parametros.json ni los resultados ya exportados.

# Organización del código

| Función | Responsabilidad |
| --- | --- |
| cargar_configuracion | Lee parámetros, valida entradas, calcula G de IN, k y b de referencia. |
| integrar / resolver | Integra el sistema con DOP853 y obtiene x, v, a y energía. |
| analitica | Calcula soluciones independientes para comprobar la integración. |
| tiempo_al_5 | Obtiene la cota exponencial o resuelve el cruce monótono. |
| senal_sintetica | Añade ruido gaussiano y aplica cuantización. |
| ajustar_estatica | Estima k e intercepto con mínimos cuadrados. |
| ajustar_dinamica | Estima parámetros de la oscilación a partir de la señal sintética. |
| generar | Ejecuta el protocolo, organiza tablas y guarda las series. |
| comparar_informe | Contrasta las salidas con cifras históricas sin usarlas en el cálculo. |
| verificar | Comprueba 18 combinaciones de material y amortiguamiento. |
| graficar | Produce las cinco figuras principales. |
| crear_interfaz | Construye los controles y la animación local. |
| exportar_gif | Guarda una animación con PillowWriter. |
| main | Procesa opciones, elige el modo y registra la ejecución. |

NumPy administra vectores, ruido, FFT y ajustes estáticos. SciPy resuelve la ecuación diferencial, ajusta las señales y busca raíces. Pandas organiza y exporta los datos. Matplotlib dibuja las gráficas y los controles. Pillow permite guardar el GIF.

# Verificación y reproducibilidad

pruebas.py contiene nueve pruebas. Una comprueba las rigideces; otra revisa el caso sin amortiguamiento en un cuarto, medio y un período. También se comprueban los regímenes no oscilatorios, el tiempo al 5 %, las 18 comparaciones analíticas, la recuperación de parámetros desde una señal limpia, la repetibilidad del ruido, las magnitudes nominales y el rechazo de entradas inválidas.

La comparación numérica exige error máximo de posición menor de 10⁻⁸ m y de velocidad menor de 10⁻⁷ m/s. La energía sin amortiguamiento debe tener una deriva relativa menor de 10⁻⁷. Con amortiguamiento se permite solamente un aumento relativo por paso inferior a 10⁻⁸ para tolerar el redondeo numérico. Los resultados reales se escriben al ejecutar; no se sustituye el error calculado por el 8,16 × 10⁻¹² m citado en el informe.

No se dispuso del script original que generó las series del informe. Aunque se conserva la semilla 2250191, la secuencia depende también del generador, del orden de consumo del ruido, del modelo de ajuste y de las versiones. Por ello no se garantiza igualdad literal de las últimas cifras de los resultados sintéticos históricos.

referencia_informe.json contiene datos redondeados transcritos del informe y solo se abre desde comparar_informe. Los registros nuevos siempre se generan desde las ecuaciones. Las diferencias se entregan en un archivo aparte. Si se cambia la configuración, esa comparación continúa usando los valores históricos y debe interpretarse como tal.

Para repetir una ejecución conserve parametros_ejecucion.json, el código y las versiones instaladas. Puede registrar las versiones exactas con python -m pip freeze > versiones_instaladas.txt. Los intervalos de requirements.txt facilitan la instalación, pero no fijan un entorno idéntico en todas las fechas.

# Problemas frecuentes

| Situación | Qué revisar |
| --- | --- |
| Python no se reconoce | Instale Python y abra de nuevo la terminal; en Windows puede usar el lanzador py. |
| Falta numpy, scipy u otra biblioteca | Ejecute instalar.bat o instale requirements.txt con el intérprete de .venv. |
| La interfaz no se abre | Use Python de escritorio con Tcl/Tk; pruebe python -m tkinter. En servidores sin pantalla utilice la exportación. |
| FileExistsError al guardar | Use otra carpeta con --salida o deje que se cree una carpeta automática. |
| La exportación tarda | El GIF guarda muchos fotogramas; use --sin-gif si solo necesita tablas y gráficas. |
| Las últimas cifras difieren del informe | Revise comparacion_con_informe.csv y las notas de reproducibilidad. No reemplace los datos calculados por las cifras esperadas. |
| Una prueba falla | Revise el mensaje, los parámetros y las versiones antes de presentar la ejecución como verificada. |

# Alcance físico y procedencia

La simulación reemplaza la adquisición prevista con CASSY-Lab porque, según el informe, no fue posible coordinar el acceso al laboratorio con el señor Yebrail. Todos los archivos producidos se identifican como resultados de simulación.

Se omiten la masa del resorte, el rozamiento seco, la precarga, la plasticidad, la variación térmica y las resonancias internas. La comparación con una solución analítica comprueba la implementación del modelo. Una validación experimental requeriría caracterizar resortes y sujeciones reales.

Las propiedades se conservan como aparecen en el informe; no se realizó una nueva verificación de los valores de los artículos originales para esta implementación.

| Referencia del informe | Uso |
| --- | --- |
| Gu, Hou y Ye (2021). DOI 10.1177/09544062211010210 | Módulo de corte de AC y relación con la geometría. |
| Zani y Solazzi (2024). DOI 10.3390/eng5030089 | Módulo de Young de IN; ν = 0,30 sigue siendo un supuesto. |
| Li et al. (2024). DOI 10.3390/aerospace11070514 | Antecedente bibliográfico de G de B120VCA citado por el informe. |
| Yokoyama (2023). DOI 10.1007/s00419-023-02414-3 | Orientación del análisis de amortiguamiento. |
| Lee y Thompson (2001). DOI 10.1006/jsvi.2000.3169 | Límites del modelo de resorte sin masa. |

# Documentación oficial de las bibliotecas

- [SciPy solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html).

- [SciPy curve_fit](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html).

- [NumPy Generator](https://numpy.org/doc/stable/reference/random/generator.html).

- [Matplotlib animaciones](https://matplotlib.org/stable/api/animation_api.html).

- [Matplotlib controles](https://matplotlib.org/stable/api/widgets_api.html).


