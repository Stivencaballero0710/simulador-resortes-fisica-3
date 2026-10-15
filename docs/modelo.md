# Modelo físico

[Volver al inicio](../README.md)

El simulador representa una masa suspendida de un resorte vertical. La posición $x$ se mide desde el equilibrio, hacia abajo, y la masa $m$ incluye el accesorio. Por eso el peso ya está compensado en la ecuación de movimiento.

## Ecuaciones

La rigidez del resorte y su elongación estática son:

$$
k = \frac{Gd^4}{8D^3n}, \qquad \Delta\ell = \frac{mg}{k}
$$

$G$ es el módulo de corte, $d$ el diámetro del alambre, $D$ el diámetro medio y $n$ las espiras activas. Para el inoxidable se calcula $G=E/[2(1+\nu)]$.

El movimiento cumple:

$$
m\ddot{x}+b\dot{x}+kx=0, \qquad x(0)=A_0, \quad \dot{x}(0)=0
$$

Se integra con `solve_ivp` y DOP853: tolerancias `rtol=1e-10`, `atol=1e-12` y paso interno máximo de 0,03 s. Los resultados se guardan cada 0,005 s.

Se usan $\omega_0=\sqrt{k/m}$, $\gamma=b/(2m)$ y $\zeta=b/(2\sqrt{km})$. Si $\zeta<1$, entonces $\omega_d=\sqrt{\omega_0^2-\gamma^2}$ y el período es $T=2\pi/\omega_d$. Para comprobar el cálculo se utiliza:

$$
x(t)=A_0e^{-\gamma t}\left[\cos(\omega_dt)+\frac{\gamma}{\omega_d}\sin(\omega_dt)\right]
$$

En el caso crítico, $x(t)=A_0(1+\omega_0t)e^{-\omega_0t}$. En el sobreamortiguado se combinan dos exponenciales reales. Estos dos casos no tienen período.

La aceleración se calcula como $a=-(bv+kx)/m$. La energía es $E=mv^2/2+kx^2/2$ y disminuye según $dE/dt=-bv^2$. La holgura es el recorrido disponible menos $mg/k+A_0$.

## Configuración inicial

Los valores completos y su procedencia están en [parametros.json](../parametros.json).

| Dato | Valor |
| --- | --- |
| Gravedad | 9,81 m/s² |
| Geometría | Alambre de 1 mm, diámetro medio de 25 mm y 20 espiras |
| Materiales | AC: $G=79,60$ GPa; IN: $E=185$ GPa y $\nu=0,30$; TI: $G=41,40$ GPa |
| Masas para estática | 0,150; 0,200; 0,250 y 0,350 kg |
| Masas para dinámica | Referencia: 0,250 kg; comparación: 0,300 kg |
| Accesorio | 0,080 kg, ya incluido en la masa total |
| Separación inicial y recorrido | 10 mm y 220 mm |
| Registro | 15 s a 200 Hz: 3001 muestras |
| Ruido y cuantización | Desviación de 0,05 mm y paso de 0,1 mm |
| Semilla | 2250191 |

El amortiguamiento representa al sistema completo. El valor de $G$ del titanio corresponde al antecedente B120VCA citado en el informe.

## Cómo se obtienen las tablas

En estática se generan tres elongaciones por masa y material: 36 valores. A la elongación ideal se añade ruido gaussiano y luego se aplica la cuantización. Con las medias se ajusta $mg=k\Delta\ell+c$, conservando el intercepto.

En dinámica se generan tres registros de referencia (REF) y tres de comparación (VAL) por material: 18 series. Se fija $b=2\zeta\sqrt{km_{ref}}$ con $\zeta=0,020$ y se conserva ese $b$ al pasar a 0,300 kg. Solo la posición destinada al ajuste lleva ruido; velocidad, aceleración y energía se obtienen de la trayectoria limpia.

`curve_fit` estima amplitud, amortiguamiento, frecuencia y desplazamiento constante. La frecuencia inicial se obtiene con una FFT. La predicción para 0,300 kg utiliza la rigidez del ajuste estático y el promedio de $b$ de REF. Los registros VAL sirven para comparar esa predicción.

El generador es NumPy PCG64. Primero consume el ruido de estática en orden AC, IN, TI, masas ascendentes y repeticiones. Después genera la dinámica por material, REF antes de VAL y repeticiones ascendentes. La secuencia, las versiones y la configuración quedan guardadas en cada ejecución.

## Amortiguamiento y sensibilidad

El barrido usa TI con 0,300 kg y $\zeta=0;\ 0,02;\ 0,10;\ 0,30;\ 1;\ 1,5$. En cada caso se recalcula $b$.

Para $0<\zeta<1$, el tiempo al 5 % es la cota de la envolvente:

$$
t_5=\frac{\ln[1/(0,05\sqrt{1-\zeta^2})]}{\zeta\omega_0}
$$

Con $\zeta=0$ no hay estabilización permanente. Para $\zeta\geq1$ se busca el cruce de $x/A_0=0,05$ con `brentq`. El programa calcula la aceleración máxima de la trayectoria; más amortiguamiento no garantiza un máximo menor.

La sensibilidad combina $G\pm5\%$, $d\pm1\%$ y $D\pm1\%$, con $n$ fijo. Son variaciones supuestas, no intervalos de confianza.

## Alcance y relación con el informe

El modelo omite la masa del resorte, el rozamiento seco, la plasticidad, los cambios de temperatura y las resonancias internas. La holgura describe el espacio disponible; no certifica la resistencia del componente. Las dispersiones de los ajustes describen las señales sintéticas, no incertidumbres de sensores.

Los materiales y sus referencias proceden del [informe](../informe/Proyecto_Final_Fase_2.pdf). No se contó con su script original: mantener la semilla no asegura repetir las mismas series si cambia el orden de generación o el ajuste. [referencia_informe.json](../referencia_informe.json) conserva sus cifras solo para compararlas; no interviene en la generación de datos.
