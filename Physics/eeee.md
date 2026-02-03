# Black Hole Simulator / Simulador de Agujeros Negros

## Simple Description / Descripción Simple

A simple program that shows how light bends around a black hole. Watch photons (particles of light) travel through space and see how gravity affects them.

Un programa simple que muestra cómo la luz se curva alrededor de un agujero negro. Observa cómo los fotones (partículas de luz) viajan por el espacio y cómo la gravedad los afecta.

---

## Descripción del Proyecto

Visualizador en tiempo real de órbitas de fotones alrededor de un agujero negro de Schwarzschild con interfaz gráfica PyQt5. Incluye cálculo exacto de curvas isoradiales del disco de acreción usando funciones elípticas de Jacobi y renderizado de imágenes de primera y segunda orden.

## Características Principales

### 🎯 Interfaz Gráfica Interactiva (PyQt5)
- Controles deslizables y spin boxes para todos los parámetros
- Visualización en dos paneles simultáneos
- Botones para ejecutar/detener simulaciones
- Configuraciones predefinidas para casos críticos

### 🌌 Física Exacta Implementada

#### Métrica de Schwarzschild
$$ds^2 = -(1-2M/r) \, dt^2 + (1-2M/r)^{-1}dr^2 + r^2d\Omega^2$$

#### Geodésicas Nulas
- Integración numérica mediante `scipy.integrate.odeint`
- Ecuaciones de movimiento: $\ddot{r} = \frac{L^2(r - 3M)}{r^4}$
- Conservación del momento angular: $L = r^2\dot{\phi}$

#### Curvas Isoradiales del Disco de Acreción
Implementación usando **funciones elípticas de Jacobi** para calcular las curvas exactas del disco vistas desde diferentes ángulos:

1. **Parámetro Q**: $Q(P,M) = \sqrt{(P-2M)(P+6M)}$
2. **Módulo elíptico**: $k^2 = \frac{Q-P+6M}{2Q}$
3. **Funciones**: `ellipj`, `ellipkinc`, `ellipk` de scipy.special

**Primera orden** (fotones directos):
$$u(\alpha) = -A_1 + A_2 \, \text{sn}^2\left(\frac{g}{2}\sqrt{\frac{Q}{P}} + F(\zeta_\infty, k)\right)$$

**Segunda orden** (fotones que dan la vuelta):
$$u(\alpha) = -A_1 + A_2 \, \text{sn}^2\left(\frac{g-2\pi}{2}\sqrt{\frac{Q}{P}} + 2K(k) - F(\zeta_\infty, k)\right)$$

Donde:
- $u = 1/r$ (inverso del radio)
- $\text{sn}$ = función elíptica seno de Jacobi
- $F$ = integral elíptica incompleta de primera especie
- $K$ = integral elíptica completa de primera especie

### 📊 Panel Izquierdo: Geodésicas

Muestra las trayectorias de los rayos de luz:
- **Líneas azules**: Fotones que impactan el disco de acreción
- **Líneas grises**: Fotones capturados o que escapan al infinito
- **Disco rojo**: Línea del disco de acreción inclinado
- **Círculo negro**: Horizonte de eventos (r = 2M)

### 🌠 Panel Derecho: Vista del Observador

Simula cómo un observador distante vería el disco:
- **Gradiente de colores cálidos**: Temperatura del disco (más caliente = más amarillo)
- **Curvas de primera orden**: Imágenes directas del disco
- **Curvas de segunda orden**: Imágenes secundarias (luz que rodea el agujero)
- **Círculo negro central**: Sombra del agujero negro
- **Círculo semitransparente**: Esfera fotónica (r = 3M)
- **Círculo blanco exterior**: Radio crítico ($r = 3\sqrt{3}M$)

## Parámetros Ajustables

### Masa del Agujero Negro
- Rango: 1 - 10 masas solares (unidades geométricas)
- Afecta el tamaño del horizonte y la curvatura

### Parámetros de Impacto
- **Mínimo/Máximo**: Define el rango de trayectorias
- $b_{\text{crítico}} = 3\sqrt{3}M \approx 5.196M$

### Número de Rayos de Luz
- Rango: 1 - 50 rayos
- Más rayos = mejor cobertura del espacio de parámetros

### Posición Inicial
- Distancia desde la cual parten los fotones
- Típicamente -40 (desde la izquierda)

### Ángulo de Observación (Theta)
- Rango: 0° - 90°
- 0°: Vista ecuatorial (disco de canto)
- 90°: Vista polar (disco de frente)
- Afecta dramáticamente la forma aparente del disco

## Uso del Programa

### Instalación de Dependencias

```powershell
pip install PyQt5 scipy seaborn matplotlib numpy
```

## How to use it / Cómo usarlo

### Install requirements / Instalar requisitos
```powershell
pip install PyQt5 scipy seaborn matplotlib numpy
```

### Run the program / Ejecutar el programa
```powershell
python Physics/blackhole.py
```

### Controls / Controles

The interface is super simple with only 2 settings / La interfaz es muy simple con solo 2 ajustes:

1. **Number of Photons / Número de Fotones** (1-50)
   - How many light particles to show / Cuántas partículas de luz mostrar
   - More photons = more interesting patterns / Más fotones = patrones más interesantes
   - Try 10-20 for best results / Prueba con 10-20 para mejores resultados

2. **View Angle / Ángulo de Vista** (0-85 degrees / grados)
   - How you look at the black hole / Cómo miras el agujero negro
   - 0° = edge view (disk looks flat) / vista de borde (el disco se ve plano)
   - 85° = top view (disk looks round) / vista superior (el disco se ve redondo)
   - Try different angles to see how it changes! / ¡Prueba diferentes ángulos para ver cómo cambia!

### What you see / Lo que ves

The program shows two panels / El programa muestra dos paneles:

**Left panel / Panel izquierdo**: Light paths around the black hole / Trayectorias de luz alrededor del agujero negro
- Blue lines = light that hits the accretion disk / Líneas azules = luz que impacta el disco de acreción
- Gray lines = light that gets captured or escapes / Líneas grises = luz capturada o que escapa
- Red disk = the accretion disk (matter orbiting the black hole) / Disco rojo = el disco de acreción (materia orbitando el agujero negro)
- Black circle = the event horizon (point of no return) / Círculo negro = el horizonte de eventos (punto sin retorno)

**Right panel / Panel derecho**: What an observer would see / Lo que un observador vería
- Colorful rings = the hot accretion disk / Anillos coloridos = el disco de acreción caliente
- Dark center = the black hole's shadow / Centro oscuro = la sombra del agujero negro
- Notice how light bends around it! / ¡Nota cómo la luz se curva a su alrededor!

## What's happening? / ¿Qué está pasando?

### The Physics (simple version) / La Física (versión simple)

Black holes have super strong gravity that bends space itself. When light passes near a black hole:

Los agujeros negros tienen gravedad súper fuerte que curva el espacio mismo. Cuando la luz pasa cerca de un agujero negro:

- **Too close / Muy cerca**: Gets sucked in forever / Es absorbida para siempre
- **Just right / En el punto justo**: Orbits around a few times before falling in or escaping / Orbita algunas veces antes de caer o escapar
- **Far enough / Suficientemente lejos**: Bends but escapes to space / Se curva pero escapa al espacio

The critical distance is about 5.2 times the black hole's radius. Light at this distance does crazy spirals!

La distancia crítica es aproximadamente 5.2 veces el radio del agujero negro. ¡La luz a esta distancia hace espirales locas!

### Cool things to try / Cosas geniales para probar

1. **Few photons (5-10) at 20°** / Pocos fotones (5-10) a 20°: See individual paths clearly / Ver trayectorias individuales claramente
2. **Many photons (30-40) at 45°** / Muchos fotones (30-40) a 45°: Beautiful patterns emerge / Emergen patrones hermosos
3. **Any number at 80°** / Cualquier número a 80°: Top-down view, very symmetrical / Vista superior, muy simétrica
4. **15 photons at 10°** / 15 fotones a 10°: Edge view, very dramatic / Vista de borde, muy dramática

---

## Física Avanzada: Funciones Elípticas / Advanced Physics: Elliptic Functions

Las funciones elípticas son fundamentales para resolver exactamente las trayectorias de luz en el espacio curvo de Schwarzschild.

Elliptic functions are fundamental to solve exactly the light trajectories in Schwarzschild's curved spacetime.

### ¿Por qué funciones elípticas? / Why elliptic functions?

Para una órbita fotónica, la ecuación de trayectoria es:

For a photon orbit, the trajectory equation is:

$$\frac{du}{d\phi} = \pm\sqrt{1 - u^2b^2(1 - 2Mu)}$$

Esta ecuación diferencial no tiene solución en funciones elementales, pero **sí en funciones elípticas de Jacobi**.

This differential equation has no solution in elementary functions, but **it does have solutions in Jacobi elliptic functions**.

### Métrica de Schwarzschild / Schwarzschild Metric

$$ds^2 = -(1-2M/r) \, dt^2 + (1-2M/r)^{-1}dr^2 + r^2d\Omega^2$$

Donde / Where:
- $M$ = masa del agujero negro / mass of the black hole (unidades geométricas / geometric units $G=c=1$)
- $r$ = coordenada radial / radial coordinate (Schwarzschild coordinate)
- $r_s = 2M$ = radio de Schwarzschild / Schwarzschild radius (horizonte de eventos / event horizon)
- $d\Omega^2 = d\theta^2 + \sin^2\theta \, d\phi^2$ = parte angular / angular part

### Geodésicas Nulas (Trayectorias de Fotones) / Null Geodesics (Photon Trajectories)

Los fotones siguen geodésicas nulas donde $ds^2 = 0$. Las ecuaciones de movimiento son:

Photons follow null geodesics where $ds^2 = 0$. The equations of motion are:

$$\frac{dr}{dt} = \pm f\sqrt{1 - \frac{L^2}{r^2}f}$$

$$\frac{d\phi}{dt} = \frac{L}{r^2}$$

Donde / Where:
- $f = 1 - r_s/r = 1 - 2M/r$ es la función métrica / is the metric function
- $L = r^2\dot{\phi}$ es el momento angular conservado / is the conserved angular momentum
- $b = L$ es el parámetro de impacto / is the impact parameter

La aceleración radial es / The radial acceleration is:

$$\frac{d^2r}{dt^2} = \frac{L^2(r - 3M)}{r^4}$$

### Parámetro de Impacto Crítico / Critical Impact Parameter

El parámetro de impacto crítico para la captura de fotones es:

The critical impact parameter for photon capture is:

$$b_{\text{crit}} = 3\sqrt{3} \, M \approx 5.196 \, M$$

**Comportamiento del fotón / Photon behavior**:
- $b < b_{\text{crit}}$: fotón capturado / photon is captured (cae al agujero negro / falls into black hole)
- $b = b_{\text{crit}}$: fotón orbita en la esfera fotónica inestable / photon orbits on the unstable photon sphere en $r = 3M$
- $b > b_{\text{crit}}$: fotón desviado pero escapa / photon is deflected but escapes al infinito / to infinity

### Radios Clave / Key Radii

1. **Horizonte de Eventos / Event Horizon**: $r_h = 2M$
   - Radio de Schwarzschild, punto sin retorno / Schwarzschild radius, point of no return
   - Nada puede escapar desde dentro / Nothing can escape from inside

2. **Esfera Fotónica / Photon Sphere**: $r_{ph} = 3M$
   - Órbitas circulares inestables para fotones / Unstable circular orbits for photons
   - Cualquier pequeña perturbación causa captura o escape / Any small perturbation causes capture or escape

3. **Órbita Estable Más Interna (ISCO)**: $r_{ISCO} = 6M$
   - Borde interior del disco de acreción / Inner edge of accretion disk
   - Órbita estable más cercana para partículas masivas / Closest stable orbit for massive particles

4. **Radio de la Sombra / Shadow Radius**: $r_{shadow} \approx 5.2M = 3\sqrt{3}M$
   - Tamaño aparente del agujero negro para un observador distante / Apparent size of black hole to distant observer
   - Lo que ves en el panel derecho / What you see in the right panel

### Curvas Isoradiales del Disco / Accretion Disk Isoradial Curves

El programa calcula curvas isoradiales exactas usando **funciones elípticas de Jacobi**:

The program calculates exact isoradial curves using **Jacobi elliptic functions**:

$$u(\alpha) = -A_1 + A_2 \, \text{sn}^2\left(\frac{g}{2}\sqrt{\frac{Q}{P}} + F(\zeta_\infty, k)\right)$$

Donde / Where:
- $u = 1/r$ (radio inverso / inverse radius)
- $\text{sn}$ = función elíptica seno de Jacobi / Jacobi elliptic sine function
- $F$ = integral elíptica incompleta de primera especie / incomplete elliptic integral of the first kind
- $Q(P,M) = \sqrt{(P-2M)(P+6M)}$
- $k^2 = \frac{Q-P+6M}{2Q}$ (módulo elíptico / elliptic modulus)
- $A_1 = \frac{Q-P+2M}{4MP}$, $A_2 = \frac{Q-P+6M}{4MP}$

### Imágenes Múltiples / Multiple Images

El programa renderiza / The program renders:
- **Imágenes primarias / Primary images** (n=0): Trayectorias de luz directas / Direct light paths
- **Imágenes secundarias / Secondary images** (n=1): Luz que rodea el agujero negro una vez / Light that loops around the black hole once

Para imágenes secundarias / For secondary images:

$$u(\alpha) = -A_1 + A_2 \, \text{sn}^2\left(\frac{g-2\pi}{2}\sqrt{\frac{Q}{P}} + 2K(k) - F(\zeta_\infty, k)\right)$$

Donde $K(k)$ es la integral elíptica completa de primera especie / Where $K(k)$ is the complete elliptic integral of the first kind.

### Integración Numérica / Numerical Integration

Las ecuaciones geodésicas se integran usando / The geodesic equations are integrated using:
- **Método / Method**: Runge-Kutta de 4to orden via `scipy.integrate.odeint`
- **Paso de tiempo / Step size**: $\Delta t = 0.01$ (tiempo coordenado / coordinate time)
- **Tiempo total / Total time**: $t_{\text{max}} = 100M$

### Sistemas de Coordenadas / Coordinate Systems

- **Coordenadas de Schwarzschild** $(t, r, \theta, \phi)$: Usadas para cálculos / Used for calculations
- **Proyección Cartesiana** $(x, y)$: Usadas para visualización / Used for display
  - $x = r\cos\phi$
  - $y = r\sin\phi$

### Efectos del Ángulo de Observación / Viewing Angle Effects

El ángulo de observación $\theta_{\text{obs}}$ afecta la apariencia / The viewing angle affects the appearance:

- **$\theta = 0°$** (vista ecuatorial / edge-on): Disco es una línea delgada / Disk is a thin line, máxima asimetría Doppler / maximum Doppler asymmetry
- **$\theta = 45°$** (intermedio / intermediate): Estructura 3D compleja visible / Complex 3D structure visible
- **$\theta = 90°$** (vista polar / face-on): Simetría casi circular / Nearly circular symmetry, mínimo efecto Doppler / minimal Doppler effect

La transformación de ángulo es / The angle transformation is:

$$\gamma(\alpha, \theta_0) = \arccos\left(\frac{\cos\alpha}{\sqrt{\cos^2\alpha + \cot^2\theta_0}}\right)$$

---

## Referencias Científicas / Scientific References

1. **Chandrasekhar, S.** (1983). *The Mathematical Theory of Black Holes*. Oxford University Press.
   - Tratamiento completo de la geometría de Schwarzschild / Comprehensive treatment of Schwarzschild geometry

2. **Luminet, J.P.** (1979). "Image of a spherical black hole with thin accretion disk". *Astronomy and Astrophysics*, 75, 228-235.
   - Primeras simulaciones de apariencia de discos de acreción / First simulations of black hole accretion disk appearance

3. **Event Horizon Telescope Collaboration** (2019). "First M87 Event Horizon Telescope Results". *The Astrophysical Journal Letters*, 875:L1.
   - Primera imagen directa de la sombra de un agujero negro / First direct image of a black hole shadow

4. **Gralla, S.E. & Lupsasca, A.** (2020). "Lensing by Kerr Black Holes". *Physical Review D*, 101, 044031.
   - Teoría moderna de lentes gravitacionales por agujeros negros / Modern theory of gravitational lensing by black holes

---

## Detalles de Implementación / Implementation Details

### Software Stack
- **Python**: 3.14.3
- **GUI**: PyQt5 (interfaz multiplataforma / cross-platform interface)
- **Numerics**: NumPy, SciPy (odeint, fsolve, funciones elípticas / elliptic functions)
- **Visualización / Visualization**: Matplotlib (gráficas 2D / 2D plotting)
- **Colores / Colors**: Seaborn (paletas térmicas / thermal color palettes)

### Performance
- Cálculo de curvas isoradiales / Isoradial curve calculation: 10-30 segundos / seconds (depende del ángulo / depends on angle)
- Trazado de rayos de luz / Light ray tracing: 0.1-2 segundos / seconds (depende del número de fotones / depends on number of photons)
- Usa arrays NumPy optimizados para vectorización / Uses optimized NumPy arrays for vectorization

### Limitaciones / Limitations
- Solo agujeros negros no rotantes / Non-rotating black holes only (Schwarzschild, not Kerr)
- Aproximación de disco delgado / Thin disk approximation (sin grosor / no thickness)
- Sin impulso Doppler relativista en colores / No relativistic Doppler boosting in colors
- Sin visualización de corrimiento al rojo gravitacional / No gravitational redshift visualization
- Imágenes estáticas / Static images (sin evolución temporal / no time evolution)

---

## Troubleshooting

**El programa no inicia / Program won't start**: Asegúrate de instalar todos los requisitos / Make sure you installed all requirements
```powershell
pip install PyQt5 scipy seaborn matplotlib numpy
```

**Muy lento / Too slow**: 
- Reduce el número de fotones a 10 o menos / Reduce number of photons to 10 or less
- ¡El programa necesita calcular física compleja, ten paciencia! / The program needs to calculate complex physics, be patient!

**No aparece nada / Nothing shows up**: Haz clic en "RUN SIMULATION" después de ajustar tus parámetros / Click "RUN SIMULATION" button after setting your parameters

**Las ecuaciones no se renderizan / Equations not rendering**: Si estás viendo en GitHub, las ecuaciones deben mostrarse automáticamente. Para visualización local, usa un visor de Markdown que soporte LaTeX/MathJax. / If viewing on GitHub, equations should display automatically. For local viewing, use a Markdown viewer that supports LaTeX/MathJax.

---

## Créditos / Credits

**Física / Physics**: Albert Einstein (Relatividad General / General Relativity, 1915)  
**Matemáticas / Mathematics**: Carl Gustav Jacobi (Funciones Elípticas / Elliptic Functions, 1829)  
**Implementación / Implementation**: Basado en / Based on McGill Physics Hackathon 2022  
**Software**: Python 3.14.3 con / with PyQt5, matplotlib, scipy, numpy, seaborn

**¡Diviértete explorando agujeros negros! / Have fun exploring black holes!** 🌌

#### 4. **Optimización de Velocidad**
- `DT = 0.01`: Paso de integración 10× más grande
- `substeps = 30`: 30 pasos de integración por frame
- `interval = 10ms`: Actualización rápida de la animación
- Resultado: velocidad ~60× más rápida que la versión inicial

#### 5. **Estética Mejorada**
- Fondo negro espacial (`#1a1a1a`)
- Texto blanco para máximo contraste
- Mapas de colores científicos (`plasma`, `hot`)
- Bordes blancos en los fotones para visibilidad
- Grid sutil para referencia espacial

## Uso del Programa

### Requisitos
```bash
pip install numpy matplotlib
```

### Ejecución
```bash
python Physics/blackhole.py
```

1. El programa preguntará: "¿Cuántos fotones quieres simular?"
2. Ingresa un número (recomendado: 5-20 fotones)
3. Se abrirá una ventana con dos paneles mostrando la simulación en tiempo real

### Controles
- Cerrar la ventana para terminar la simulación
- Los fotones desaparecen automáticamente al caer en el horizonte de eventos

## Resultados Esperados

Verás diferentes comportamientos según el parámetro de impacto de cada fotón:

1. **Captura directa**: Fotones que caen rápidamente en espiral hacia el agujero negro
2. **Órbitas críticas**: Fotones que dan varias vueltas alrededor de la esfera fotónica antes de caer
3. **Dispersión**: Fotones que se curvan pero escapan hacia el infinito

El panel derecho muestra cómo se vería el disco de acreción, similar a las famosas imágenes del Event Horizon Telescope.

## Física Avanzada (Referencia)

El código de referencia proporcionado implementa características adicionales más avanzadas:

- Funciones elípticas de Jacobi para curvas isoradiales exactas
- Vista del disco desde diferentes ángulos de inclinación
- Interfaz interactiva con PyQt5
- Cálculo de segunda orden para imágenes secundarias del disco
- Detección de impactos con el disco de acreción

Estas características podrían agregarse en futuras versiones para una simulación aún más realista.

## Referencias

- Chandrasekhar, S. (1983). *The Mathematical Theory of Black Holes*
- Misner, Thorne, Wheeler (1973). *Gravitation*
- Event Horizon Telescope Collaboration (2019). First M87 Black Hole Image

## Autor

Simulación creada con fines educativos para visualizar la relatividad general en acción.
