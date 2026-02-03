# Black Hole Simulator

## Simple Description

A simple program that shows how light bends around a black hole. Watch photons (particles of light) travel through space and see how gravity affects them.

---

## Project Description

Real-time visualizer of photon orbits around a Schwarzschild black hole with PyQt5 graphical interface. Includes exact calculation of isoradial curves of the accretion disk using Jacobi elliptic functions and rendering of first and second order images.

## Main Features

### 🎯 Interactive Graphical Interface (PyQt5)
- Sliders and spin boxes for all parameters
- Simultaneous two-panel visualization
- Buttons to run/stop simulations
- Predefined configurations for critical cases

### 🌌 Exact Physics Implementation

#### Schwarzschild Metric
$$ds^2 = -(1-2M/r) \, dt^2 + (1-2M/r)^{-1}dr^2 + r^2d\Omega^2$$

#### Null Geodesics
- Numerical integration using `scipy.integrate.odeint`
- Equations of motion: $\ddot{r} = \frac{L^2(r - 3M)}{r^4}$
- Conservation of angular momentum: $L = r^2\dot{\phi}$

#### Accretion Disk Isoradial Curves
Implementation using **Jacobi elliptic functions** to calculate exact disk curves viewed from different angles:

1. **Parameter Q**: $Q(P,M) = \sqrt{(P-2M)(P+6M)}$
2. **Elliptic modulus**: $k^2 = \frac{Q-P+6M}{2Q}$
3. **Functions**: `ellipj`, `ellipkinc`, `ellipk` from scipy.special

**First order** (direct photons):
$$u(\alpha) = -A_1 + A_2 \, \text{sn}^2\left(\frac{g}{2}\sqrt{\frac{Q}{P}} + F(\zeta_\infty, k)\right)$$

**Second order** (photons that loop around):
$$u(\alpha) = -A_1 + A_2 \, \text{sn}^2\left(\frac{g-2\pi}{2}\sqrt{\frac{Q}{P}} + 2K(k) - F(\zeta_\infty, k)\right)$$

Where:
- $u = 1/r$ (inverse radius)
- $\text{sn}$ = Jacobi elliptic sine function
- $F$ = incomplete elliptic integral of the first kind
- $K$ = complete elliptic integral of the first kind

### 📊 Left Panel: Geodesics

Shows the trajectories of light rays:
- **Blue lines**: Photons that hit the accretion disk
- **Gray lines**: Photons captured or escaping to infinity
- **Red disk**: Inclined accretion disk line
- **Black circle**: Event horizon (r = 2M)

### 🌠 Right Panel: Observer View

Simulates how a distant observer would see the disk:
- **Warm color gradient**: Disk temperature (hotter = more yellow)
- **First order curves**: Direct images of the disk
- **Second order curves**: Secondary images (light wrapping around the hole)
- **Central black circle**: Black hole shadow
- **Semi-transparent circle**: Photon sphere (r = 3M)
- **Outer white circle**: Critical radius ($r = 3\sqrt{3}M$)

## Adjustable Parameters

### Black Hole Mass
- Range: 1 - 10 solar masses (geometric units)
- Affects horizon size and curvature

### Impact Parameters
- **Minimum/Maximum**: Defines the range of trajectories
- $b_{\text{critical}} = 3\sqrt{3}M \approx 5.196M$

### Number of Light Rays
- Range: 1 - 50 rays
- More rays = better coverage of parameter space

### Initial Position
- Distance from which photons start
- Typically -40 (from the left)

### Observation Angle (Theta)
- Range: 0° - 90°
- 0°: Equatorial view (edge-on disk)
- 90°: Polar view (face-on disk)
- Dramatically affects the apparent shape of the disk

## Program Usage

### Installing Dependencies

```powershell
pip install PyQt5 scipy seaborn matplotlib numpy
```

## How to use it

### Install requirements
```powershell
pip install PyQt5 scipy seaborn matplotlib numpy
```

### Run the program
```powershell
python Physics/blackhole.py
```

### Controls

The interface is super simple with only 2 settings:

1. **Number of Photons** (1-50)
   - How many light particles to show
   - More photons = more interesting patterns
   - Try 10-20 for best results

2. **View Angle** (0-85 degrees)
   - How you look at the black hole
   - 0° = edge view (disk looks flat)
   - 85° = top view (disk looks round)
   - Try different angles to see how it changes!

### What you see

The program shows two panels:

**Left panel**: Light paths around the black hole
- Blue lines = light that hits the accretion disk
- Gray lines = light that gets captured or escapes
- Red disk = the accretion disk (matter orbiting the black hole)
- Black circle = the event horizon (point of no return)

**Right panel**: What an observer would see
- Colorful rings = the hot accretion disk
- Dark center = the black hole's shadow
- Notice how light bends around it!

## What's happening?

### The Physics (simple version)

Black holes have super strong gravity that bends space itself. When light passes near a black hole:

- **Too close**: Gets sucked in forever
- **Just right**: Orbits around a few times before falling in or escaping
- **Far enough**: Bends but escapes to space

The critical distance is about 5.2 times the black hole's radius. Light at this distance does crazy spirals!

### Cool things to try

1. **Few photons (5-10) at 20°**: See individual paths clearly
2. **Many photons (30-40) at 45°**: Beautiful patterns emerge
3. **Any number at 80°**: Top-down view, very symmetrical
4. **15 photons at 10°**: Edge view, very dramatic

---

## Advanced Physics: Elliptic Functions

Elliptic functions are fundamental to solve exactly the light trajectories in Schwarzschild's curved spacetime.

### Why elliptic functions?

For a photon orbit, the trajectory equation is:

$$\frac{du}{d\phi} = \pm\sqrt{1 - u^2b^2(1 - 2Mu)}$$

This differential equation has no solution in elementary functions, but **it does have solutions in Jacobi elliptic functions**.

### Schwarzschild Metric

$$ds^2 = -(1-2M/r) \, dt^2 + (1-2M/r)^{-1}dr^2 + r^2d\Omega^2$$

Where:
- $M$ = mass of the black hole (geometric units $G=c=1$)
- $r$ = radial coordinate (Schwarzschild coordinate)
- $r_s = 2M$ = Schwarzschild radius (event horizon)
- $d\Omega^2 = d\theta^2 + \sin^2\theta \, d\phi^2$ = angular part

### Null Geodesics (Photon Trajectories)

Photons follow null geodesics where $ds^2 = 0$. The equations of motion are:

$$\frac{dr}{dt} = \pm f\sqrt{1 - \frac{L^2}{r^2}f}$$

$$\frac{d\phi}{dt} = \frac{L}{r^2}$$

Where:
- $f = 1 - r_s/r = 1 - 2M/r$ is the metric function
- $L = r^2\dot{\phi}$ is the conserved angular momentum
- $b = L$ is the impact parameter

The radial acceleration is:

$$\frac{d^2r}{dt^2} = \frac{L^2(r - 3M)}{r^4}$$

### Critical Impact Parameter

The critical impact parameter for photon capture is:

$$b_{\text{crit}} = 3\sqrt{3} \, M \approx 5.196 \, M$$

**Photon behavior**:
- $b < b_{\text{crit}}$: photon is captured (falls into black hole)
- $b = b_{\text{crit}}$: photon orbits on the unstable photon sphere at $r = 3M$
- $b > b_{\text{crit}}$: photon is deflected but escapes to infinity

### Key Radii

1. **Event Horizon**: $r_h = 2M$
   - Schwarzschild radius, point of no return
   - Nothing can escape from inside

2. **Photon Sphere**: $r_{ph} = 3M$
   - Unstable circular orbits for photons
   - Any small perturbation causes capture or escape

3. **Innermost Stable Circular Orbit (ISCO)**: $r_{ISCO} = 6M$
   - Inner edge of accretion disk
   - Closest stable orbit for massive particles

4. **Shadow Radius**: $r_{shadow} \approx 5.2M = 3\sqrt{3}M$
   - Apparent size of black hole to distant observer
   - What you see in the right panel

### Accretion Disk Isoradial Curves

The program calculates exact isoradial curves using **Jacobi elliptic functions**:

$$u(\alpha) = -A_1 + A_2 \, \text{sn}^2\left(\frac{g}{2}\sqrt{\frac{Q}{P}} + F(\zeta_\infty, k)\right)$$

Where:
- $u = 1/r$ (inverse radius)
- $\text{sn}$ = Jacobi elliptic sine function
- $F$ = incomplete elliptic integral of the first kind
- $Q(P,M) = \sqrt{(P-2M)(P+6M)}$
- $k^2 = \frac{Q-P+6M}{2Q}$ (elliptic modulus)
- $A_1 = \frac{Q-P+2M}{4MP}$, $A_2 = \frac{Q-P+6M}{4MP}$

### Multiple Images

The program renders:
- **Primary images** (n=0): Direct light paths
- **Secondary images** (n=1): Light that loops around the black hole once

For secondary images:

$$u(\alpha) = -A_1 + A_2 \, \text{sn}^2\left(\frac{g-2\pi}{2}\sqrt{\frac{Q}{P}} + 2K(k) - F(\zeta_\infty, k)\right)$$

Where $K(k)$ is the complete elliptic integral of the first kind.

### Numerical Integration

The geodesic equations are integrated using:
- **Method**: 4th-order Runge-Kutta via `scipy.integrate.odeint`
- **Step size**: $\Delta t = 0.01$ (coordinate time)
- **Total time**: $t_{\text{max}} = 100M$

### Coordinate Systems

- **Schwarzschild coordinates** $(t, r, \theta, \phi)$: Used for calculations
- **Cartesian projection** $(x, y)$: Used for display
  - $x = r\cos\phi$
  - $y = r\sin\phi$

### Viewing Angle Effects

The viewing angle $\theta_{\text{obs}}$ affects the appearance:

- **$\theta = 0°$** (edge-on): Disk is a thin line, maximum Doppler asymmetry
- **$\theta = 45°$** (intermediate): Complex 3D structure visible
- **$\theta = 90°$** (face-on): Nearly circular symmetry, minimal Doppler effect

The angle transformation is:

$$\gamma(\alpha, \theta_0) = \arccos\left(\frac{\cos\alpha}{\sqrt{\cos^2\alpha + \cot^2\theta_0}}\right)$$

---

## Scientific References

1. **Chandrasekhar, S.** (1983). *The Mathematical Theory of Black Holes*. Oxford University Press.
   - Comprehensive treatment of Schwarzschild geometry

2. **Luminet, J.P.** (1979). "Image of a spherical black hole with thin accretion disk". *Astronomy and Astrophysics*, 75, 228-235.
   - First simulations of black hole accretion disk appearance

3. **Event Horizon Telescope Collaboration** (2019). "First M87 Event Horizon Telescope Results". *The Astrophysical Journal Letters*, 875:L1.
   - First direct image of a black hole shadow

4. **Gralla, S.E. & Lupsasca, A.** (2020). "Lensing by Kerr Black Holes". *Physical Review D*, 101, 044031.
   - Modern theory of gravitational lensing by black holes

---

## Implementation Details

### Software Stack
- **Python**: 3.14.3
- **GUI**: PyQt5 (cross-platform interface)
- **Numerics**: NumPy, SciPy (odeint, fsolve, elliptic functions)
- **Visualization**: Matplotlib (2D plotting)
- **Colors**: Seaborn (thermal color palettes)

### Performance
- Isoradial curve calculation: 10-30 seconds (depends on angle)
- Light ray tracing: 0.1-2 seconds (depends on number of photons)
- Uses optimized NumPy arrays for vectorization

### Limitations
- Non-rotating black holes only (Schwarzschild, not Kerr)
- Thin disk approximation (no thickness)
- No relativistic Doppler boosting in colors
- No gravitational redshift visualization
- Static images (no time evolution)

---

## Troubleshooting

**Program won't start**: Make sure you installed all requirements
```powershell
pip install PyQt5 scipy seaborn matplotlib numpy
```

**Too slow**: 
- Reduce number of photons to 10 or less
- The program needs to calculate complex physics, be patient!

**Nothing shows up**: Click "RUN SIMULATION" button after setting your parameters

**Equations not rendering**: If viewing on GitHub, equations should display automatically. For local viewing, use a Markdown viewer that supports LaTeX/MathJax.

---

## Credits

**Physics**: Albert Einstein (General Relativity, 1915)  
**Mathematics**: Carl Gustav Jacobi (Elliptic Functions, 1829)  
**Implementation**: Based on McGill Physics Hackathon 2022  
**Software**: Python 3.14.3 with PyQt5, matplotlib, scipy, numpy, seaborn

**Have fun exploring black holes!** 🌌
