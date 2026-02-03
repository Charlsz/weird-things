
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
#  Physical constants (SI)
# ----------------------------------------------------------------------
G = 6.67430e-11            # m³ kg⁻¹ s⁻²
c = 299_792_458.0          # m s⁻¹

# ----------------------------------------------------------------------
#  Helper functions
# ----------------------------------------------------------------------
def schwarzschild_radius(M_kg: float) -> float:
    """
    Schwarzschild radius rs = 2 G M / c²  (meters)

    Parameters
    ----------
    M_kg : float
        Mass of the black hole in kilograms.

    Returns
    -------
    float
        rs in meters.
    """
    return 2.0 * G * M_kg / c**2


# ----------------------------------------------------------------------
#  Schwarzschild black hole class
# ----------------------------------------------------------------------
class SchwarzschildBlackHole:
    """
    Encapsulates a static, spherically symmetric (Schwarzschild) BH.
    All motion is restricted to the equatorial plane θ = π/2.
    """

    def __init__(self, mass_kg: float):
        """
        Parameters
        ----------
        mass_kg : float
            Black‑hole mass in kilograms.
        """
        self.M = mass_kg                            # kg
        self.rs = schwarzschild_radius(mass_kg)     # meters

    # ------------------------------------------------------------------
    #  Conserved quantities for a circular timelike orbit
    # ------------------------------------------------------------------
    def circular_orbit_E_L(self, r: float):
        """
        Exact energy per unit rest mass (E) and angular momentum per unit
        rest mass (L) of a *stable* circular orbit at radius r (> 3 rs).

        The formulas are (geometrised units G = c = 1, M = rs/2)

            L² = M r / (1 - 3M/r)                (1)
            E² = (r - 2M)² / ( r (r - 3M) )      (2)

        Converting back to SI gives the same dimensionless numbers;
        we multiply by the appropriate powers of c to obtain J kg⁻¹ (energy)
        and m² s⁻¹ (angular momentum).

        Parameters
        ----------
        r : float
            Radius of the circular orbit (meters). Must be > 3 rs.

        Returns
        -------
        E : float
            Conserved energy per unit rest mass (J kg⁻¹).  For a particle at
            rest at infinity, E = c².
        L : float
            Conserved angular momentum per unit rest mass (m² s⁻¹).
        """
        if r <= 3 * self.rs:
            raise ValueError("Circular timelike orbits only exist for r > 3 rs.")

        # Geometrised mass M_geo = G M / c² (has dimensions of length)
        M_geo = G * self.M / c**2
        # Formulas (1)–(2) in geometrised units
        L2_geo = (M_geo * r) / (1.0 - 3.0 * M_geo / r)
        L_geo = np.sqrt(L2_geo)                      # meters
        E2_geo = ((r - 2.0 * M_geo) ** 2) / (r * (r - 3.0 * M_geo))
        E_geo = np.sqrt(E2_geo)                       # dimensionless ( = E / c² )

        # Convert to SI per‑unit‑mass quantities
        L = L_geo * c                               # m² s⁻¹
        E = E_geo * c**2                            # J kg⁻¹
        return E, L

    # ------------------------------------------------------------------
    #  Timelike (massive) geodesic integrator
    # ------------------------------------------------------------------
    def timelike_geodesic(self,
                         r0: float,
                         phi0: float = 0.0,
                         t0: float = 0.0,
                         E: float = None,
                         L: float = None,
                         dr0: float = None,
                         direction: int = -1,
                         tau_max: float = 1e5,
                         max_step: float = None):
        """
        Integrate a massive test‑particle trajectory.

        The equations of motion (using proper time τ as the independent
        variable) are from the Schwarzschild metric with conserved E and L:

            dt/dτ =  E / [ (1 - rs/r) c² ]               (3)
            dφ/dτ =  L / r²                              (4)
            (dr/dτ)² =  E²/c² - (1 - rs/r)(c² + L²/r²)   (5)

        For a given (E, L) the sign of dr/dτ is chosen by ``direction`` or,
        if ``dr0`` is supplied, by its sign.

        Parameters
        ----------
        r0 : float
            Initial radial coordinate (metres).
        phi0, t0 : float, optional
            Initial azimuthal angle (rad) and Schwarzschild coordinate time
            (seconds).  Default to 0.
        E, L : float, optional
            Conserved energy and angular momentum per unit rest mass.
            If omitted the function assumes a *circular* orbit at r0 and
            computes them from ``circular_orbit_E_L``.
        dr0 : float, optional
            Initial proper‑radial velocity dr/dτ (m s⁻¹).  If given, its sign
            overrides ``direction``.
        direction : int, optional
            +1 = outward, -1 = inward.  Used only if ``dr0`` is ``None``.
        tau_max : float, optional
            Maximum proper time to integrate (seconds).  The solver stops
            earlier if the particle crosses the horizon.
        max_step : float, optional
            Upper bound on the integration step‑size (seconds).  If omitted
            it is set automatically to tau_max/5000.

        Returns
        -------
        dict
            Keys ``'tau'``, ``'t'``, ``'r'``, ``'phi'``, ``'E'``, ``'L'``.
        """
        # ------------------------------------------------------------------
        #  Set default conserved quantities (circular orbit) if they are missing
        # ------------------------------------------------------------------
        if (E is None) or (L is None):
            E, L = self.circular_orbit_E_L(r0)
            # For a perfect circular orbit the radial proper velocity is zero
            if dr0 is None:
                dr0 = 0.0

        # ------------------------------------------------------------------
        #  ODE system in proper time τ
        # ------------------------------------------------------------------
        def ode(tau, y):
            """
            y = [t, r, phi]
            Returns dy/dτ.
            """
            t, r, phi = y
            # Guard against entering the singularity
            if r <= self.rs * 1.000001:
                # Freeze the solution – the integrator will stop shortly
                return [0.0, 0.0, 0.0]

            # (1 - rs/r) factor appears often
            factor = 1.0 - self.rs / r

            # ---- dt/dτ (eq. 3) ----
            dt_dtau = E / (factor * c**2)

            # ---- dφ/dτ (eq. 4) ----
            dphi_dtau = L / (r**2)

            # ---- (dr/dτ)² (eq. 5) ----
            radicand = (E**2) / (c**2) - factor * (c**2 + L**2 / (r**2))

            # If the radicand becomes negative the particle has reached a
            # turning point.  We set dr/dτ = 0 there – the sign will be
            # determined by the next step (the solver is adaptive enough).
            if radicand < 0.0:
                dr_dtau = 0.0
            else:
                # Sign choice – from user supplied dr0 or from generic direction
                sign = np.sign(dr0) if dr0 is not None else direction
                dr_dtau = sign * np.sqrt(radicand)

            return [dt_dtau, dr_dtau, dphi_dtau]

        # ------------------------------------------------------------------
        #  Integration
        # ------------------------------------------------------------------
        y0 = [t0, r0, phi0]
        if max_step is None:
            max_step = tau_max / 5000.0

        sol = solve_ivp(
            fun=ode,
            t_span=(0.0, tau_max),
            y0=y0,
            method='RK45',
            max_step=max_step,
            rtol=1e-9,
            atol=1e-12,
        )

        return {
            'tau': sol.t,
            't': sol.y[0],
            'r': sol.y[1],
            'phi': sol.y[2],
            'E': E,
            'L': L,
            'r0': r0,
            'phi0': phi0,
            't0': t0,
        }

    # ------------------------------------------------------------------
    #  Null (photon) geodesic integrator
    # ------------------------------------------------------------------
    def null_geodesic(self,
                      r0: float,
                      impact_parameter: float,
                      phi0: float = 0.0,
                      t0: float = 0.0,
                      direction: int = -1,
                      tau_max: float = 1e5,
                      max_step: float = None):
        """
        Integrate a photon trajectory with a given impact parameter b = L/E.

        In the null case the normalization condition g_{μν}u^μu^ν = 0 leads to

            (dr/dτ)² =  E²/c² - (1 - rs/r) L²/r²         (6)

        We fix a convenient scale by setting E = c (units of velocity); then
        L = b E = b c, where ``b`` is the physical impact parameter
        (metres).

        Parameters
        ----------
        r0 : float
            Starting radius (usually many rs, e.g. 100 rs).
        impact_parameter : float
            b = L/E (metres).  b < 3√3 M (≈ 1.5 rs) leads to capture, larger
            values are deflected.
        direction : int, optional
            -1 = photon moving inward toward the hole, +1 outward.
        tau_max, max_step : float, optional
            Same meaning as in ``timelike_geodesic``.

        Returns
        -------
        dict with the same keys as ``timelike_geodesic`` (plus ``b``).
        """
        # Scale choices (E = c, so L = b c)
        E = c
        L = impact_parameter * c
        b = impact_parameter

        # ------------------------------------------------------------------
        #  ODE system (null)
        # ------------------------------------------------------------------
        def ode(tau, y):
            t, r, phi = y
            if r <= self.rs * 1.000001:
                return [0.0, 0.0, 0.0]

            factor = 1.0 - self.rs / r
            dt_dtau = E / (factor * c**2)            # same as massive case
            dphi_dtau = L / (r**2)

            radicand = (E**2) / (c**2) - factor * (L**2) / (r**2)
            dr_dtau = direction * np.sqrt(radicand) if radicand > 0 else 0.0

            return [dt_dtau, dr_dtau, dphi_dtau]

        y0 = [t0, r0, phi0]
        if max_step is None:
            max_step = tau_max / 5000.0

        sol = solve_ivp(
            fun=ode,
            t_span=(0.0, tau_max),
            y0=y0,
            method='RK45',
            max_step=max_step,
            rtol=1e-9,
            atol=1e-12,
        )

        return {
            'tau': sol.t,
            't': sol.y[0],
            'r': sol.y[1],
            'phi': sol.y[2],
            'E': E,
            'L': L,
            'b': b,
            'r0': r0,
        }

    # ------------------------------------------------------------------
    #  Plotting utilities
    # ------------------------------------------------------------------
    def plot_trajectory(self,
                       traj: dict,
                       ax=None,
                       show_horizon: bool = True,
                       label: str = None,
                       **plot_kwargs):
        """
        Plot a trajectory (r, φ) in Cartesian coordinates.

        Parameters
        ----------
        traj : dict
            Output dictionary from ``timelike_geodesic`` or ``null_geodesic``.
            Must contain ``'r'`` and ``'phi'`` arrays.
        ax : matplotlib.axes.Axes, optional
            Axis to draw on.  If omitted, a new figure+axis is created.
        show_horizon : bool, optional
            Draw the event horizon as a filled black disc.
        label : str, optional
            Legend entry.
        plot_kwargs : dict
            Additional arguments passed to ``ax.plot``.

        Returns
        -------
        matplotlib.axes.Axes
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))

        # Convert polar → Cartesian
        x = traj['r'] * np.cos(traj['phi'])
        y = traj['r'] * np.sin(traj['phi'])

        ax.plot(x, y, label=label, **plot_kwargs)

        if show_horizon:
            horizon = plt.Circle((0, 0), self.rs, color='black', zorder=10)
            ax.add_artist(horizon)

        ax.set_aspect('equal')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_title('Geodesic around a Schwarzschild Black Hole')
        if label is not None:
            ax.legend()
        return ax


# ----------------------------------------------------------------------
#  Example usage (run the file directly)
# ----------------------------------------------------------------------
if __name__ == '__main__':
    # ------------------------------------------------------------------
    #  Choose a black‑hole mass
    # ------------------------------------------------------------------
    M_sun = 1.98847e30                    # kg
    M_bh = 10 * M_sun                     # 10‑solar‑mass BH (typical stellar BH)
    bh = SchwarzschildBlackHole(M_bh)

    print(f"Mass = {M_bh:.3e} kg")
    print(f"Schwarzschild radius = {bh.rs / 1e3:.3f} km")

    # ------------------------------------------------------------------
    #  1) Massive particle on a slightly perturbed circular orbit
    # ------------------------------------------------------------------
    r_circ = 12.0 * bh.rs                # radius of (approx) circular orbit

    # Exact (E, L) for a perfect circle at r_circ
    E_circ, L_circ = bh.circular_orbit_E_L(r_circ)
    print("\nCircular orbit at {:.0f} rs:".format(r_circ / bh.rs))
    print(f"  Energy per unit mass = {E_circ / c**2:.6f} c²")
    print(f"  Angular momentum per unit mass = {L_circ:.3e} m² s⁻¹")

    # Perturb L a little to make the orbit elliptical (precessing)
    L_pert = 0.95 * L_circ

    # Integrate for a proper time of ~ 1 s (for a 10‑M⊙ BH the dynamical time is ≈ 10⁻⁴ s,
    # so 1 s gives many orbits)
    traj_massive = bh.timelike_geodesic(
        r0=r_circ,
        phi0=0.0,
        E=E_circ,
        L=L_pert,
        dr0=0.0,
        direction=-1,
        tau_max=2.0,               # seconds of proper time
    )
    # Plot the orbit
    bh.plot_trajectory(traj_massive, color='tab:blue',
                        label='massive particle', linewidth=1.5)
    plt.show()

    # ------------------------------------------------------------------
    #  2) Photon trajectories with different impact parameters
    # ------------------------------------------------------------------
    # Starting far away (100 rs)
    r_start = 100.0 * bh.rs

    # Critical impact parameter for capture is b_crit = 3*sqrt(3) M_geo = 1.5*sqrt(3) rs ≈ 2.598 rs
    b_crit = 1.5 * np.sqrt(3) * bh.rs
    print("\nCritical photon impact parameter (capture) = {:.3f} rs".format(b_crit / bh.rs))

    # Choose three impact parameters: below, at, and above the critical one
    impact_params = [0.8 * b_crit, 1.0 * b_crit, 1.5 * b_crit]

    fig, ax = plt.subplots(figsize=(6, 6))
    colors = ['red', 'orange', 'green']
    for b, col in zip(impact_params, colors):
        traj_photon = bh.null_geodesic(
            r0=r_start,
            impact_parameter=b,
            direction=-1,
            tau_max=10.0,               # enough to see the turn‑around or capture
        )
        label = f"b = {b / bh.rs:.2f} rs"
        bh.plot_trajectory(traj_photon, ax=ax, label=label,
                            color=col, linewidth=1.2, show_horizon=False)
    # Re‑draw the horizon once (so it sits on top)
    bh.plot_trajectory(traj_photon, ax=ax, show_horizon=True, color='k', linewidth=0)

    ax.set_title('Photon deflection – impact parameter vs. capture')
    ax.legend()
    plt.show()
