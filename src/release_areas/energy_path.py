"""Per-cell quantities for path-energy crack arrest (hypothesis, not validated).

A crack carries slab tension N along its path; it advances into a cell while
the energy that tension can release, G_avail, covers the weak layer's G_c.
See docs/release_area_methods.md §3.

Units: SI (Pa, m, N/m, J/m2); angles in degrees.
"""
import numpy as np

from release_areas import config
from release_areas.arrest_indices import (
    G_ACC,
    gravitational_shear,
    plane_stress_modulus,
    residual_shear,
    slab_energy_cap,
    weak_layer_fracture_energy,
)


def tau_net(rho, h, psi_deg, phi_deg=None):
    """Net driving shear rho g h (sin psi - cos psi tan phi) (Pa); <= 0 at psi <= phi.
    Written as rho g h sin(psi - phi) / cos(phi), which is exactly 0 at psi = phi."""
    phi_deg = config.PHI_DEG if phi_deg is None else phi_deg
    return (rho * G_ACC * h * np.sin(np.radians(psi_deg - phi_deg))
            / np.cos(np.radians(phi_deg)))


def g_c(tau_p, K_wl, delta=None, tau_r=0.0):
    """Weak-layer fracture energy (tau_p - tau_r)^2 (1 + delta) / (2 K_wl) (J/m2)."""
    delta = config.ENERGY_DELTA if delta is None else delta
    return weak_layer_fracture_energy(tau_p, K_wl, delta, tau_r)


def g_slab(sigma_t, h, E, nu=None):
    """Slab energy cap sigma_t^2 h / (2 E') (J/m2)."""
    nu = config.NU if nu is None else nu
    return slab_energy_cap(sigma_t, h, E, nu)


def n_max(sigma_t, h):
    """Slab tension cap sigma_t h (N/m)."""
    return sigma_t * h


def g_avail(N, E, h, nu=None, R=None):
    """Energy available from slab tension N: R N^2 / (2 E' h) (J/m2)."""
    nu = config.NU if nu is None else nu
    R = config.R_FIT if R is None else R
    return R * np.asarray(N, float)**2 / (2.0 * plane_stress_modulus(E, nu) * h)


def e_column(relation=None):
    """Feature column holding E_eff under `relation`."""
    relation = config.E_RELATION if relation is None else relation
    return 'E_eff' if relation == config.E_RELATION else f'E_eff__{relation}'


def cell_energy(row, psi_deg, sigma_t_col=None, e_relation=None, delta=None,
                G_wl=None, use_residual=None):
    """tau_net, G_c, G_slab_path, N_max for one cluster row; NaN where inputs are missing.

    Defaults come from config. K_wl = G_wl / wl_thickness (D_WL_FALLBACK if absent).
    """
    sigma_t_col = config.SIGMA_T_COLUMN if sigma_t_col is None else sigma_t_col
    e_relation = config.E_RELATION if e_relation is None else e_relation
    delta = config.ENERGY_DELTA if delta is None else delta
    G_wl = config.G_WL if G_wl is None else G_wl
    use_residual = config.ENERGY_USE_RESIDUAL if use_residual is None else use_residual

    def get(k):
        try:
            v = float(row.get(k, np.nan))
        except (TypeError, ValueError):
            return np.nan
        return v

    rho, h = get('slab_density'), get('slab_thickness')
    E, sig = get(e_column(e_relation)), get(sigma_t_col)
    tau_p, D_wl = get('wl_shear_strength'), get('wl_thickness')
    if not (np.isfinite(D_wl) and D_wl > 0):
        D_wl = config.D_WL_FALLBACK
    K_wl = G_wl / D_wl
    psi = float(psi_deg)
    tr = float(residual_shear(rho, h, psi, config.PHI_DEG)) if use_residual else 0.0
    return {
        'tau_net': float(tau_net(rho, h, psi)),
        'G_c': float(g_c(tau_p, K_wl, delta, tr)),
        'G_slab_path': float(g_slab(sig, h, E)),
        'N_max': float(n_max(sig, h)),
    }
