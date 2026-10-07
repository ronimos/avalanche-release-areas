"""Crack-arrest indices for dry-slab avalanches.

Sources: Meloche et al. (2025, JGR Earth Surface), Meloche et al. (ISSW 2026),
Meloche et al. (2026, TARP report). The energy-cap functions are a hypothesis
(see docs/release_area_methods.md §3), not published results.

Units: SI (Pa, m, kg/m3, J/m2); angles in degrees.
"""
import numpy as np

G_ACC = 9.81


# ---------------------------------------------------------------------------
# Material and loading helpers
# ---------------------------------------------------------------------------
def plane_stress_modulus(E, nu=0.3):
    """E' = E / (1 - nu^2)."""
    return E / (1.0 - nu**2)


def shear_modulus(E, nu=0.3):
    """G = E / (2 (1 + nu))."""
    return E / (2.0 * (1.0 + nu))


def weak_layer_stiffness(G_wl, D_wl):
    """K_wl = G_wl / D_wl (Pa/m)."""
    return G_wl / D_wl


def elastic_length(E, h, D_wl, G_wl, nu=0.3):
    """Upslope (mode II) Lambda = sqrt(E h D_wl / ((1 - nu^2) G_wl))."""
    return np.sqrt(E * h * D_wl / ((1.0 - nu**2) * G_wl))


def elastic_length_cross(E, h, D_wl, G_wl, nu=0.3):
    """Cross-slope (mode III) Lambda using slab shear modulus G instead of E'.
    Hypothesis: Lambda_III / Lambda_II = sqrt((1 - nu) / 2)."""
    return np.sqrt(shear_modulus(E, nu) * h * D_wl / G_wl)


def gravitational_shear(rho, h, psi_deg):
    """tau_g = rho g h sin(psi)."""
    return rho * G_ACC * h * np.sin(np.radians(psi_deg))


def residual_shear(rho, h, psi_deg, phi_deg=27.0):
    """tau_r = rho g h cos(psi) tan(phi)."""
    return rho * G_ACC * h * np.cos(np.radians(psi_deg)) * np.tan(np.radians(phi_deg))


def wave_speeds(E, rho, nu=0.3):
    """(c_s, c_p) with c_p = sqrt(E'/rho)."""
    return np.sqrt(shear_modulus(E, nu) / rho), np.sqrt(plane_stress_modulus(E, nu) / rho)


# ---------------------------------------------------------------------------
# Strength-based indices (published scalings)
# ---------------------------------------------------------------------------
def tension_gradient(rho, psi_deg, phi_deg=27.0):
    """k_f = rho g sin(psi) (1 - tan(phi)/tan(psi)); 0 when psi <= phi."""
    psi, phi = np.radians(psi_deg), np.radians(phi_deg)
    k = rho * G_ACC * np.sin(psi) * (1.0 - np.tan(phi) / np.tan(psi))
    return np.maximum(k, 0.0)


def dynamic_factor(E, rho, nu=0.3, speed_ratio=1.6):
    """k_x / k_f = c_p^2 / (c_p^2 + adot^2), adot = speed_ratio * c_s."""
    c_s, c_p = wave_speeds(E, rho, nu)
    adot = speed_ratio * c_s
    return c_p**2 / (c_p**2 + adot**2)


def dynamic_gradient(k_f, E, rho, nu=0.3, speed_ratio=1.6):
    """k_x = k_f c_p^2 / (c_p^2 + adot^2)."""
    return k_f * dynamic_factor(E, rho, nu, speed_ratio)


# Crack-speed regimes. Upslope (mode II) cracks run supershear; mode III
# (antiplane) cracks cannot exceed the shear wave speed.
#
# On the 1.6: Meloche et al. (2025) write the supershear speed as 1.6 c_s but
# identify it with c_p -- "A horizontal dashed line is set at 1.6 c_s ~= sqrt(E'/rho),
# which is the supershear speed (Trottet et al., 2022)" -- and plot "until
# c_p ~= 1.6 c_s". So 1.6 is their rounding of c_p/c_s = sqrt(2/(1-nu)) = 1.690
# at nu = 0.3. Taking adot = c_p exactly gives k_x/k_f = 1/2, which is what the
# paper's own Fig. 3 caption states; the literal 1.6 gives 0.527.
# We keep the literal 1.6 so published figures reproduce; cp_over_cs() below
# gives the exact alternative. The choice moves mode3_length_ratio by ~5%
# (0.712 -> 0.675).
MODE2_SPEED_RATIO = 1.6   # adot / c_s upslope, Meloche et al. (2025) as written
MODE3_SPEED_RATIO = 1.0   # adot / c_s cross-slope cap, Broberg (1989)


def cp_over_cs(nu=0.3):
    """c_p / c_s = sqrt(2/(1-nu)); 1.690 at nu=0.3.

    The exact supershear bound that Meloche et al. round to 1.6. Pass as
    speed_ratio_along to evaluate() for the self-consistent k_x/k_f = 1/2.
    """
    return np.sqrt(2.0 / (1.0 - nu))


def mode3_length_ratio(E, rho, nu=0.3,
                       speed_ratio_along=MODE2_SPEED_RATIO,
                       speed_ratio_cross=MODE3_SPEED_RATIO):
    """L_dyn,III / L_dyn,II = k_x,II / k_x,III.

    A slower crack builds slab tension faster per unit advance, so the mode III
    speed cap shortens the distance to first slab fracture. Since
    c_s^2/c_p^2 = (1-nu)/2, E and rho cancel and this reduces to a function of
    nu and the two speed ratios alone: 0.712 at nu=0.3.
    """
    return (dynamic_factor(E, rho, nu, speed_ratio_along)
            / dynamic_factor(E, rho, nu, speed_ratio_cross))


def tensile_length(sigma_t, k):
    """sigma_t / k; inf if k <= 0."""
    k = np.asarray(k, dtype=float)
    with np.errstate(divide="ignore"):
        return np.where(k > 0, sigma_t / k, np.inf)


def arrest_length(tau_g, theta, Lam, sigma_t, delta=1.0, L_t=None, C=0.045):
    """JGR Eq. 20: A_ca / L_t = C * X * sqrt(sigma_t / tau_g),
    X = tau_g / (theta Lam sqrt(1 + delta)). C = 0.045 is a two-run fit (Fig. 8).
    Returns (A_ca or A_ca/L_t, X)."""
    X = tau_g / (theta * Lam * np.sqrt(1.0 + delta))
    rel = C * X * np.sqrt(sigma_t / tau_g)
    return (rel * L_t if L_t is not None else rel), X


def release_size_index(a_t, a_c, a_sc, arrest_ratio=6.5):
    """ISSW 2026 preliminary RSI: 1 arrest, 2 sub-Rayleigh en-echelon, 3 supershear."""
    if a_t / a_sc >= 1.0:
        return 3
    return 1 if a_t / a_c < arrest_ratio else 2


def shear_gradient(x, tau_p):
    """theta (Pa/m): linear-fit slope of tau_p along an upslope transect."""
    return float(np.polyfit(np.asarray(x, float), np.asarray(tau_p, float), 1)[0])


# ---------------------------------------------------------------------------
# Energy-cap framework (hypothesis)
# ---------------------------------------------------------------------------
R_FIT = 0.48  # fitted from Meloche et al. 2025 Fig. 8 (two runs); see docs/release_area_methods.md §3


def slab_energy_cap(sigma_t, h, E, nu=0.3):
    """Energy a slab segment loaded to sigma_t can deliver per unit crack advance:
    G_slab = sigma_t^2 h / (2 E')  (J/m2)."""
    return sigma_t**2 * h / (2.0 * plane_stress_modulus(E, nu))


def slab_energy_cap_cross(tau_flank, h, E, nu=0.3):
    """Cross-slope analogue: G_slab,III = tau_flank^2 h / (2 G).
    tau_flank is a slab flank (shear) strength; no published parameterization."""
    return tau_flank**2 * h / (2.0 * shear_modulus(E, nu))


def weak_layer_fracture_energy(tau_p, K_wl, delta=1.0, tau_r=0.0):
    """G_c = (tau_p - tau_r)^2 (1 + delta) / (2 K_wl)  (J/m2)."""
    return (np.asarray(tau_p, float) - tau_r)**2 * (1.0 + delta) / (2.0 * K_wl)


def energy_ratio(tau_p, K_wl, sigma_t, h, E, delta=1.0, tau_r=0.0, nu=0.3):
    """R = G_c(tau_p) / G_slab. Hypothesis: arrest when R >= R_crit."""
    return weak_layer_fracture_energy(tau_p, K_wl, delta, tau_r) / slab_energy_cap(sigma_t, h, E, nu)


def critical_strength(sigma_t, h, E, K_wl, delta=1.0, tau_r=0.0, R=R_FIT, nu=0.3):
    """tau_p* where G_c(tau_p*) = R * G_slab.
    R must be a number; pass R='dynamic' to evaluate() to resolve it via k_x/k_f first."""
    G_cap = R * slab_energy_cap(sigma_t, h, E, nu)
    return tau_r + np.sqrt(2.0 * K_wl * G_cap / (1.0 + delta))


def arrest_length_energy(tau_p0, tau_p_star, theta):
    """A_ca = (tau_p* - tau_p0) / theta for a linear strength ramp.
    0 if tau_p0 >= tau_p*; inf if theta <= 0."""
    if theta <= 0:
        return np.inf
    return max(tau_p_star - tau_p0, 0.0) / theta


# ---------------------------------------------------------------------------
# Combined evaluation
# ---------------------------------------------------------------------------
def evaluate(rho, h, psi_deg, E, sigma_t, D_wl, G_wl, theta, a_c=None, a_sc=None,
             nu=0.3, phi_deg=27.0, delta=1.0, C=0.045, tau_p0=None, R=R_FIT,
             use_residual=False, speed_ratio_along=MODE2_SPEED_RATIO,
             speed_ratio_cross=MODE3_SPEED_RATIO):
    """All indices for one slab / weak-layer configuration.
    R: critical energy ratio (number) or 'dynamic' to use k_x/k_f.
    speed_ratio_along / speed_ratio_cross: crack speed as a multiple of c_s for
    mode II (supershear) and mode III (capped at c_s) propagation."""
    Lam = elastic_length(E, h, D_wl, G_wl, nu)
    K_wl = weak_layer_stiffness(G_wl, D_wl)
    tau_g = gravitational_shear(rho, h, psi_deg)
    tau_r = residual_shear(rho, h, psi_deg, phi_deg)
    k_f = float(tension_gradient(rho, psi_deg, phi_deg))
    k_x = float(dynamic_gradient(k_f, E, rho, nu, speed_ratio_along))
    k_x_cross = float(dynamic_gradient(k_f, E, rho, nu, speed_ratio_cross))
    L_t = float(tensile_length(sigma_t, k_f))
    L_dyn = float(tensile_length(sigma_t, k_x))
    L_dyn_cross = float(tensile_length(sigma_t, k_x_cross))

    out = dict(Lambda=Lam, Lambda_cross=elastic_length_cross(E, h, D_wl, G_wl, nu),
               K_wl=K_wl, tau_g=tau_g, tau_r=tau_r, k_f=k_f, k_x=k_x,
               k_x_cross=k_x_cross, L_t=L_t, L_dyn=L_dyn,
               L_dyn_cross=L_dyn_cross,
               mode3_length_ratio=float(mode3_length_ratio(
                   E, rho, nu, speed_ratio_along, speed_ratio_cross)),
               sustained=psi_deg > phi_deg)

    if np.isfinite(L_t) and theta > 0:
        out["A_ca"], out["X"] = arrest_length(tau_g, theta, Lam, sigma_t, delta, L_t, C)

    R_val = dynamic_factor(E, rho, nu) if R == "dynamic" else R
    tr = tau_r if use_residual else 0.0
    out["G_slab"] = slab_energy_cap(sigma_t, h, E, nu)
    out["R_used"] = R_val
    out["tau_p_star"] = critical_strength(sigma_t, h, E, K_wl, delta, tr, R_val, nu)
    if tau_p0 is not None:
        out["G_c0"] = float(weak_layer_fracture_energy(tau_p0, K_wl, delta, tr))
        out["R0"] = out["G_c0"] / out["G_slab"]
        out["A_ca_energy"] = arrest_length_energy(tau_p0, out["tau_p_star"], theta)

    if a_c and a_sc and np.isfinite(L_dyn):
        out["RSI"] = release_size_index(L_dyn, a_c, a_sc)
    return out
