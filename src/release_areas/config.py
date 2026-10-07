"""
config.py — Default physical and operational parameters for release_areas.

This module is the single source of truth for every tunable constant in the
package. Other modules import from here rather than redefining literals, so
changing a value in one place changes it everywhere.

Function keyword arguments and generate_scenarios CLI flags override these
defaults; the values below are only the fallbacks.
"""

from __future__ import annotations

# -----------------------------------------------------------------------
# Coordinate reference system
# -----------------------------------------------------------------------
# The Little Professor rasters are NAD83(2011) / UTM 13N + NAVD88 height,
# whose horizontal component is EPSG:6342 (NOT EPSG:32613 / WGS84 UTM 13N).
# Vector inputs are reprojected into this CRS so they register to the grid.
RASTER_EPSG = 6342

# -----------------------------------------------------------------------
# Snowpack / weak-layer material
# -----------------------------------------------------------------------
G_WL     = 0.2e6   # weak-layer shear modulus (Pa); Meloche et al. (2025) Table 1
NU       = 0.3     # slab Poisson's ratio
PHI_DEG  = 27.0    # snow friction angle (degrees)
DELTA    = 1.0     # Meloche softening coefficient δ; see methods §6 on δ sensitivity
D_WL_FALLBACK = 0.04   # m; weak-layer thickness used only when wl_thickness is NaN
L_SS     = 20.0    # m; steady-state propagation length, Meloche et al. (2025)

# -----------------------------------------------------------------------
# Feature extraction (compute_meloche_features)
# -----------------------------------------------------------------------
K_NEIGHBORS        = 6      # k-nearest cluster centroids used for the θ gradient
TAU_G_FEATURE_FLOOR = 50.0  # Pa; below this the arrest scaling laws are not evaluated
THETA_MIN          = 1e-6   # Pa/m; θ below this is treated as no usable gradient

# -----------------------------------------------------------------------
# Trigger selection filter chain (generate_scenarios)
# -----------------------------------------------------------------------
MIN_TAU_G          = 40.0   # Pa; minimum gravitational driving shear stress
MAX_SK38           = 1.0    # Sk38 must be below this to count as unstable
MIN_SLAB_THICKNESS = 0.5    # m
MAX_SLAB_THICKNESS = 2.0    # m; use 1.5 for skier-triggered scenarios
SLOPE_MARGIN_DEG   = 2.0    # degrees added to stauchwall angle for the trigger gate
ELEVATION_PCTILE   = 50     # keep clusters at or above this elevation percentile

# -----------------------------------------------------------------------
# BFS crack propagation
# -----------------------------------------------------------------------
STAUCHWALL_DEG        = 28.0   # slope threshold for downslope arrest
GAUME_ASPECT_CAP      = 2.5    # max cross-slope width / A_ca
MODE3_SCALE           = 1.5    # mode III lateral arrest multiplier
MIN_POLYGON_AREA      = 200.0  # m²; discard smaller polygons
TAU_G_ABS_FLOOR       = 350.0  # Pa; absolute driving-stress floor
MIN_PROPAGATION_SLAB  = 0.50   # m; thin-slab arrest threshold
MIN_PROPAGATION_LAMBDA = 0.1   # m; compliant-slab arrest threshold
# Λ continuity thresholds, along-slope (mode II: upslope + downslope).
LAMBDA_DROP_FACTOR    = 0.20   # fractional softening arrest threshold
LAMBDA_RISE_FACTOR    = 0.30   # fractional stiffening arrest threshold

# Λ continuity thresholds, cross-slope (mode III / flanks). Separate knobs so
# the flanks can be tuned independently once Λ_III is known; identical to the
# along-slope values for now, so splitting the code changed no results.
LAMBDA_CROSS_DROP_FACTOR = 0.20
LAMBDA_CROSS_RISE_FACTOR = 0.30

THICKNESS_DROP_FACTOR = 0.20   # fractional thinning arrest threshold
THICKNESS_RISE_FACTOR = 0.30   # fractional thickening arrest threshold

# ---- Mode III (cross-slope) elastic length -----------------------------
# The cross-slope elastic length is a different and SMALLER quantity than the
# along-slope (mode II) Λ. Johan Gaume's antiplane formula is not published
# yet, so propagation currently uses the mode II Λ in every direction and the
# cross-slope branch is a no-op.
#
# When the formula arrives:
#   1. update `arrest_indices.elastic_length_cross()` with the real expression;
#   2. flip USE_MODE3_LAMBDA to True.
# Nothing else should need to change — if it does, the directional split in
# release_geometry._qualifies() has drifted and needs fixing instead.
#
# The `Lambda_cross` column is already written to the meloche CSV from the
# current (unpublished) derivation, so the data side is ready; this flag only
# controls whether the BFS reads it. Note the shipped Jan 18 CSV predates the
# column — regenerate it before the flag can do anything.
#
# IMPORTANT, measured: flipping this flag has NO effect if Λ_III is a constant
# multiple of Λ_II. The continuity test compares the *relative* change
# (Λ_nbr − Λ_cur)/Λ_cur, which is scale-invariant, and the absolute
# MIN_PROPAGATION_LAMBDA floor is never reached (min Λ_II = 0.54 m, so
# min Λ_III = 0.32 m, both far above 0.1 m). The current derivation
# Λ_III/Λ_II = √((1−ν)/2) IS such a constant, so it provably changes nothing.
# A real Λ_III only matters if it has different *spatial structure* — i.e.
# depends on material properties that vary across clusters differently than E′
# does. Verified: a spatially varying Λ_cross moves every release area.
USE_MODE3_LAMBDA = False

# ---- Mode III (cross-slope) crack-speed cap ----------------------------
# Upslope (mode II) cracks run supershear at ~1.6 c_s; mode III (antiplane)
# cracks cannot exceed c_s (Broberg 1989). Since the dynamic tension gradient
# is k_x = k_f c_p^2/(c_p^2 + adot^2), a slower crack builds slab tension
# FASTER per unit advance, so the distance to first slab fracture is shorter
# cross-slope: L_dyn,III / L_dyn,II = 0.712 at nu = 0.3.
#
# Unlike USE_MODE3_LAMBDA this is not scale-invariant — it enters as an
# absolute distance cap on lateral propagation, so it does change results.
# The ratio depends only on nu and the two speed ratios (E and rho cancel).
MODE2_SPEED_RATIO   = 1.6    # adot / c_s upslope, Meloche et al. (2025)
MODE3_SPEED_RATIO   = 1.0    # adot / c_s cross-slope cap, Broberg (1989)
USE_MODE3_SPEED_CAP = True   # apply the ratio to the lateral distance cap
BFS_K_NEIGHBOURS      = 8      # cluster adjacency degree for the flood-fill
MAX_BFS_CLUSTERS      = 500    # safety cap on region size; see methods §5

# Meloche per-direction arrest test inside the BFS. Disabled by default: the
# flood-fill uses the Λ/thickness discontinuity heuristics plus TAU_G_ABS_FLOOR
# instead. See methods §5 — enabling this changes every release polygon.
USE_MELOCHE_ARREST = False

# -----------------------------------------------------------------------
# Scenario sweep
# -----------------------------------------------------------------------
N_TOP_TRIGGERS   = 5       # number of trigger clusters to evaluate
SIZE_FACTORS     = [1.0]   # release size multipliers; one polygon per trigger
