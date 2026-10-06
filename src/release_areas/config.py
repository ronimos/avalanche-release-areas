"""
config.py — Default physical and operational parameters for release_areas.

All values can be overridden via keyword arguments to individual functions
or via a TOML config file passed to run_scenarios.py.
"""

from __future__ import annotations

# -----------------------------------------------------------------------
# Snowpack / weak-layer material
# -----------------------------------------------------------------------
G_WL     = 0.2e6   # weak-layer shear modulus (Pa)
NU       = 0.3     # slab Poisson's ratio
PHI_DEG  = 27.0    # snow friction angle (degrees)
DELTA    = 1.0     # Meloche softening coefficient δ

# -----------------------------------------------------------------------
# BFS crack propagation
# -----------------------------------------------------------------------
STAUCHWALL_DEG        = 28.0   # slope threshold for downslope arrest
GAUME_ASPECT_CAP      = 2.5    # max cross-slope width / A_ca
MIN_POLYGON_AREA      = 200.0  # m²; discard smaller polygons
TAU_G_ABS_FLOOR       = 350.0  # Pa; absolute driving-stress floor
MIN_PROPAGATION_SLAB  = 0.50   # m; thin-slab arrest threshold
MIN_PROPAGATION_LAMBDA = 0.1   # m; compliant-slab arrest threshold
LAMBDA_DROP_FACTOR    = 0.20   # fractional softening arrest threshold
LAMBDA_RISE_FACTOR    = 0.30   # fractional stiffening arrest threshold
THICKNESS_DROP_FACTOR = 0.20   # fractional thinning arrest threshold
THICKNESS_RISE_FACTOR = 0.30   # fractional thickening arrest threshold

# -----------------------------------------------------------------------
# Scenario sweep
# -----------------------------------------------------------------------
N_TOP_TRIGGERS   = 5       # number of trigger clusters to evaluate
SIZE_FACTORS     = [0.75, 0.875, 1.0, 1.125, 1.25]
DEPTH_PCTS       = [0.25, 0.50, 0.75]   # slab-depth percentiles for density

# -----------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------
SCENARIO_ID_PREFIX = "scenario"
