"""
compute_indices_from_xsnow.py — Example: generate release_areas input CSVs from xsnow data.

This script shows the complete path from a SNOWPACK xsnow Dataset (zarr or .pro files)
to the two CSVs consumed by generate_scenarios:

    all_start_zone_features_<DATE>.csv   — per-cluster WL/slab features
    meloche_features_all_<DATE>.csv      — per-cluster arrest indices (Π₁, A_ca, θ, …)

ILLUSTRATIVE ONLY — the paths below are placeholders and do not exist. It writes
to OUT_DIR, which defaults to a scratch directory and deliberately NOT to
data/little_prof/features: those are the committed reference CSVs, and this
script uses a different MIN_DEPTH_CM and does not pin the cluster set or the
`group` labels, so it would not reproduce them. To regenerate the reference
pair, use examples/regenerate_jan18_reference_csvs.py instead.

Run from the avalanche-release-areas repo root after editing the five paths at the top.

Requirements (not in release-areas pyproject.toml; installed separately in avachain env):
    xsnow, xarray, dask
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import xarray as xr

# release_areas package (pip install -e . from repo root)
from release_areas.snowpack_features import profile_features, compute_meloche_features

# ---------------------------------------------------------------------------
# Paths — edit these
# ---------------------------------------------------------------------------
PRO_DIR     = Path("/home/snowpath/avachain/pro/little_prof")   # cluster_*.pro files
ZARR_PATH   = Path("/home/snowpath/avachain/zarr/little_prof")  # zarr store (faster)
CLUSTER_MAP = Path("data/little_prof/spatial/cluster_map.tif")
DEM_PATH    = Path("data/little_prof/dem_1m.tif")
# NOT data/little_prof/features — see the module docstring.
OUT_DIR     = Path("outputs/xsnow_example_features")

EVENT_DATE  = pd.Timestamp("2026-01-18 12:00")   # noon on the analysis date

MIN_DEPTH_CM = 10.0   # minimum WL burial depth to count (cm)

# ---------------------------------------------------------------------------
# 1. Load xsnow Dataset
# ---------------------------------------------------------------------------
# The Dataset has three dimensions:
#   location  — one entry per cluster_*.pro file; name is the station ID string
#               e.g. "cluster_0001_cluster_0001"; cluster integer ID = int(name.split('_')[-1])
#   time      — hourly or noon snapshots
#   layer     — SNOWPACK finite-element layers (bottom → top)
#
# Key per-layer variables:
#   z               height of layer mid-point (m, negative from snow surface; basal WL is most negative)
#   grain_type      SNOWPACK grain code  (4xx = FC, 5xx = DH, 2xx = DF, 3xx = RG, 7xx+ = crust/ice)
#   density         kg/m³
#   hand_hardness   SNOWPACK F/4F/1F/… scale as float
#   grain_size      m
#   shear_strength  Pa   (SNOWPACK variable 0508)
#   sk38            Sk38 skier stability index (dimensionless)
#   ssi             structural stability index
#   sn38            SN38 natural stability index
#   critical_cut_length   r_c (m); SNOWPACK variable 0606

try:
    import xsnow
    if ZARR_PATH.exists():
        print(f"Loading from zarr: {ZARR_PATH}")
        _dr = xr.open_zarr(str(ZARR_PATH))
        ds = xsnow.xsnowDataset(_dr)
    else:
        print(f"Loading from .pro files: {PRO_DIR}")
        ds = xsnow.read(str(PRO_DIR), lazy=False, n_cpus_use=4)
    print(f"  dims: { {k: v for k, v in ds.sizes.items()} }")
except ImportError:
    raise SystemExit(
        "xsnow is not installed in this environment.\n"
        "Run this script from the avachain conda/venv environment, or\n"
        "install xsnow there and then: pip install -e /path/to/avalanche-release-areas"
    )

location_names = ds.coords["location"].values

# ---------------------------------------------------------------------------
# 2. Select the analysis timestep and extract per-cluster features
# ---------------------------------------------------------------------------
# Select noon on EVENT_DATE (nearest available timestep).
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ds_t = ds.sel(time=EVENT_DATE, method="nearest").compute()

n_locs = ds_t.sizes.get("location", 1)
print(f"Extracting features for {n_locs} locations at {EVENT_DATE.date()}")

rows = []
for i in range(n_locs):
    # Squeeze to a single-location Dataset (dims: layer only)
    ds_loc = ds_t.isel(location=i) if "location" in ds_t.dims else ds_t

    try:
        feat = profile_features(ds_loc, min_depth_cm=MIN_DEPTH_CM)
    except Exception as exc:
        warnings.warn(f"profile_features failed for location {i}: {exc}")
        feat = {}

    # Recover the integer cluster ID from the location name
    loc_name = str(location_names[i])
    try:
        cid = int(loc_name.split("_")[-1])
    except ValueError:
        cid = i
    feat["cluster_id"] = cid
    rows.append(feat)

features_df = pd.DataFrame(rows).set_index("cluster_id")
print(f"  Features: {len(features_df)} clusters, {features_df.notna().sum(axis=1).mean():.0f} non-null cols avg")

# ---------------------------------------------------------------------------
# 3. Load the cluster map and DEM
# ---------------------------------------------------------------------------
with rasterio.open(CLUSTER_MAP) as src:
    cluster_map = src.read(1).astype(float)
    transform   = src.transform

with rasterio.open(DEM_PATH) as src:
    dem = src.read(1).astype(float)
    dem[dem <= -9999] = np.nan

# ---------------------------------------------------------------------------
# 4. Compute Meloche et al. (2025) crack-arrest indices
# ---------------------------------------------------------------------------
# compute_meloche_features() adds:
#   slope_angle   — per-cluster mean slope from DEM (degrees)
#   tau_g         — gravitational shear stress (Pa)
#   theta         — WL shear-strength gradient (Pa/m), by config.THETA_ESTIMATOR:
#                   'plane_fit' (default) fits a plane to τp within L_ss;
#                   'knn' takes the mean |Δτp|/d over k=6 nearest neighbours
#   theta_grad_east / theta_grad_north — the θ gradient vector ('plane_fit' only)
#   Lambda_cross  — cross-slope (mode III) elastic length (m)
#   Pi1_elastic   — Π₁ = τ_g / (θ Λ √(1+δ))
#   Pi2_brittle   — Π₂ = Π₁ √(σ_t/τ_g)
#   Lambda        — upslope elastic length (m)
#   L_t           — quasi-static tensile length (m)
#   tau_p         — WL shear strength carried through from profile_features (Pa)
#   A_ca_brittle  — brittle-slab arrest length from Eq. 20, C=0.045 (m)
#   G_slab        — slab energy cap (J/m²)
#   tau_p_star    — critical WL strength for energy-cap arrest (Pa)
#   R0            — G_c(τ_p0) / G_slab; < R_FIT (0.48) → crack propagates
#   A_ca_energy   — energy-cap arrest length (m)

meloche_df = compute_meloche_features(
    snap_data   = {"all": features_df},
    cluster_map = cluster_map,
    dem         = dem,
    transform   = transform,
)
print(f"  Meloche features: {len(meloche_df)} clusters with valid indices")

# ---------------------------------------------------------------------------
# 5. Save CSVs (generate_scenarios.py reads these directly)
# ---------------------------------------------------------------------------
OUT_DIR.mkdir(parents=True, exist_ok=True)
date_str = EVENT_DATE.strftime("%Y-%m-%d")

feat_out = OUT_DIR / f"all_start_zone_features_{date_str}.csv"
mel_out  = OUT_DIR / f"meloche_features_all_{date_str}.csv"

features_df.to_csv(feat_out)
meloche_df.to_csv(mel_out)
print(f"Wrote {feat_out}")
print(f"Wrote {mel_out}")

# ---------------------------------------------------------------------------
# 6. Quick sanity check — show the filter-chain candidate counts
# ---------------------------------------------------------------------------
# These should match the counts printed by generate_scenarios.py.
valid = meloche_df.dropna(subset=["tau_g", "A_ca_brittle"])
n_tau = (valid["tau_g"] >= 40).sum()
sk_col = "min_sk38" if "min_sk38" in features_df.columns else "sk38_min"
if sk_col in features_df.columns:
    n_sk = (features_df.loc[features_df.index.isin(valid.index), sk_col] < 1.0).sum()
    print(f"  tau_g >= 40 Pa: {n_tau} | Sk38 < 1.0: {n_sk}")
r0_valid = meloche_df["R0"].dropna()
print(f"  R0 median (release side should be < 0.48): {r0_valid.median():.3f}")
print(f"  A_ca_brittle median: {meloche_df['A_ca_brittle'].dropna().median():.1f} m")
