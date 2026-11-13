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
G_WL     = 0.2e6   # weak-layer shear modulus (Pa). Meloche et al. (2025) Table 1
                   # / p.7 attribute this to cold-lab shear tests by Reiweger et
                   # al. (2010); that attribution is not independently verified
                   # here (not in the arXiv:2406.01360v1 preprint's bibliography;
                   # the final JGR text was not checked) -- see methods §4, §6, §9.
NU       = 0.3     # slab Poisson's ratio
PHI_DEG  = 27.0    # snow friction angle (degrees)
DELTA    = 1.0     # Meloche softening coefficient δ. Matches the paper's own
                   # full-slope-scale campaign (Table 1); their Fig. 8 brittle
                   # runs, which calibrate C and R_FIT, used δ = 0.
D_WL_FALLBACK = 0.04   # m; weak-layer thickness used only when wl_thickness is NaN
L_SS     = 20.0    # m; steady-state propagation length, Meloche et al. (2025)

# -----------------------------------------------------------------------
# Feature extraction (compute_meloche_features)
# -----------------------------------------------------------------------
K_NEIGHBORS        = 6      # k-nearest cluster centroids used for the θ gradient
TAU_G_FEATURE_FLOOR = 50.0  # Pa; below this the arrest scaling laws are not
                            # evaluated. UNSUPPORTED: Meloche et al. never varied
                            # tau_g (h=0.5 m, rho=250, psi=35 fixed in all four
                            # campaigns, so tau_g = 703 Pa is the only simulated
                            # value). The paper's size caveat is "size 3 or more",
                            # not D2+. See THETA_VALID_* for the thresholds the
                            # paper does give.
THETA_MIN          = 1e-6   # Pa/m; θ below this is treated as no usable gradient

# How θ is estimated from the per-cluster τp field. See methods §5, "θ is a
# function of the lag it is measured at".
#   'knn'        mean |Δτp|/d over the K_NEIGHBORS nearest cluster centroids.
#                This is the shipped behaviour and what every committed CSV and
#                every published number was generated with.
#   'plane_fit'  local least-squares plane through τp over all centroids within
#                THETA_FIT_RADIUS_M, θ = |∇τp|. Estimates the slope-scale ramp
#                that Meloche's θ represents, instead of the local noise.
#
# Measured on Jan 18: start-zone clusters are 3.0 m across and the k = 6
# neighbourhood averages 3.3 m, so 'knn' samples θ inside the 0.5-10 m band the
# paper treats as local noise. The τp variogram gives θ(d) = 3.87 + 72/d Pa/m,
# i.e. at 3.3 m θ is 15% slope-scale trend and 85% local noise. A_ca ∝ 1/θ, so
# the two estimators differ by ~3.4x in A_ca.
#
# Default changed to 'plane_fit' on 2026-10-09, argued from the propagation
# scale the scaling law was calibrated at: Meloche's θ is the gradient of a
# linear ramp over L_ss, so the estimator's neighbourhood should be L_ss, not
# whatever the cluster spacing happens to be. 'knn' is retained to reproduce
# every pre-2026-10-09 CSV and published number.
#
# Jan 18 corroborates but did not decide it — see methods §6. Trigger-matched
# mean IoU is 0.606 against 0.598 for the old configuration, a margin far too
# thin on n = 1 to carry the choice. What Jan 18 does show decisively is the
# *interaction*: 'plane_fit' alone scores 0.536 and the corrected E/D_wl
# features alone 0.462, because the two move A_ca in opposite directions
# (~x2.9 and ~x0.5). Each correction looks wrong in isolation and they are
# right together, so single-factor IoU selection would have rejected all of
# them. The independent check is scale, not overlap: mean area/observed lands
# at 0.82 against the 0.83 a perfect mask-limited model implies, and only this
# configuration passes it (0.77 / 0.62 / 1.18 for the others).
THETA_ESTIMATOR          = 'plane_fit'   # 'knn' | 'plane_fit'
THETA_FIT_RADIUS_M       = L_SS    # m; neighbourhood radius for 'plane_fit'
THETA_FIT_MIN_NEIGHBOURS = 8       # minimum centroids in radius, else θ = NaN

# Validity range of the arrest scaling laws in θ, from Meloche et al. (2025).
# These are the only quantitative propagation thresholds the paper gives:
#   p.12 "Shear strength gradients of more than 20 Pa m-1 were needed to obtain
#         crack arrest within the PST simulated length."
#   p.13 runs with θ < 20 Pa/m were REMOVED from their analysis (no arrest).
#   p.10 "When the shear strength gradient is larger than around 5,000 Pa m-1,
#         the gradient is too large, causing crack arrest after L_ss."
# On Jan 18, 36% of start-zone clusters sit below 20 Pa/m, i.e. outside the
# calibrated regime. Not currently enforced — see methods §5.
THETA_VALID_MIN = 20.0      # Pa/m
THETA_VALID_MAX = 5000.0    # Pa/m

# -----------------------------------------------------------------------
# Layered-slab aggregation of sigma_t and E (arrest_indices.aggregate_slab)
# -----------------------------------------------------------------------
# Faceted snow is about half as strong in tension as rounded snow at equal
# density (Jamieson & Johnston 1990): sigma_t,i = sigma_t,RG(rho_i) (1 - a f_i).
FACET_STRENGTH_FACTOR = 0.5

# How f is derived. 'grain_class' (Ron, 2026-10-27) sets f = 1 for J&J's
# group II, "layers of faceted grains and of layers of grains which showed
# both faceting and rounding": primary class in FACETED_GRAIN_CLASSES
# (SNOWPACK 4 FC, 5 DH, 9 FCxr), or a MIXED_PRIMARY_CLASSES primary (PP/DF/RG,
# their group I) with a faceted secondary, e.g. 341; 0 otherwise. Melt forms
# stay 0 (J&J tested none). 'grain_class_primary' ignores the secondary.
# 'sphericity' is the earlier continuous f = 1 - sphericity, which
# under-counts FCxr: over the Jan 18 start-zone slabs it gave f 0.35 inside
# the release / 0.52 outside, against 0.50 / 1.00 by primary class. The Eq. 20 arrest reads the plain sigma_t
# column, which has no knockdown, so this moves sigma_t_mean/_wl/_ligament
# (and A_ca_* built on them), the energy mode and the slab tensile gate.
# 2026-11-13 (Ron, both repos): 'grain_class_primary' is the default. The mixed
# RG/FC forms keep RG strength: they are 51% of the top 0.45 m of the Jan 18
# release slab (pure RG 6%), so the secondary rule knocked the shallow slab
# down as if it were facets (S_t 0.44x a uniformly thinned slab; 1.00x with the
# primary rule). Start-zone median f 0.97 -> 0.79.
FACET_SOURCE = 'grain_class_primary'
FACETED_GRAIN_CLASSES = (4, 5, 9)
MIXED_PRIMARY_CLASSES = (1, 2, 3)

# ============================================================================
# WEAK-LAYER DETECTION (snowpack_features.split_wl_slab)
# ============================================================================
# Basal FCxr (9xx, rounding facets) under the FC/DH sequence counts as part of
# the basal weak layer (Ron, 2026-10-25). Without it, 251 of 732 Jan 18 crown
# clusters had no weak layer: 1-4 basal 990 elements (median 4.7 cm) sat under
# FC and the first-element test failed. Only a run of FCxr AT THE GROUND is
# admitted; FCxr above the first FC/DH element still ends the block, so every
# profile with an FC/DH base splits exactly as before.
BASAL_FCXR_IS_WL = True

# Which weak layer the release model uses (snowpack_features.select_wl_slab).
# 'basal' (default, operational): the basal FC/DH block above. It stays the
# operational choice because the forecast cares about the larger avalanches
# (Ron, 2026-11-12). 'threshold_sum' is an ADDITION, the layer of concern by
# the threshold-sum approach (Schweizer & Jamieson 2007, Cold Reg. Sci. Technol.
# 47, 50-59; Monti et al. 2014): the interface with the most "lemons". On the
# Jan 18 re-referenced replay it put the bed at a median 0.53 m against the
# observed 0.44 m (basal: 0.82 m). avachain CLI: analyze --wl-selection.
WL_SELECTION = 'basal'
WL_SELECTIONS = ('basal', 'threshold_sum')

# Threshold-sum lemons (Schweizer & Jamieson 2007, Table 3). One lemon each:
TSA_GS_DIFF_MM   = 0.75    # grain size difference across the interface >= (mm)
TSA_HH_DIFF      = 1.7     # hand-hardness index difference >= (1 F, 2 4F, 3 1F, 4 P, 5 K)
TSA_DEPTH_M      = (0.18, 0.94)  # weak-layer depth within (m)
TSA_PERSISTENT   = (4, 5, 6, 9)  # weak-layer primary grain class FC, DH, SH, FCxr
TSA_WL_GS_MM     = 1.25    # weak-layer grain size >= (mm)
TSA_WL_HH        = 1.3     # weak-layer hand-hardness index <=

# Extra sigma_t(rho) relations emitted by aggregate_slab as
# sigma_t_mean__<name> / sigma_t_wl__<name> (arrest_indices.ST_RELATIONS).
ST_RELATIONS_EXTRA = ('jamieson1990',)

# FACET_SOURCE = 'sphericity' only. None = f is 1 - sphericity outright. Our rounded-grain layers sit at
# sphericity ~0.86, so they carry f ~0.14 and lose ~7% of sigma_t; measured
# over the Jan 18 slabs this lowers sigma_t_mean to 0.85x the bulk value.
# Setting this to 0.86 instead anchors f so typical RG keeps the unmodified
# sigma_t,RG(rho) and only genuinely faceted layers are penalised.
FACET_SP_REF = None

# Slab elastic modulus relation. 'vanherwijnen2016' is E = 0.93 rho^2.8 Pa
# (van Herwijnen et al. 2016, J. Glaciol. 62(236), Eq. 8). 'project_fit' is the
# older in-house hand-fit 4 MPa (rho/300)^2.5, kept as an option.
E_RELATION         = 'vanherwijnen2016'
E_RELATIONS_EXTRA  = ('project_fit',)   # also emitted as E_eff__<name>

# Edge-crack (ligament) bound — HYPOTHESIS. Computed and emitted, but never
# fed into the BFS and never validated against an observed crown.
USE_LIGAMENT_BOUND = True
# SLAB_K_IC overrides the relation with a constant (Pa m^0.5); None = use
# K_IC_RELATION. 'schweizer2004' and 'kirchner2000' are implemented; only
# 'borstad2013' is still gated, its published equation never having been
# sourced and quoted back for checking.
SLAB_K_IC          = None
K_IC_RELATION      = 'schweizer2004'
LIGAMENT_MODEL     = 'elhaddad'   # or 'lefm' for the uncorrected form
SELF_CONSISTENT_A0 = False        # True evaluates a0 at F(a0/h) iteratively
DMAX_FACTOR        = 1.0          # d_max = DMAX_FACTOR * grain_size (both m)

# Extra ligament bounds computed in the same pass, as
# (label, K_Ic relation, dmax_factor) -> keys suffixed __<label>.
LIGAMENT_VARIANTS = (
    ('kirchner2000',     'kirchner2000',  1.0),
    ('schweizer2004_d1', 'schweizer2004', 1.0),
    ('schweizer2004_d2', 'schweizer2004', 2.0),
)

# ---- Arrest-index fits and slab-fracture choices (arrest_indices) ----------
# Here since 2026-11-11 (both repos; avachain model_params). Read at call time
# as config.X, never bound as a default argument.
C_FIT = 0.045  # Eq. 20 prefactor, fitted from Meloche et al. 2025 Fig. 8 (two runs); see docs/release_area_methods.md §3
R_FIT = 0.48   # energy ratio, fitted from Meloche et al. 2025 Fig. 8 (two runs); see docs/release_area_methods.md §3

# Floor on r_c in slab_tensile_support (m). SNOWPACK reports r_c = 0 where the
# slab's own weight overloads the weak layer (no cut needed): every such cluster
# on the 2025-26 Little Professor CSVs has tau_g > tau_p, and they are the deep
# slabs (97-100% of clusters with h > 3 m; 14 of 1608 on Jan 18, ~300 from late
# Feb). r_c -> 0 drives the bending load to 0 and S_t to infinity, which is the
# right limit (the slab cannot break before a zero-length crack), but it was a
# divide by zero (891 RuntimeWarnings per 11-date run). 1 mm keeps S_t finite
# (~1e5) without changing any gate outcome: the smallest positive r_c seen is
# ~1 mm and S_t there is far above 1. Ron, 2026-11-10; same in both repos.
RC_MIN_M = 1.0e-3

# The density ranges in K_IC_RELATIONS are a GUIDELINE, not a limit (Ron,
# 2026-11-09). They are what each relation was fitted over; Schweizer et al.
# (2004) is 80-300 kg/m^3 from Swiss beam tests. Little Professor is a Colorado
# high-alpine, wind-affected snowpack. On Jan 18 2026, 82% of the start-zone
# clusters (avachain) have a mean slab density above 300, and every one has at
# least one intact slab layer outside 80-300. Clamping evaluated those layers
# at the range edge and understated K_Ic of dense slabs: the median ligament
# K_Ic is 1318 clamped against 1650 Pa m^0.5 extrapolated (+25%). Kirchner et
# al. (2000), fitted to 540, has 3% of clusters out of range. So since
# 2026-11-10 (Ron) the power law is EXTRAPOLATED past its range by default
# (K_IC_RANGE_MODE); n_rho_clamped now counts the extrapolated layers. The
# trigger toughness gate (config.USE_TRIGGER_TOUGHNESS) reads the
# ligament K_Ic, so this moves which triggers can start a release.
K_IC_RANGE_MODE = 'extrapolate'   # 'extrapolate' | 'clamp' | 'nan'

# Trigger fracture toughness (Ron, 2026-11-09/10). A slab fracture at the
# trigger starts at one point and propagates, so the slab's tensile capacity
# there is also limited by its toughness: an edge crack of depth a in a slab of
# toughness K_Ic fails at K_Ic / (F sqrt(pi a)). F = 1.122 is the
# semi-infinite edge-crack factor and a = h/2 the crack depth; both are the
# values the 2026-11-09 sweep evaluated (ski non-triggers 13/46 -> 0/46, Jan 18
# triggers unchanged, robust to the K source). The finite-thickness factor
# (edge_crack_strength's polynomial) would be ~2.8 at a/h = 0.5, i.e. stricter;
# that was not evaluated.
TOUGHNESS_EDGE_F = 1.122
TOUGHNESS_A_FRAC = 0.5

# -----------------------------------------------------------------------
# Trigger selection filter chain (generate_scenarios)
# -----------------------------------------------------------------------
MIN_TAU_G          = 40.0   # Pa; minimum gravitational driving shear stress
MAX_SK38           = 1.0    # Sk38 must be below this to count as unstable
MIN_SLAB_THICKNESS = 0.3    # m; 0.5 until 2026-11-13 (Ron, both repos): Jan 18 failed at ~0.45 m
# m; the deepest slab explosive control is expected to trigger at the basal
# WL. This is an operational choice, not a physical result, and it binds hard.
# Sk38 at the WL falls almost monotonically with slab thickness (Spearman
# -0.93 in the avachain 2026-01-18 pool), so the Sk38 ranking picks the
# thickest slabs under this cap. In avachain, 2.0 -> 2.5 m replaced all five
# triggers on both 2026-01-18 and 2026-03-07, and 1.5 m left only two
# candidates on 2026-03-07. Use 1.5 for skier-triggered scenarios.
MAX_SLAB_THICKNESS = 2.0
SLOPE_MARGIN_DEG   = 2.0    # degrees added to stauchwall angle for the trigger gate
ELEVATION_PCTILE   = 50     # keep clusters at or above this elevation percentile

# -----------------------------------------------------------------------
# BFS crack propagation
# -----------------------------------------------------------------------
STAUCHWALL_DEG        = 29.0   # slope threshold for downslope arrest; Ron 2026-10-27, was 28.0
GAUME_ASPECT_CAP      = 2.5    # max cross-slope width / A_ca
MIN_POLYGON_AREA      = 200.0  # m²; discard smaller polygons
TAU_G_ABS_FLOOR       = 350.0  # Pa; absolute driving-stress floor
# Distance caps (used only while USE_DISTANCE_CAPS): d_up = max(A_ca,
# D_UP_MIN_M) (D_CAP_DEFAULT_M without an A_ca), d_down = max(Stauchwall
# distance, D_DOWN_MIN_M), d_lat = max(Gaume width, D_LAT_MIN_M). Were literals
# in propagate_crack until 2026-11-10 (same values as ~/avachain).
D_UP_MIN_M      = 20.0   # m
D_DOWN_MIN_M    = 20.0   # m
D_LAT_MIN_M     = 15.0   # m
D_CAP_DEFAULT_M = 50.0   # m
# Load-to-strength floor (Ron, 2026-11-10; same as ~/avachain; 0 = off):
# arrest a cluster whose tau_g / wl_shear_strength is below it. Dimensionless,
# unlike TAU_G_ABS_FLOOR. avachain sets it per slope via [calibration].
TAU_RATIO_FLOOR       = 0.0
MIN_PROPAGATION_SLAB  = 0.30   # m; thin-slab arrest threshold (0.50 until 2026-11-13, Ron: Jan 18 failed at ~0.45 m)
MIN_PROPAGATION_LAMBDA = 0.1   # m; compliant-slab arrest threshold
# Λ continuity thresholds, along-slope (mode II: upslope + downslope).
LAMBDA_DROP_FACTOR    = 0.20   # fractional softening arrest threshold
LAMBDA_RISE_FACTOR    = 0.30   # fractional stiffening arrest threshold

# Λ continuity thresholds, cross-slope (mode III / flanks). Separate knobs so
# the flanks can be tuned independently once Λ_III carries spatial structure
# (its formula is confirmed, but is a constant multiple of Λ_II under this
# pipeline's single global NU); identical to the along-slope values for now,
# so splitting the code changed no results.
LAMBDA_CROSS_DROP_FACTOR = 0.20
LAMBDA_CROSS_RISE_FACTOR = 0.30

# Cross-slope θ sampling in estimate_cross_slope_width(). The lateral sample
# has to sit near the true cross-slope axis (aspect ± 90°). Selecting it by
# raster row/col instead — "same row ± 50 px, column offset > 5 px" — silently
# assumes the fall line runs north-south. On the Jan 18 ESE start zone
# (aspect 113-120°, cross axis 23-30°) that put the chosen neighbours a mean
# 30° off the cross axis, with some of them within 1.2° of the fall line, so
# theta_cross was contaminated by the along-slope gradient.
THETA_CROSS_SECTOR_DEG = 30.0  # ± half-angle accepted about the cross axis
THETA_CROSS_MIN_SEP_M  = 5.0   # m; closer pairs make |Δτp|/d noise-dominated
THETA_CROSS_MAX_SEP_M  = 50.0  # m; beyond this the pair is not abeam

THICKNESS_DROP_FACTOR = 0.20   # fractional thinning arrest threshold
THICKNESS_RISE_FACTOR = 0.30   # fractional thickening arrest threshold

# ---- Mode III (cross-slope) elastic length -----------------------------
# The cross-slope elastic length is a different and SMALLER quantity than the
# along-slope (mode II) Λ. `arrest_indices.elastic_length_cross()` now matches
# Johan Gaume's own antiplane-shear derivation exactly (personal
# communication, 2026-10-18; not yet published in the literature) -- but
# whether antiplane shear is even the right cross-slope regime, vs.
# bending-dominated deformation, is explicitly still open on Gaume's own
# account (see docs/release_area_methods.md §3). Propagation currently uses
# the mode II Λ in every direction regardless, and the cross-slope branch is
# a no-op, for the reason in the IMPORTANT note below.
#
# Before flipping USE_MODE3_LAMBDA to True:
#   1. confirm the regime question above has moved (antiplane vs. bending);
#   2. make Λ_III vary spatially -- a constant ratio changes nothing
#      regardless of which formula turns out right, see below.
# Nothing else should need to change — if it does, the directional split in
# release_geometry._qualifies() has drifted and needs fixing instead.
#
# The `Lambda_cross` column is already written to the meloche CSV from the
# confirmed derivation, so the data side is ready; this flag only controls
# whether the BFS reads it. Note the shipped Jan 18 CSV predates the column —
# regenerate it before the flag can do anything.
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
#
# DEFAULT FLIPPED True -> False on 2026-10-10 (Ron's call). The argument is
# COMPLETENESS, not fit: only the *restrictive* half of mode III is coded. The
# speed cap shortens lateral propagation, but because G < E' an equal-strength
# slab also has a LARGER energy cap cross-slope — G_slab,III = tau_flank^2 h/(2G)
# against sigma_t^2 h/(2E') — which pushes the other way. That term needs a slab
# flank strength tau_flank for which no parameterisation exists (blocked on a
# Cam-Clay beta ambiguity, methods §3), so it is absent. Applying one side of a
# two-sided effect biases the model toward over-arrest, and the honest default
# for a half-implemented term is off.
#
# Jan 18 is consistent with that reading but did NOT decide it: cap off gives
# mean IoU 0.622 against 0.562 and area/observed 0.78 against 0.63, i.e. the
# cap-on model under-predicts area and the cap is what makes it do so. Per the
# no-tuning rule, read that as corroboration of the one-sided-implementation
# argument, not as the justification — on n = 1 the metric cannot carry a
# decision (see the §6 2x2, where single-factor IoU selection would have
# rejected all three corrections that are jointly the best configuration).
#
# TO FLIP IT BACK: once tau_flank lands, turn BOTH halves on together —
# slab_energy_cap_cross plus this cap. Turning this one back on alone would
# restore the same one-sided bias.
USE_MODE3_SPEED_CAP = False  # apply the ratio to the lateral distance cap

# ---- Anisotropic BFS rings (cross-slope growth rate) -------------------
# The same speed ratio, applied to the *growth* of the flood-fill rather than
# only to a terminal distance cap. If one BFS ring stands for one increment of
# crack-front advance in TIME, then a front running at 1.6 c_s along-slope and
# c_s cross-slope should already be ~1.6x further along-slope when it arrives.
# The isotropic ring loop instead advances one cluster hop in every direction
# per ring and only clips afterwards with d_up / d_down / d_lat, which leaves
# the front shape an artifact of cluster geometry.
#
# With rings on, a cross-slope candidate is admitted at ring k only once
#     dist <= k * RING_STEP_M * (MODE3_SPEED_RATIO / MODE2_SPEED_RATIO)
# and is otherwise DEFERRED to a later ring — arrest is still evaluated when
# the front arrives, never skipped because it arrived early. The event has a
# finite duration: when nothing else can be admitted, clusters still outside
# the cross-slope budget are dropped, so the flanks are cut at ~0.625x the
# along-slope reach.
#
# RING_STEP_M = None derives the step from the median nearest-neighbour cluster
# centroid separation, i.e. the natural along-slope advance of one hop, so the
# fast axis grows at the unmodified BFS rate and only the slow axis is held
# back. Set a number to override.
#
# Judge this against the observed crown's aspect (Jan 18 major axis 119 deg,
# L/W 3.51), which is a sharper test than IoU, and do NOT tune the ratio — 1.6
# is Broberg's, not ours. Note it reuses the SAME 1.6 as USE_MODE3_SPEED_CAP,
# so running both risks double-counting one physical effect; see methods §3.
USE_ANISOTROPIC_RINGS = False
RING_STEP_M           = None
BFS_K_NEIGHBOURS      = 8      # cluster adjacency degree for the flood-fill
MAX_BFS_CLUSTERS      = 2000   # safety cap on region size; see methods §5

# Meloche per-direction arrest test inside the BFS. ON since 2026-10-23 (Ron),
# matching ~/avachain. When on, the flood-fill stops at the directional Eq. 20
# test (arrest_indices.arrest_length, C = 0.045) and SKIPS the Λ/thickness
# discontinuity heuristics; TAU_G_ABS_FLOOR, thin slab, low Λ and the distance
# caps still apply. See methods §5.
USE_MELOCHE_ARREST = True

# Slab tensile-failure gate (Ron, 2026-10-26; ON by default since 2026-10-27,
# Ron, matching ~/avachain; --no-slab-tensile-gate turns it off). Arrest a
# candidate whose slab would break in tension before its weak-layer crack
# reaches critical length: S_t = sigma_t / (3 rho g r_c^2 / h) < 1, Schweizer et
# al. (2014), via arrest_indices.slab_tensile_support. Per-cluster, after the
# tau_g / thin-slab / low-Lambda gates. sigma_t is Jamieson & Johnston (1990)
# with the facet knockdown (f by grain class since 2026-10-27, FACET_SOURCE),
# the strength relation the 2014 paper uses. The threshold is 1, not a fit. Same gate as ~/avachain release_geometry.
# Why on: S_t beat every random-arrest control of the same count, and with
# grain-class sigma_t it costs no trigger it lets through. Against: where S_t = 1
# falls depends on how r_c is scaled. A CSV without SLAB_TENSILE_SIGMA_COL turns
# it into a no-op; generate_scenarios warns when that happens.
USE_SLAB_TENSILE_GATE  = True
SLAB_TENSILE_SIGMA_COL = 'sigma_t_mean__jamieson1990'

# Trigger fracture toughness (Ron, 2026-11-10; same as ~/avachain
# release_geometry; CLI --no-trigger-toughness). At the TRIGGER only, sigma_t in
# S_t is capped by the slab's edge-crack strength K_Ic / (1.122 sqrt(pi h/2))
# (arrest_indices.toughness_limited_strength), K_Ic the ligament column
# (Schweizer 2004, extrapolated past its density range). The BFS per-cluster
# tensile gate is unchanged. On 2026-11-09 (avachain) it took the Jan 18
# ski-track non-triggers from 13/46 releasing to 0/46 and left the five Jan 18
# triggers and IoU unchanged; it does not cut the no-event false alarms. A CSV
# without TOUGHNESS_K_COL leaves the plain S_t.
USE_TRIGGER_TOUGHNESS = True
TOUGHNESS_K_COL       = 'K_Ic_ligament'

# BFS distance caps (up: A_ca, down: Stauchwall distance, cross: Gaume width).
# Default on, both repos. Off (CLI --no-distance-caps; avachain: a slope's
# [calibration] distance_caps = false) leaves the release bounded by the
# per-cluster gates, the start zone and MAX_BFS_CLUSTERS, and the size factors
# then scale nothing in the BFS.
USE_DISTANCE_CAPS = True

# Running slab cracks (Ron, 2026-11-10; same as ~/avachain release_geometry;
# CLI --no-slab-crack-lines).
# A slab that the per-cluster tensile gate breaks (slab_tensile_failure) starts
# a crack IN THE SLAB: a crown if the BFS reached it by an upslope step (runs
# cross-slope), a flank if by a cross-slope step (runs along the fall line);
# downslope steps start none (Stauchwall, compression). After the BFS each crack
# runs as a line, one best-aligned neighbour per step (within
# SLAB_CRACK_ALIGN_COS of the line), while the slab there is loaded: inside the
# admitted region, with detached slab below a crown / beside a flank within
# SLAB_CRACK_BAND_M, and the band-mean slope above PHI_DEG. The area the
# lines sever from the trigger (4-connected pixels) is removed; the line
# clusters stay (they are the boundary). Order-independent by construction.
# No fracture toughness in the run test: with any static load of a detached
# slab (full weight, or Stauchwall- and flank-supported, q*min(L,W)/2), K was
# 10-75x the slab K_Ic on Jan 18 (0 of 2220 steps stopped), because a crack a
# few metres long needs < ~0.5 kPa to arrest. Jan 18 in avachain (calibration 2026-11-10.2):
# IoU 0.601 -> 0.610, area 7429 -> 6831 m2, AUC 0.801 -> 0.848, ski 0/46, same
# five triggers; 419 m2 overshoot removed for 179 m2 of release; 9/10 no-event
# days unchanged (the tensile gate never fires there).
USE_SLAB_CRACK_LINES = True
SLAB_CRACK_BAND_M    = 5.0     # m; half-width of the strip that sets L and W
SLAB_CRACK_ALIGN_COS = 0.707   # a line step must lie within 45 deg of the line

# Severed slab slides (Ron, 2026-11-13; both repos, same as ~/avachain). The
# crack lines above drop every piece they cut off from the trigger. But the
# weak layer under that piece has already fractured (the BFS reached it), so
# once the slab breaks around it nothing holds it but crack-face friction. With
# this on, a severed piece (connected in the BFS graph) stays in the release
# when its mean cluster slope exceeds atan(SEVERED_SLIDE_MU); only a piece on
# gentler ground is dropped. 0.62 is the crack-face friction where slabs came
# to rest, van Herwijnen & Heierli 2009 GRL (31.8 deg). On a 0.45 m Jan 18 slab
# one line across the neck at the trigger otherwise dropped ~410 of ~650
# reached clusters at facet_a 0.5 and none at 0.6 (a knife edge).
USE_SEVERED_SLIDE = True
SEVERED_SLIDE_MU  = 0.62

# -----------------------------------------------------------------------
# Path-energy arrest (energy_path.py) — hypothesis, see methods §3
# -----------------------------------------------------------------------
# Energy-mode delta in G_c = (tau_p - tau_r)^2 (1 + delta) / (2 K_wl), separate
# from the Eq. 20 DELTA (1.0): R = 0.48 was fitted to the JGR runs at delta = 0,
# so the R-based test uses 0 (Ron, 2026-11-11; avachain ENERGY_DELTA, same value).
# It was config.DELTA = 1.0 here until 2026-11-11. CLI --energy-delta.
ENERGY_DELTA        = 0.0
ARREST_MODE         = 'legacy'        # BFS arrest rule: 'legacy' (propagate_crack as before) or 'energy' (path energy)
SIGMA_T_COLUMN      = 'sigma_t_mean'  # thickness-weighted slab σt; 'sigma_t_wl' = weakest-link bound (methods §4)
ENERGY_USE_RESIDUAL = False           # τr = 0 in G_c, as arrest_indices.evaluate(use_residual=False)
MASK_MODE           = 'gate'          # start-zone KML: 'gate' (BFS rejects, polygon clipped) or 'diagnostic' (reporting only); brief Step 4
ENERGY_R            = None           # energy ratio R; None = R_FIT (Meloche Fig. 8 fit). Set only for sensitivity runs

# -----------------------------------------------------------------------
# Scenario sweep
# -----------------------------------------------------------------------
# Number of trigger clusters to evaluate. Dropped 5 → 3 on 2026-10-09: five
# overlapping polygons made the comparison figure unreadable, and the two
# dropped scenarios were the two worst-ranked by Sk38. Note the cut now lands
# inside a three-way Sk38 tie at 0.72 (clusters 2859, 5656, 6191 on Jan 18), so
# which three you get is decided by the ordering of a tie, not by a margin —
# see methods §8 item 5. Override with --n-triggers.
N_TOP_TRIGGERS   = 3
SIZE_FACTORS     = [1.0]   # release size multipliers; one polygon per trigger
