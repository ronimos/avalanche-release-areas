"""
release_areas — Avalanche release zone delineation from SNOWPACK features.

Modules
-------
arrest_indices   Meloche et al. (2025) crack-arrest scaling laws
features         Per-cluster SNOWPACK feature extraction
release_geometry BFS crack-propagation and polygon construction
scenario_writer  Scenario-to-GeoJSON / CSV output
config           Default physical and operational parameters
"""

from release_areas.arrest_indices import evaluate, R_FIT
from release_areas.release_geometry import (
    make_release_polygon_2d,
    propagate_release,
    fill_polygon_holes,
    plot_release_comparison,
)
from release_areas.features import (
    geojson_to_mask,
    compute_meloche_features,
)

__all__ = [
    "evaluate",
    "R_FIT",
    "make_release_polygon_2d",
    "propagate_release",
    "fill_polygon_holes",
    "plot_release_comparison",
    "geojson_to_mask",
    "compute_meloche_features",
]
