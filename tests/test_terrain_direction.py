"""Tests for the terrain direction helpers in release_geometry.

Two invariants matter here, and both are about *direction*, not magnitude:

  1. compute_slope_aspect() must return the true compass bearing of steepest
     descent on a north-up raster, where row increases southward.
  2. find_stauchwall() must walk along that bearing. An integer
     dr=sign(-cos), dc=sign(sin) step reaches only 45/135/225/315 deg, so on
     an ESE slope it veers ~20 deg clockwise toward south and the release
     polygon rotates with it.
"""
import numpy as np
import pytest
from affine import Affine

from release_areas.release_geometry import (
    compute_slope_aspect,
    find_stauchwall,
    pixel_to_utm,
    project_along_aspect,
)

# North-up 1 m transform: col -> +easting, row -> -northing.
T = Affine(1.0, 0.0, 400000.0, 0.0, -1.0, 4000000.0)


def planar_dem(bearing_deg, slope_deg, n=201):
    """DEM of a perfect plane whose steepest descent runs along bearing_deg."""
    r, c = np.mgrid[0:n, 0:n]
    x = c * 1.0                      # easting, +col
    y = -r * 1.0                     # northing, -row
    b = np.radians(bearing_deg)
    # unit downslope direction in map frame
    e_hat, n_hat = np.sin(b), np.cos(b)
    # descend along (e_hat, n_hat)
    return -np.tan(np.radians(slope_deg)) * (x * e_hat + y * n_hat)


def ang_diff(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


class TestComputeSlopeAspect:
    @pytest.mark.parametrize("bearing", [0, 45, 90, 114, 135, 180, 225, 270, 315])
    def test_recovers_plane_bearing(self, bearing):
        dem = planar_dem(bearing, 30.0)
        slope, aspect = compute_slope_aspect(dem, 1.0)
        mid = (slice(5, -5), slice(5, -5))
        assert abs(ang_diff(float(np.mean(aspect[mid])), bearing)) < 1e-6
        assert float(np.mean(slope[mid])) == pytest.approx(30.0, abs=1e-6)

    def test_east_facing_slope_is_90_deg(self):
        # elevation falls toward +easting only
        r, c = np.mgrid[0:50, 0:50]
        _, aspect = compute_slope_aspect(-0.5 * c.astype(float), 1.0)
        assert float(aspect[10, 10]) == pytest.approx(90.0)

    def test_south_facing_slope_is_180_deg(self):
        # elevation falls toward -northing, i.e. toward +row
        r, c = np.mgrid[0:50, 0:50]
        _, aspect = compute_slope_aspect(-0.5 * r.astype(float), 1.0)
        assert float(aspect[10, 10]) == pytest.approx(180.0)


class TestFindStauchwall:
    @pytest.mark.parametrize("bearing", [114.0, 100.0, 120.0, 160.0, 20.0, 290.0])
    def test_walk_follows_aspect_not_the_diagonals(self, bearing):
        """The realised trigger->stauchwall bearing must match the aspect.

        The old sign()-quantised step could only produce 45/135/225/315, so
        this fails by up to 45 deg on any non-diagonal slope.
        """
        n = 201
        dem = planar_dem(bearing, 35.0, n)
        slope, aspect = compute_slope_aspect(dem, 1.0)
        # uniform 35 deg plane: never drops below threshold, so the walk runs
        # until it leaves the grid -> pure direction test
        r0 = c0 = n // 2
        r1, c1 = find_stauchwall(r0, c0, slope, aspect,
                                 threshold_deg=10.0, max_steps=60)
        x0, y0 = pixel_to_utm(r0, c0, T)
        x1, y1 = pixel_to_utm(r1, c1, T)
        realised = np.degrees(np.arctan2(x1 - x0, y1 - y0)) % 360.0
        assert abs(ang_diff(realised, bearing)) < 2.0

    def test_ese_slope_does_not_veer_south(self):
        """Regression: the Little Professor case. Aspect 114 deg (ESE) used to
        walk at exactly 135 deg (SE), a +21 deg clockwise error."""
        n = 201
        dem = planar_dem(114.0, 35.0, n)
        slope, aspect = compute_slope_aspect(dem, 1.0)
        r0 = c0 = n // 2
        r1, c1 = find_stauchwall(r0, c0, slope, aspect,
                                 threshold_deg=10.0, max_steps=60)
        x0, y0 = pixel_to_utm(r0, c0, T)
        x1, y1 = pixel_to_utm(r1, c1, T)
        realised = np.degrees(np.arctan2(x1 - x0, y1 - y0)) % 360.0
        assert abs(ang_diff(realised, 135.0)) > 15.0, "still snapped to the SE diagonal"
        assert abs(ang_diff(realised, 114.0)) < 2.0

    def test_stops_at_slope_threshold(self):
        slope = np.full((50, 50), 40.0)
        slope[:, 30:] = 10.0              # flat beyond col 30
        aspect = np.full((50, 50), 90.0)  # due east
        r, c = find_stauchwall(10, 10, slope, aspect, threshold_deg=28.0)
        assert (r, c) == (10, 30)

    def test_returns_immediately_when_already_below_threshold(self):
        slope = np.full((20, 20), 5.0)
        aspect = np.full((20, 20), 135.0)
        assert find_stauchwall(7, 8, slope, aspect, threshold_deg=28.0) == (7, 8)

    def test_stops_at_grid_edge_without_error(self):
        slope = np.full((20, 20), 40.0)
        aspect = np.full((20, 20), 90.0)   # east, will run off the right edge
        r, c = find_stauchwall(10, 10, slope, aspect, threshold_deg=28.0)
        assert 0 <= r < 20 and 0 <= c < 20

    def test_stops_at_nodata(self):
        slope = np.full((20, 20), 40.0)
        slope[:, 15:] = np.nan
        aspect = np.full((20, 20), 90.0)
        r, c = find_stauchwall(10, 10, slope, aspect, threshold_deg=28.0)
        assert (r, c) == (10, 14)


class TestProjectAlongAspect:
    def test_downslope_matches_aspect_bearing(self):
        x, y = project_along_aspect(0.0, 0.0, 100.0, 90.0, upslope=False)
        assert (x, y) == pytest.approx((100.0, 0.0), abs=1e-9)

    def test_upslope_is_opposite(self):
        x, y = project_along_aspect(0.0, 0.0, 100.0, 90.0, upslope=True)
        assert (x, y) == pytest.approx((-100.0, 0.0), abs=1e-9)

    def test_round_trip_on_a_plane_loses_elevation_downslope(self):
        dem = planar_dem(114.0, 30.0)
        _, aspect = compute_slope_aspect(dem, 1.0)
        r0 = c0 = 100
        a = float(aspect[r0, c0])
        # step 20 m downslope in pixel space and confirm the DEM drops
        dr = -np.cos(np.radians(a)) * 20.0
        dc = np.sin(np.radians(a)) * 20.0
        assert dem[int(r0 + dr), int(c0 + dc)] < dem[r0, c0]
