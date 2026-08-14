#!/usr/bin/env python3
"""Light that can be blocked.

The claim these tests defend is a narrow one, and the module says so in its own docstring: this is
volume sampling of a coarse density grid, not ray tracing. So they ask what that claim implies -- a
clear line delivers light, a cluster in the way takes it away, and a ray fired into the map lands
somewhere -- and they do NOT ask for exactness the model does not have.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import rays as R  # noqa: E402


def wall(n=600, seed=0):
    """A slab at x=0 as dense as a UMAP cluster, with a lone gene either side of it.

    As dense deliberately: a shadow here is transmittance through the binned density, so a diffuse
    occluder is SUPPOSED to leave most of the light through. Testing against a thin one measures
    how thin it is, not whether shadows work.
    """
    rng = np.random.default_rng(seed)
    slab = np.column_stack([rng.normal(0, 0.4, n), rng.uniform(-3, 3, n), rng.uniform(-3, 3, n)])
    return np.vstack([slab, [[-20.0, 0.0, 0.0], [20.0, 0.0, 0.0]]])


# --------------------------------------------------------------------------- the grid
def test_the_grid_holds_the_cloud_and_nothing_else():
    coords = wall()
    g = R.build(coords)
    assert g.resolution == R.RESOLUTION
    assert g.occupancy.max() == pytest.approx(1.0), "no cell reached opaque in a dense slab"
    assert float((g.occupancy > 0).mean()) < 0.2, "a thin slab filled the whole box"
    i, j, k = g.cells_at(coords)
    assert g.occupancy[i, j, k].min() > 0, "a gene landed in a cell the grid calls empty"


def test_an_empty_map_bins_to_an_empty_box_rather_than_raising():
    g = R.build(np.zeros((0, 3)))
    assert g.occupancy.shape == (1, 1, 1)
    assert g.transmittance(np.zeros((0, 3)), np.ones(3)).shape == (0, 1)


def test_points_outside_the_box_are_clamped_to_it():
    g = R.build(wall())
    i, j, k = g.cells_at(np.array([[1e6, -1e6, 0.0]]))
    assert (0 <= i).all() and (i < g.resolution).all()
    assert (0 <= j).all() and (j < g.resolution).all()


# --------------------------------------------------------------------------- shadows
def test_a_cluster_in_the_way_takes_the_light_away():
    """The whole feature: a gene lights up when the line between it and the light is clear."""
    coords = wall()
    g = R.build(coords)
    lamp = np.array([-40.0, 0.0, 0.0])
    seen = g.transmittance(coords, lamp)
    near, far = seen[-2, 0], seen[-1, 0]     # the lone gene on the lamp's side, and the one behind
    assert near > 0.9, f"a clear line only delivered {near:.2f}"
    assert far < 0.5, f"the slab let {far:.2f} through"


def test_a_gene_does_not_shadow_itself_and_the_lamp_does_not_shadow_everything():
    """Sampled strictly between the ends. A sample at the gene finds the gene, and every gene goes
    dark; a sample at the lamp finds whatever the lamp stands in, and the whole map goes dark."""
    coords = wall()
    g = R.build(coords)
    inside = coords[5]                       # a lamp standing in the middle of the slab
    seen = g.transmittance(coords, inside)
    assert seen.max() > 0.9, "a lamp inside a cluster shadowed the entire map"
    lit = g.transmittance(np.array([[-20.0, 0.0, 0.0]]), np.array([-40.0, 0.0, 0.0]))
    assert lit[0, 0] > 0.95, "a lone gene shadowed itself"


def test_more_of_the_map_in_the_way_means_less_light():
    """Soft rather than binary, which is what the density model buys and what it costs: there is no
    hard shadow edge anywhere, and one lonely gene in the way barely dims anything."""
    thin = np.vstack([np.zeros((2, 3)), [[20.0, 0.0, 0.0]]])
    dense = np.vstack([np.random.default_rng(1).normal(0, 0.3, (400, 3)), [[20.0, 0.0, 0.0]]])
    lamp = np.array([-40.0, 0.0, 0.0])
    behind = np.array([[20.0, 0.0, 0.0]])
    assert R.build(thin).transmittance(behind, lamp)[0, 0] > \
        R.build(dense).transmittance(behind, lamp)[0, 0]


# --------------------------------------------------------------------------- bounces
def test_a_ray_fired_into_the_map_lands_on_it():
    coords = wall()
    g = R.build(coords)
    hits = g.first_hit(np.array([-40.0, 1.0, 0.0]), np.array([[1.0, 0.0, 0.0]]), reach=80.0)
    assert len(hits) == 1
    assert abs(hits[0][0]) < 2.0, f"landed at x={hits[0][0]:.1f} rather than on the slab"


def test_a_ray_that_meets_nothing_bounces_off_nothing():
    """No hit rather than a hit at the end of the ray: a bounce invented in empty space is a light
    the reader cannot account for."""
    g = R.build(wall())
    assert g.first_hit(np.array([-40.0, 40.0, 0.0]), np.array([[0.0, 0.0, -1.0]]), reach=80.0) == []


def test_a_ray_does_not_bounce_off_the_cell_it_started_in():
    """A torch stands INSIDE the cloud. Marched from zero, every ray bounces where it started and
    the bounce is a copy of the source sitting on top of it."""
    coords = wall()
    g = R.build(coords)
    origin = coords[0]                       # in the thick of the slab
    hits = g.first_hit(origin, np.array([[1.0, 0.0, 0.0]]), reach=80.0)
    assert not hits or np.linalg.norm(hits[0] - origin) > g.cell


def test_a_fan_of_rays_is_the_same_fan_every_frame():
    """Random directions are a different set each frame, and the bounce flickers between whatever
    each one happened to land on."""
    a = R.cone((0.0, 0.0, 1.0), 7)
    assert len(a) == 7
    assert np.array_equal(a, R.cone((0.0, 0.0, 1.0), 7))
    assert a[0] == pytest.approx([0.0, 0.0, 1.0]), "the first ray is not the one aimed"
    spread = np.linalg.norm(a[1:] - a[0], axis=1)
    assert spread.min() > 0 and spread.max() < 1.0


def test_a_fan_aimed_along_x_still_has_two_axes_across_it():
    """The obvious way to find an axis across a ray is to cross it with x, which is zero exactly
    when the ray IS x -- and "straight along x" is what a default camera looks down."""
    a = R.cone((1.0, 0.0, 0.0), 5)
    assert np.isfinite(a).all()
    assert len({tuple(np.round(v, 6)) for v in a}) == 5, "the rays collapsed onto each other"


def test_a_fan_with_no_direction_points_at_the_viewer():
    assert R.cone((0.0, 0.0, 0.0), 3)[0] == pytest.approx([0.0, 0.0, 1.0])


def test_a_fan_of_one_is_just_the_ray():
    assert len(R.cone((0.0, 1.0, 0.0), 1)) == 1
    assert len(R.cone((0.0, 1.0, 0.0), 0)) == 1
