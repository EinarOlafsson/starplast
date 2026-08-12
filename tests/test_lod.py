"""The level-of-detail tiers, and the bug that made the coarsest one meaningless.

The galaxy tier drew one centroid per compartment and produced a knot of dots in the middle of the
screen. It was not a rendering fault: averaging points that are spread across the whole map lands you
near the middle of it. These tests pin the property that matters -- a centroid has to be closer to its
own members than to everything else -- which the compartment version failed and the spatial version
passes.
"""
import numpy as np
import pytest

from starplast import lod


def _blobs(centres, n=60, sd=0.25, seed=0):
    rng = np.random.default_rng(seed)
    return np.vstack([rng.normal(c, sd, size=(n, 3)) for c in centres])


def test_separated_blobs_become_separate_galaxies():
    xyz = _blobs([(0, 0, 0), (10, 0, 0), (0, 10, 0)], n=80)
    lab = lod.galaxies(xyz, min_members=20)
    assert len(set(lab[lab >= 0])) == 3
    # Each input blob must land wholly in one component, or the tier is splitting real structures.
    for k in range(3):
        block = lab[k * 80:(k + 1) * 80]
        assert len(set(block[block >= 0])) == 1


def test_galaxies_are_numbered_largest_first():
    """A landmark that changes number between runs is not a landmark."""
    xyz = np.vstack([_blobs([(0, 0, 0)], n=200), _blobs([(20, 0, 0)], n=60, seed=1)])
    lab = lod.galaxies(xyz, min_members=20)
    sizes = [int((lab == i).sum()) for i in sorted(set(lab[lab >= 0]))]
    assert sizes == sorted(sizes, reverse=True)


def test_a_centroid_represents_its_own_members():
    """The property the compartment version violated.

    A component's centroid must sit closer to its members than the members are spread, or the point
    drawn for it is not standing for anything.
    """
    xyz = _blobs([(0, 0, 0), (12, 0, 0), (0, 12, 0)], n=80)
    lab = lod.galaxies(xyz, min_members=20)
    pos, num, spread = lod.centroids(xyz, lab)
    span = float(np.linalg.norm(xyz - xyz.mean(0), axis=1).max())
    for p, s in zip(pos, spread):
        assert s < span * 0.5, "members are more spread out than the map, so the mean means nothing"
    # And the centroids must be far apart -- the failure mode was all of them landing together.
    d = [np.linalg.norm(a - b) for i, a in enumerate(pos) for b in pos[i + 1:]]
    assert min(d) > max(spread), "centroids are closer to each other than a component is wide"


def test_a_scattered_label_would_have_collapsed_to_the_centre():
    """Demonstrates the original bug directly, so the reason for the change cannot be lost.

    Three interleaved labels spread over the same cloud have almost identical means, near the middle.
    That is what compartment centroids were.
    """
    xyz = _blobs([(0, 0, 0), (12, 0, 0), (0, 12, 0)], n=90)
    scattered = np.arange(len(xyz)) % 3          # a label with no spatial meaning
    centre = xyz.mean(0)
    span = float(np.linalg.norm(xyz - centre, axis=1).max())
    means = np.stack([xyz[scattered == k].mean(0) for k in range(3)])
    assert all(np.linalg.norm(m - centre) < span * 0.2 for m in means)
    # Whereas the spatial tier keeps them apart.
    pos, _, _ = lod.centroids(xyz, lod.galaxies(xyz, min_members=20))
    assert max(np.linalg.norm(p - centre) for p in pos) > span * 0.3


def test_empty_input_is_not_an_error():
    assert lod.galaxies(np.zeros((0, 3))).shape == (0,)
    pos, num, spread = lod.centroids(np.zeros((0, 3)), np.zeros(0, dtype=int))
    assert len(pos) == len(num) == len(spread) == 0


def test_all_noise_returns_no_galaxies():
    """Points too sparse to fill any cell must yield nothing rather than one meaningless component."""
    rng = np.random.default_rng(3)
    xyz = rng.normal(0, 40, size=(30, 3))
    lab = lod.galaxies(xyz, min_cell=10, min_members=500)
    assert (lab == -1).all()
    pos, _, _ = lod.centroids(xyz, lab)
    assert len(pos) == 0


def test_a_flat_embedding_does_not_divide_by_zero():
    """A degenerate axis is a real possibility for a 2D embedding shown in a 3D view."""
    xyz = np.column_stack([np.linspace(0, 10, 200), np.linspace(0, 10, 200), np.zeros(200)])
    lab = lod.galaxies(xyz, min_members=10)
    assert np.isfinite(xyz).all() and len(lab) == 200


def test_dominant_reports_a_fraction_not_just_a_name():
    """"cytosol" is a claim; "cytosol, 21%" is a description, and these components are not pure."""
    xyz = _blobs([(0, 0, 0), (12, 0, 0)], n=60)
    lab = lod.galaxies(xyz, grid=4, min_members=20)
    vals = np.array(["a"] * 60 + ["b"] * 60)
    dom = lod.dominant(vals, lab)
    for name, frac in dom.values():
        assert 0.0 < frac <= 1.0
        assert name in ("a", "b")


def test_dominant_ignores_absence_labels():
    """Without this every galaxy in the real map is named "unassigned" -- true, and useless."""
    # grid=3 because a single blob normalised to its own extent spreads across a fine grid until no
    # cell is occupied. The tier is relative to the spread of what it is given.
    xyz = _blobs([(0, 0, 0)], n=90)
    lab = lod.galaxies(xyz, grid=3, min_members=20)
    vals = np.array(["unassigned"] * 60 + ["cytosol"] * 30)
    assert lod.dominant(vals, lab)[0][0] == "unassigned"
    name, frac = lod.dominant(vals, lab, ignore={"unassigned"})[0]
    assert name == "cytosol"
    # The fraction is of NAMED members, not of the whole component.
    assert frac == pytest.approx(1.0)


def test_dominant_skips_a_component_with_nothing_named():
    xyz = _blobs([(0, 0, 0)], n=60)
    lab = lod.galaxies(xyz, grid=3, min_members=20)
    vals = np.array(["unassigned"] * 60)
    assert lod.dominant(vals, lab, ignore={"unassigned"}) == {}


def test_min_frac_refuses_to_name_a_component_with_no_majority():
    xyz = _blobs([(0, 0, 0)], n=90)
    lab = lod.galaxies(xyz, grid=3, min_members=20)
    vals = np.array(["a"] * 30 + ["b"] * 30 + ["c"] * 30)
    assert lod.dominant(vals, lab, min_frac=0.5) == {}
    assert lod.dominant(vals, lab, min_frac=0.2) != {}


def test_corner_touching_blobs_stay_separate():
    """6-connectivity is deliberate: two lobes meeting at a corner share one point, not a structure."""
    a = _blobs([(0.0, 0.0, 0.0)], n=120, sd=0.15)
    b = _blobs([(3.0, 3.0, 3.0)], n=120, sd=0.15, seed=7)
    lab = lod.galaxies(np.vstack([a, b]), grid=6, min_members=20)
    assert len(set(lab[lab >= 0])) == 2
