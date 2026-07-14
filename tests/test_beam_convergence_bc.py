# -*- coding: utf-8 -*-
"""Boundary-condition index used by the beam_convergence example."""

import numpy as np


def right_edge_mid_height_index(nx, ny):
    """Return the node index at the right edge, mid-height of a rect_grid.

    Nodes are ordered with x varying fastest over a (nx+1) by (ny+1) lattice
    (see solidspy.preprocesor.rect_grid / beam_convergence example).
    """
    return (nx + 1) * (ny // 2 + 1) - 1


def test_right_edge_mid_height_index_matches_coordinates():
    """The y-constraint must land on x = L/2, y = 0 for even ny.

    The previous formula ``nx * (ny//2 + 1) - 1`` selects an interior node.
    """
    L, h = 8.0, 1.0
    for nx, ny in [(3, 1), (6, 2), (12, 4), (24, 8)]:
        n_nodes = (nx + 1) * (ny + 1)
        # Coordinates as in rect_grid: x from -L/2 to L/2, y from -h/2 to h/2
        x = np.linspace(-L / 2, L / 2, nx + 1)
        y = np.linspace(-h / 2, h / 2, ny + 1)
        # Flatten with x varying fastest (matches rect_grid / np.mgrid flatten)
        xs = np.tile(x, ny + 1)
        ys = np.repeat(y, nx + 1)

        assert xs.size == n_nodes
        idx = right_edge_mid_height_index(nx, ny)
        assert 0 <= idx < n_nodes
        assert np.isclose(xs[idx], L / 2)
        # mid-height: for even ny, ny//2 is exact mid row of ny+1 points
        mid_y = y[ny // 2]
        assert np.isclose(ys[idx], mid_y)

        wrong = nx * (ny // 2 + 1) - 1
        # Wrong index is not on the right edge for these meshes
        assert not np.isclose(xs[wrong], L / 2) or wrong == idx


def test_wrong_legacy_index_is_interior_for_typical_mesh():
    nx, ny = 6, 2
    L = 8.0
    x = np.linspace(-L / 2, L / 2, nx + 1)
    xs = np.tile(x, ny + 1)
    wrong = nx * (ny // 2 + 1) - 1
    correct = right_edge_mid_height_index(nx, ny)
    assert wrong != correct
    assert not np.isclose(xs[wrong], L / 2)
    assert np.isclose(xs[correct], L / 2)
