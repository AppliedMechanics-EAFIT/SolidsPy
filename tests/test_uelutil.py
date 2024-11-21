# -*- coding: utf-8 -*-
"""
Test cases for functions on ``uelutil`` module

"""
import numpy as np
import src.uelutil as uel


#%% Continuum elements
def test_elast_tri3():
    """Tests for 3-noded tri uel"""
    coord = np.array([
        [0, 0],
        [1, 0],
        [0, 1]])
    params = 8/3, 1/3, 1
    stiff, mass = uel.elast_tri3(coord, params)
    stiff_ex = 1/2 * np.array([
       [4, 2, -3, -1, -1, -1],
       [2, 4, -1, -1, -1, -3],
       [-3, -1, 3, 0, 0, 1],
       [-1, -1, 0, 1, 1, 0],
       [-1, -1, 0, 1, 1, 0],
       [-1, -3, 1, 0, 0, 3]])
    mass_ex = 1/24 * np.array([
        [2, 0, 1, 0, 1, 0],
        [0, 2, 0, 1, 0, 1],
        [1, 0, 2, 0, 1, 0],
        [0, 1, 0, 2, 0, 1],
        [1, 0, 1, 0, 2, 0],
        [0, 1, 0, 1, 0, 2]])
    assert np.allclose(stiff, stiff_ex)
    assert np.allclose(mass, mass_ex)


def test_elast_tri6():
    """Tests for 6-noded tri uel"""
    coord = np.array([
        [0, 0],
        [1, 0],
        [0, 1],
        [0.5, 0],
        [0.5, 0.5],
        [0, 0.5]])
    params = 8/3, 1/3, 1
    stiff, mass = uel.elast_tri6(coord, params)
    stiff_ex = 1/6 * np.array([
        [12, 6, 3, 1, 1, 1, -12, -4, 0, 0, -4, -4],
        [6, 12, 1, 1, 1, 3, -4, -4, 0, 0, -4, -12],
        [3, 1, 9, 0, 0, -1, -12, -4, 0, 4, 0, 0],
        [1, 1, 0, 3, -1, 0, -4, -4, 4, 0, 0, 0],
        [1, 1, 0, -1, 3, 0, 0, 0, 0, 4, -4, -4],
        [1, 3, -1, 0, 0, 9, 0, 0, 4, 0, -4, -12],
        [-12, -4, -12, -4, 0, 0, 32, 8, -8, -8, 0, 8],
        [-4, -4, -4, -4, 0, 0, 8, 32, -8, -24, 8, 0],
        [0, 0, 0, 4, 0, 4, -8, -8, 32, 8, -24, -8],
        [0, 0, 4, 0, 4, 0, -8, -24, 8, 32, -8, -8],
        [-4, -4, 0, 0, -4, -4, 0, 8, -24, -8, 32, 8],
        [-4, -12, 0, 0, -4, -12, 8, 0, -8, -8, 8, 32]])
    mass_ex = 1/360 * np.array([
        [6, 0, -1, 0, -1, 0, 0, 0, -4, 0, 0, 0],
        [0, 6, 0, -1, 0, -1, 0, 0, 0, -4, 0, 0],
        [-1, 0, 6, 0, -1, 0, 0, 0, 0, 0, -4, 0],
        [0, -1, 0, 6, 0, -1, 0, 0, 0, 0, 0, -4],
        [-1, 0, -1, 0, 6, 0, -4, 0, 0, 0, 0, 0],
        [0, -1, 0, -1, 0, 6, 0, -4, 0, 0, 0, 0],
        [0, 0, 0, 0, -4, 0, 32, 0, 16, 0, 16, 0],
        [0, 0, 0, 0, 0, -4, 0, 32, 0, 16, 0, 16],
        [-4, 0, 0, 0, 0, 0, 16, 0, 32, 0, 16, 0],
        [0, -4, 0, 0, 0, 0, 0, 16, 0, 32, 0, 16],
        [0, 0, -4, 0, 0, 0, 16, 0, 16, 0, 32, 0],
        [0, 0, 0, -4, 0, 0, 0, 16, 0, 16, 0, 32]])
    assert np.allclose(stiff, stiff_ex)
    assert np.allclose(mass, mass_ex)


# Quadrilaterals
def test_elast_quad4():
    """Tests for 4-noded quad uel"""
    coord = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
    params = 8/3, 1/3
    stiff, mass = uel.elast_quad4(coord, params)
    stiff_ex = 1/6 * np.array([
         [ 8,  3, -5,  0, -4, -3,  1,  0],
         [ 3,  8,  0,  1, -3, -4,  0, -5],
         [-5,  0,  8, -3,  1,  0, -4,  3],
         [ 0,  1, -3,  8,  0, -5,  3, -4],
         [-4, -3,  1,  0,  8,  3, -5,  0],
         [-3, -4,  0, -5,  3,  8,  0,  1],
         [ 1,  0, -4,  3, -5,  0,  8, -3],
         [ 0, -5,  3, -4,  0,  1, -3,  8]])
    mass_ex = 1/9 * np.array([
        [4, 0, 2, 0, 1, 0, 2, 0],
        [0, 4, 0, 2, 0, 1, 0, 2],
        [2, 0, 4, 0, 2, 0, 1, 0],
        [0, 2, 0, 4, 0, 2, 0, 1],
        [1, 0, 2, 0, 4, 0, 2, 0],
        [0, 1, 0, 2, 0, 4, 0, 2],
        [2, 0, 1, 0, 2, 0, 4, 0],
        [0, 2, 0, 1, 0, 2, 0, 4]])
    assert np.allclose(stiff, stiff_ex)
    assert np.allclose(mass, mass_ex)


#%% Structural elements
def test_spring():
    """Tests for 2-noded springs uel"""
    coord = np.array([
        [0, 0],
        [1, 0]])
    params = 1
    stiff, _ = uel.spring(coord, params)
    stiff_ex = np.array([
        [1, 0, -1, 0],
        [0, 0, 0, 0],
        [-1, 0, 1, 0],
        [0, 0, 0, 0]])
    assert np.allclose(stiff, stiff_ex)


def test_truss2D():
    """Tests for 2-noded 2D truss uel"""
    coord = np.array([
        [0, 0],
        [1, 0]])
    params = 1.0, 1.0, 1.0
    stiff, mass = uel.truss2D(coord, params)
    stiff_ex =  np.array([
        [1, 0, -1, 0],
        [0, 0, 0, 0],
        [-1, 0, 1, 0],
        [0, 0, 0, 0]])
    mass_ex =  1/6*np.array([
        [2, 0, 1, 0],
        [0, 2, 0, 1],
        [1, 0, 2, 0],
        [0, 1, 0, 2]])
    assert np.allclose(stiff, stiff_ex)
    assert np.allclose(mass, mass_ex)


def test_beam2DU():
    """Tests for 2-noded 2D beam uel"""
    coord = np.array([
        [0, 0],
        [1, 0]])
    params = 1.0, 1.0, 1.0, 1.0
    stiff, mass = uel.beam2DU(coord, params)
    stiff_ex =  np.array([
        [0, 0, 0, 0, 0, 0],
        [0, 12, 6, 0, -12, 6],
        [0, 6, 4, 0, -6, 2],
        [0, 0, 0, 0, 0, 0],
        [0, -12, -6, 0, 12, -6],
        [0, 6, 2, 0, -6, 4]])
    mass_ex = 1/420*np.array([
        [0, 0, 0, 0, 0, 0],
        [0, 156, 22, 0, 54, -13],
        [0, 22, 4, 0, 13, -3],
        [0, 0, 0, 0, 0, 0],
        [0, 54, 13, 0, 156, -22],
        [0, -13, -3, 0, -22, 4]])
    assert np.allclose(stiff, stiff_ex)
    assert np.allclose(mass, mass_ex)
