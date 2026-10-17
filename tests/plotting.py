#!/usr/bin/env python3
"""Check that conversion of shell sums preserves the integrated wave action."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from gp2d_plotting import density_in_frequency


for bin_width in (1., .5, 2. / 3.):
    k = np.arange(8) * bin_width
    shells = np.array([0., .3, .1, .8, .2, .5, .4, .1])
    density = density_in_frequency(k, shells, -1., bin_width=bin_width)
    # d omega = |d omega/d k| d k for the same Jacobian discretization.
    recovered = np.sum(density[1:] * 2. * k[1:] * bin_width)
    assert np.isclose(recovered, shells.sum(), rtol=1e-14)

print("Frequency-density conversion preserves shell-integrated wave action")
