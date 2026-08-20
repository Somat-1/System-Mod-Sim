#!/usr/bin/env python3
"""Shared color palettes for the mode-shape/strain-energy/kinetic-energy
figures, requested 2026-08-20. Derived from four given anchors

    #454040 (69,64,64)  #605B51 (96,91,81)  #D8D365 (216,211,101)  #E6F082 (230,240,130)

by walking the piecewise-linear path through those anchors in CIE Lab space
(equal-arc-length sampling, not equal-RGB-t -- the anchors themselves are
perceptually uneven: dE(#454040,#605B51)=12.8 and dE(#D8D365,#E6F082)=11.2,
but dE(#605B51,#D8D365)=66.3, so naive linear interpolation crowds new
stops into the two short segments and wastes the long one), then adding a
small alternating perpendicular offset so adjacent stops clear a ~10-25 dE
gap (minimum pairwise dE ~11) instead of ~6, which plain arc-length
sampling alone still left too close for reliable bar-to-bar distinction.

DOF_COLORS order matches build_bode_rev4.STATE_LABELS
    [theta_m, theta_c, theta_s, theta_sb, x_s, x_n]
ELEMENT_COLORS order matches the spring decomposition used in
plot_modal_strain_energy.py / plot_stiffness_root_locus.py
    [k_EM, k_d, k_c, k_s1, k_s2, k_brg, k_nut]
"""

from __future__ import annotations

DOF_COLORS = [
    "#645756",  # theta_m
    "#535647",  # theta_c
    "#938564",  # theta_s
    "#9FA254",  # theta_sb
    "#E4D671",  # x_s
    "#E1FF8D",  # x_n
]

ELEMENT_COLORS = [
    "#645756",  # k_EM
    "#4C4E42",  # k_d
    "#867861",  # k_c
    "#8C8E52",  # k_s1
    "#C9BA6D",  # k_s2
    "#CCD35E",  # k_brg
    "#FFFFA6",  # k_nut
]
