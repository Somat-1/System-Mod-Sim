# Guyan Model Reduction

This folder implements the supplied classical Guyan reduction without
alteration. The retained coordinates are theta_m and x_n; the eliminated
coordinates are theta_c, theta_s, theta_sb and x_s.

The transformation is evaluated from the supplied closed forms for beta,
kappa, k_ch, nu and mu. The reduced matrices use only the stated Galerkin
projection:

    M_r = T^T M T
    C_r = T^T C T
    K_r = T^T K T
    b_r = T^T b

The command column is b = [k_EM, 0, ..., 0]^T, exactly as specified in the
supplied derivation and in Rev 4/state_space_6dof.md. Neither the optional
diagonal mass tuning nor the IRS correction is applied.

## Files

- guyan_model.py implements the full matrices, closed-form transformation,
  projection, port Jacobians, fixed-interface modes and frequency response.
- generate_guyan_bode.py verifies every closed form against numerical static
  condensation, checks DC exactness, and writes both the standalone and
  full-frictionless overlay Bode plots, response arrays and numerical summary.
- rendered_assets/guyan_bode.png is the untuned two-master Bode and phase
  response; its complete numerical response is retained through 8 kHz.
- rendered_assets/guyan_vs_frictionless_bode.png overlays the two-master
  Guyan response with the documented full 6-DOF frictionless baseline.
- rendered_assets/guyan_model_summary.json contains the complete audit.

## Overlay convention

Both overlay responses use the physical command vector
b = [k_EM, 0, ..., 0]^T specified by the supplied derivation and
Rev 4/state_space_6dof.md. The older Rev 4/scripts/build_bode_rev4.py
uses k_EM + k_d in that command column; its pre-existing plot is therefore
not used as the baseline curve in this like-for-like comparison.
