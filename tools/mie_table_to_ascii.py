#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
# SPDX-License-Identifier: GPL-3.0-or-later
#
# This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).
"""
Convert an HDF5 grid of Mie-scattering particle properties for log-normal
size distributions into the plain-text table read by BeAR's ``mie`` cloud
model.

The HDF5 file is expected to contain the datasets

  r_g_micron                (nb_r)                 median radius of the distribution
  sigma_g                   (nb_sigma)             geometric standard deviation
  wavelength_micron         (nb_wavelength)        wavelength grid
  c_ext_cm2                 (nb_r, nb_sigma, nb_wavelength)  mean extinction cross section per particle
  single_scattering_albedo  (nb_r, nb_sigma, nb_wavelength)
  asymmetry_parameter       (nb_r, nb_sigma, nb_wavelength)

The ASCII output has the layout (lines starting with '#' are comments and
can appear anywhere):

  nb_r nb_sigma nb_wavelength
  r_g values            (one line, micron, ascending)
  sigma_g values        (one line, ascending)
  wavelength values     (one line, micron, ascending)
  for each r_g (outer loop) and each sigma_g (inner loop):
    c_ext values   (one line, cm^2)
    ssa values     (one line)
    g values       (one line)

Example:
  python3 mie_table_to_ascii.py mg2sio4_optics_grid.h5 mg2sio4_mie.dat --wavelength-min 0.3 --wavelength-max 30
"""

import argparse
import sys

import numpy as np

try:
  import h5py
except ImportError:
  sys.exit("This script requires the h5py package (pip install h5py)")


def main():
  parser = argparse.ArgumentParser(
    description="Convert an HDF5 Mie property grid into BeAR's ASCII cloud table format")
  parser.add_argument("input", help="HDF5 file with the tabulated Mie properties")
  parser.add_argument("output", help="ASCII output file")
  parser.add_argument("--wavelength-min", type=float, default=None,
    help="drop wavelengths below this value (micron) to reduce the file size")
  parser.add_argument("--wavelength-max", type=float, default=None,
    help="drop wavelengths above this value (micron) to reduce the file size")
  parser.add_argument("--precision", type=int, default=6,
    help="number of significant digits after the decimal point (default: 6)")
  args = parser.parse_args()

  with h5py.File(args.input, "r") as f:
    r_g = np.asarray(f["r_g_micron"][:], dtype=float)
    sigma_g = np.asarray(f["sigma_g"][:], dtype=float)
    wavelength = np.asarray(f["wavelength_micron"][:], dtype=float)
    c_ext = np.asarray(f["c_ext_cm2"][:], dtype=float)
    ssa = np.asarray(f["single_scattering_albedo"][:], dtype=float)
    g = np.asarray(f["asymmetry_parameter"][:], dtype=float)
    description = f.attrs.get("description", "")
    refractive_index = f.attrs.get("refractive_index", "")
    size_distribution = f.attrs.get("size_distribution", "")

  expected_shape = (r_g.size, sigma_g.size, wavelength.size)

  for name, data in (("c_ext_cm2", c_ext), ("single_scattering_albedo", ssa), ("asymmetry_parameter", g)):
    if data.shape != expected_shape:
      sys.exit(f"Dataset {name} has shape {data.shape}, expected {expected_shape}")

  #BeAR requires strictly ascending axes
  for name, axis in (("r_g_micron", r_g), ("sigma_g", sigma_g), ("wavelength_micron", wavelength)):
    if np.any(np.diff(axis) <= 0):
      sys.exit(f"Axis {name} is not strictly ascending")

  #optional wavelength crop
  mask = np.ones(wavelength.size, dtype=bool)

  if args.wavelength_min is not None:
    mask &= wavelength >= args.wavelength_min

  if args.wavelength_max is not None:
    mask &= wavelength <= args.wavelength_max

  if mask.sum() < 2:
    sys.exit("Fewer than two wavelengths left after cropping")

  wavelength = wavelength[mask]
  c_ext = c_ext[:, :, mask]
  ssa = ssa[:, :, mask]
  g = g[:, :, mask]

  if np.any(c_ext <= 0):
    sys.exit("The extinction cross sections must be positive (they are interpolated in log space)")

  if not np.all(np.isfinite(c_ext)) or not np.all(np.isfinite(ssa)) or not np.all(np.isfinite(g)):
    sys.exit("The table contains non-finite values")

  fmt_e = f"%.{args.precision}e"
  fmt_f = f"%.{args.precision}f"

  def line(values, fmt):
    return " ".join(fmt % v for v in values) + "\n"

  with open(args.output, "w") as out:
    out.write("# BeAR Mie cloud particle property table\n")
    out.write(f"# source file: {args.input}\n")

    if description:
      out.write(f"# {description}\n")
    if refractive_index:
      out.write(f"# refractive index: {refractive_index}\n")
    if size_distribution:
      out.write(f"# size distribution: {size_distribution}\n")

    out.write("# layout: nb_r nb_sigma nb_wavelength; r_g [micron]; sigma_g; wavelength [micron];\n")
    out.write("#         then for each r_g and each sigma_g: c_ext [cm^2], single scattering albedo, asymmetry parameter\n")
    out.write(f"{r_g.size} {sigma_g.size} {wavelength.size}\n")
    out.write("# r_g [micron]\n")
    out.write(line(r_g, fmt_e))
    out.write("# sigma_g\n")
    out.write(line(sigma_g, fmt_e))
    out.write("# wavelength [micron]\n")
    out.write(line(wavelength, fmt_e))

    for i in range(r_g.size):
      for j in range(sigma_g.size):
        out.write(f"# r_g = {r_g[i]:.6e} sigma_g = {sigma_g[j]:.6e}\n")
        out.write(line(c_ext[i, j, :], fmt_e))
        out.write(line(ssa[i, j, :], fmt_f))
        out.write(line(g[i, j, :], fmt_f))

  print(f"Wrote {args.output}: {r_g.size} radii, {sigma_g.size} sigma_g values, {wavelength.size} wavelengths")


if __name__ == "__main__":
  main()
