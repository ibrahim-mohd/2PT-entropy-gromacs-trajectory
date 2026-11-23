# Mohd Ibrahim, Technical University of Munich
# Claude-code Opus-5.5 (xhigh)
"""
entropy_2pt_method.py -- 2PT entropy of rigid water from one GROMACS .tpr + .trr pair
=====================================================================================

Method: S.-T. Lin, P. K. Maiti, W. A. Goddard III, J. Phys. Chem. B, Vol. 114, No. 24, 2010
        https://doi.org/10.1021/jp103120q
        "eq N" below refers to the equations of that paper.

What it does
  The .trr must contain velocities (here they were written every 8 fs). The script
  1. splits the motion of every water molecule into translation (centre of mass) and rotation,
  2. Fourier-transforms these velocities -> density of states (DoS) S(nu),
  3. splits each DoS into a gas-like and a solid-like part (2PT),
  4. weights each part with its entropy per mode and integrates -> S in J/(mol K) per molecule.
  Rigid water (TIP3P with SETTLE) has no intramolecular vibration, so S = S_trans + S_rot.

How to use
  python3 entropy_2pt_method.py -s gap.tpr -f gap.trr -t 303.15 -o gap
    -s  topology of the atoms in the .trr      (default gap.tpr)
    -f  trajectory with positions + velocities (default gap.trr)
    -t  temperature of the simulation in K     (default 303.15)
    -o  basename of the output files           (default: none)
  Prints the result and writes (with -o gap: gap_entropy.dat and gap_dos.dat)
    entropy.dat  S_total, S_trans, S_rot (J/(mol K) per molecule), f_trans, f_rot
    dos.dat      nu (1/ps), S_trans(nu)/N, S_rot(nu)/N (density of states per molecule, ps)
  The .tpr and .trr must contain the same atoms, all of them water (O, H, H per molecule).
  To make such a pair for one region of one segment, see README.txt in this folder.
  Needs numpy, scipy and MDAnalysis. GROMACS 2026 .tpr files need MDAnalysis >= 2.11;
  with an older MDAnalysis, give a .gro of the same atoms as TPR instead (the masses are
  then guessed from the atom names: O 15.999, H 1.008 -- the result changes by ~0.001).

Units: nm, ps, amu, kJ/mol (the GROMACS units). Entropy is printed in J/(mol K) per molecule.
"""
import argparse

import numpy as np
import MDAnalysis as mda
from scipy.optimize import brentq

# trapezoidal integration: np.trapezoid in NumPy >= 2.0, np.trapz in older versions
trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

# ----------------------------------------------------------------- INPUT --------
parser = argparse.ArgumentParser(description="2PT entropy of rigid water from one GROMACS .tpr + .trr pair")
parser.add_argument("-s", default="gap.tpr", help="topology (atom names, masses) of the atoms in the .trr (default: gap.tpr)")
parser.add_argument("-f", default="gap.trr", help="positions + velocities of the same atoms (default: gap.trr)")
parser.add_argument("-t", type=float, default=303.15, help="temperature of the simulation in K (default: 303.15)")
parser.add_argument("-o", default="", help="basename of the output files: <o>_entropy.dat and <o>_dos.dat "
                                           "(default: entropy.dat and dos.dat)")
args = parser.parse_args()
TPR = args.s                          # topology (atom names, masses) of the atoms in TRR
TRR = args.f                          # positions + velocities of the same atoms
T = args.t                            # temperature of the simulation (K)
OUT_ENTROPY = f"{args.o}_entropy.dat" if args.o else "entropy.dat"   # output files
OUT_DOS = f"{args.o}_dos.dat" if args.o else "dos.dat"
V_PER_MOLECULE = 5.41332**3 / 5322    # volume per water molecule (nm^3) used in eqs 13, 16:
                                      # bulk TIP3P at 303 K (a 5322-water box of 5.41332 nm).
                                      # For a pure water box you would use box volume / N.
# --------------------------------------------------------------------------------

kB = 0.0083144626      # Boltzmann constant, kJ/(mol K)     (= amu nm^2 ps^-2 K^-1)
h = 0.399031271        # Planck constant x Avogadro, kJ/mol ps
R = 8.314462618        # gas constant, J/(mol K): converts "entropy in units of k" to J/(mol K)


# ------------------------------------------------------------------------------
# 1. Read all frames: x = positions, v = velocities, shape (frame, molecule, atom, xyz)
# ------------------------------------------------------------------------------
u = mda.Universe(TPR, TRR)
water = u.atoms
assert water.n_atoms % 3 == 0, "expected 3 atoms (O, H, H) per molecule"
N = water.n_atoms // 3                       # number of water molecules
nframes = len(u.trajectory)

x = np.zeros((nframes, N, 3, 3))
v = np.zeros((nframes, N, 3, 3))
for i, ts in enumerate(u.trajectory):
    x[i] = water.positions.reshape(N, 3, 3) / 10.0      # Angstrom -> nm
    v[i] = water.velocities.reshape(N, 3, 3) / 10.0     # Angstrom/ps -> nm/ps
box = u.trajectory[0].dimensions[:3] / 10.0             # box lengths (nm)
dt = u.trajectory[1].time - u.trajectory[0].time        # time between frames (ps)

m = water.masses[:3].astype(float)    # masses of O, H1, H2 (amu)
M = m.sum()                           # mass of one water molecule
print(f"{N} water molecules, {nframes} frames, dt = {1000 * dt:.1f} fs, masses {m}")


# ------------------------------------------------------------------------------
# 2. Translation and rotation of every molecule in every frame (eqs 27, 28)
# ------------------------------------------------------------------------------
# O->H bond vectors; a molecule split by the periodic box is made whole (minimum image)
d1 = x[:, :, 1] - x[:, :, 0] 
d1 -= box * np.round(d1 / box)
d2 = x[:, :, 2] - x[:, :, 0]
d2 -= box * np.round(d2 / box)

# atom positions relative to the centre of mass (COM) of their molecule
com = (m[1] * d1 + m[2] * d2) / M                        # COM seen from the O atom
r = np.stack([-com, d1 - com, d2 - com], axis=2)         # (frame, molecule, atom, xyz)

# translation: centre-of-mass velocity of each molecule
v_com = (m[:, None] * v).sum(axis=2) / M                 # (frame, molecule, xyz)

# rotation: angular momentum L = sum over atoms of m r x (v - v_com)          (eq 28)
L = (m[:, None] * np.cross(r, v - v_com[:, :, None, :])).sum(axis=2)

# principal axes of a rigid water molecule:
#   e1 = H-O-H bisector, e2 = along H1-H2, e3 = perpendicular to the molecular plane
e1 = d1 + d2
e1 /= np.linalg.norm(e1, axis=-1, keepdims=True)
e2 = d1 - d2
e2 /= np.linalg.norm(e2, axis=-1, keepdims=True)
e3 = np.cross(e1, e2)
axes = [e1, e2, e3]

# principal moments of inertia I_k = sum over atoms of m (|r|^2 - (r.e_k)^2);
# constant for rigid water, so the average over all frames and molecules is used
r2 = (r * r).sum(-1)                                     # |r|^2, shape (frame, molecule, atom)
I = np.zeros(3)
for k, e in enumerate(axes):
    r_along_e = (r * e[:, :, None, :]).sum(-1)           # r . e_k, shape (frame, molecule, atom)
    I[k] = (m * (r2 - r_along_e ** 2)).sum(axis=2).mean()

# angular velocity about each principal axis: omega_k = (L . e_k) / I_k        (eq 28)
omega = np.stack([(L * e).sum(-1) / I[k] for k, e in enumerate(axes)], axis=-1)   # (frame, molecule, 3)
print(f"moments of inertia I = {I} amu nm^2")


# ------------------------------------------------------------------------------
# 3. Density of states S(nu) from the Fourier transform of the velocities (eqs 1, 2, 20)
#    S(nu) = 2/(kT) * sum over molecules and directions of  weight * |FT of q(t)|^2 / t_total
#    weight = molecule mass M for translation, moment of inertia I_k for rotation.
#    The integral of S(nu) over nu is the number of degrees of freedom, 3N (eq 4).
# ------------------------------------------------------------------------------
nu = np.fft.rfftfreq(nframes, dt)                  # frequencies 0 ... 1/(2 dt), in 1/ps
t_total = nframes * dt


def density_of_states(q, weight):
    q_hat = np.fft.rfft(q, axis=0) * dt            # approximates  integral q(t) exp(-2 pi i nu t) dt
    return 2.0 / (kB * T) * (weight * np.abs(q_hat) ** 2).sum(axis=(1, 2)) / t_total

S_trans = density_of_states(v_com, np.array([M, M, M]))
S_rot = density_of_states(omega, I)
print(f"check eq 4: integral of S_trans = {trapezoid(S_trans, nu) / N:.3f} N, "
      f"of S_rot = {trapezoid(S_rot, nu) / N:.3f} N   (both should be close to 3)")


# ------------------------------------------------------------------------------
# 4. Fluidicity f (eqs 12, 13): the fraction of the 3N modes that is gas-like.
#    It follows from S(0) (diffusivity) and the density N/V. Translation and rotation
#    each get their own f, with s0 = S_trans(0) or S_rot(0) (text after eq 22).
# ------------------------------------------------------------------------------

def fluidicity(s0):
    # eq 13: dimensionless diffusivity Delta
    Delta = (2 * s0 / (9 * N)) * np.sqrt(np.pi * kB * T / M) * (1 / V_PER_MOLECULE) ** (1 / 3) * (6 / np.pi) ** (2 / 3)

    # eq 12: f is the root between 0 and 1 of
    def eq12(f):
        return (2 * Delta ** -4.5 * f ** 7.5 - 6 * Delta ** -3 * f ** 5 - Delta ** -1.5 * f ** 3.5
                + 6 * Delta ** -1.5 * f ** 2.5 + 2 * f - 2)
    return brentq(eq12, 1e-10, 1.0, xtol=1e-14), Delta


# ------------------------------------------------------------------------------
# 5. Entropy weighting functions (entropy per mode, in units of k)
# ------------------------------------------------------------------------------
# solid-like modes = quantum harmonic oscillators (eq 8b); u = h nu / kT
W_solid = np.zeros_like(nu)
u_ho = h * nu[1:] / (kB * T)
W_solid[1:] = u_ho / np.expm1(u_ho) - np.log(1 - np.exp(-u_ho))
# (W_solid diverges at nu = 0, but the solid-like DoS is exactly 0 there, so that point adds nothing)


# ------------------------------------------------------------------------------
# 6. Translation: gas-like part = hard-sphere gas (eqs 11, 15b, 16, 17)
# ------------------------------------------------------------------------------
s0 = S_trans[0]
f_trans, Delta_trans = fluidicity(s0)
S_trans_gas = s0 / (1 + (np.pi * s0 * nu / (6 * f_trans * N)) ** 2)     # eq 11 (a Lorentzian)
S_trans_solid = S_trans - S_trans_gas                                   # eq 10

y = f_trans ** 2.5 / Delta_trans ** 1.5                                  # hard-sphere packing fraction
z = (1 + y + y ** 2 - y ** 3) / (1 - y) ** 3                             # eq 17, Carnahan-Starling
S_HS = (2.5 + np.log((2 * np.pi * M * kB * T / h ** 2) ** 1.5 * V_PER_MOLECULE / f_trans * z)
        + y * (3 * y - 4) / (1 - y) ** 2)                                # eq 16, hard-sphere entropy / k
W_gas_trans = S_HS / 3                                                   # eq 15b

# eq 24b:  S/k = integral(S_solid W_solid) + integral(S_gas W_gas)
# W_gas is a constant and the Lorentzian S_gas integrates to exactly 3 f N over 0..infinity,
# so the gas term is 3 f N W_gas (the numerical spectrum only reaches 1/(2 dt)).

s_trans = (trapezoid(S_trans_solid * W_solid, nu) + 3 * f_trans * N * W_gas_trans) / N

# ------------------------------------------------------------------------------
# 7. Rotation: gas-like part = free rigid rotor (eqs 21, 22, 25b)
# ------------------------------------------------------------------------------
s0 = S_rot[0]
f_rot, Delta_rot = fluidicity(s0)
S_rot_gas = s0 / (1 + (np.pi * s0 * nu / (6 * f_rot * N)) ** 2)        # eq 22
S_rot_solid = S_rot - S_rot_gas                                         # eq 21

theta = h ** 2 / (8 * np.pi ** 2 * I * kB)                              # rotational temperatures (K)
sigma = 2                                                               # symmetry number of water
S_R = np.log(np.sqrt(np.pi) * np.exp(1.5) / sigma * np.sqrt(T ** 3 / np.prod(theta)))   # rigid rotor S/k
W_gas_rot = S_R / 3                                                     # eq 25b

s_rot = (trapezoid(S_rot_solid * W_solid, nu) + 3 * f_rot * N * W_gas_rot) / N


# ------------------------------------------------------------------------------
# 8. Result (eq 23b): per molecule, in J/(mol K)
# ------------------------------------------------------------------------------
print(f"translation: Delta = {Delta_trans:.4f}  f = {f_trans:.4f}  S_trans = {R * s_trans:.3f} J/(mol K)")
print(f"rotation:    Delta = {Delta_rot:.4f}  f = {f_rot:.4f}  S_rot   = {R * s_rot:.3f} J/(mol K)")
print(f"total:       S = {R * (s_trans + s_rot):.3f} J/(mol K) per water molecule")

# entropy.dat: one line with the entropies (J/(mol K) per molecule) and the fluidicities
with open(OUT_ENTROPY, "w") as fh:
    fh.write(f"# 2PT entropy per water molecule; tpr = {TPR}, trr = {TRR}, T = {T} K, N = {N}\n")
    fh.write("# S_total (J/(mol K))  S_trans (J/(mol K))  S_rot (J/(mol K))  f_trans  f_rot\n")
    fh.write(f"{R * (s_trans + s_rot):.4f}  {R * s_trans:.4f}  {R * s_rot:.4f}  {f_trans:.6f}  {f_rot:.6f}\n")

# dos.dat: density of states per molecule (each integrates to 3 over nu)
np.savetxt(OUT_DOS, np.column_stack([nu, S_trans / N, S_rot / N]), fmt="%.6e",
           header=f"density of states per water molecule, S(nu)/N in ps; tpr = {TPR}, trr = {TRR}, T = {T} K\n"
                  "nu (1/ps; x 33.356 for cm^-1)  S_trans(nu)/N  S_rot(nu)/N")
print(f"wrote {OUT_ENTROPY} and {OUT_DOS}")
