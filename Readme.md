# entropy_2pt_method

Absolute entropy of rigid water from a GROMACS trajectory with the two-phase thermodynamic (2PT)
model of Lin, Maiti and Goddard,
[J. Phys. Chem. B, Vol. 114, No. 24, 2010](https://doi.org/10.1021/jp103120q). The script reads one `.tpr` + `.trr` pair that contains
only water, and returns the translational, rotational and total entropy per molecule, the two
fluidicities and the density of states (DoS).

It is a single, short Python file that follows the equations of the paper step by step, with the
paper's equation numbers in the comments.

## Requirements

- Python 3
- NumPy (1.x or 2.x), SciPy, MDAnalysis
- MDAnalysis ≥ 2.11 to read `.tpr` files written by GROMACS 2026. With an older MDAnalysis, pass a
  `.gro` of the same atoms to `-s` instead; the masses are then guessed from the atom names
  (O 15.999, H 1.008), which changes the entropy by about 0.001 J mol⁻¹ K⁻¹.

```bash
pip install numpy scipy MDAnalysis
```

## Usage

```bash
python3 entropy_2pt_method.py -s water.tpr -f water.trr -t 303.15 -o water
```

| Option | Meaning | Default |
| --- | --- | --- |
| `-s` | topology (`.tpr`, or `.gro`) of the atoms in the `.trr` | `gap.tpr` |
| `-f` | trajectory with positions **and velocities** | `gap.trr` |
| `-t` | temperature of the simulation (K) | `303.15` |
| `-o` | basename of the output files | none |

`python3 entropy_2pt_method.py -h` lists the options.

## Input

- The `.tpr` and `.trr` must contain the same atoms, all of them water, in the order O, H, H for
  every molecule.
- The water must be rigid (e.g. TIP3P or SPC/E with SETTLE): there is no vibrational term.
- Every frame must contain velocities, equally spaced in time. Frames every 8 fs (Nyquist frequency
  2085 cm⁻¹) cover the librational band of water; 20-ps trajectories (2501 frames) are enough.
- The box must be rectangular and constant (NVT); it is read from the first frame.

To cut a group of water out of a larger system, write a `.tpr` and a `.trr` of the same index group:

```bash
echo Water | gmx convert-tpr -s topol.tpr -n index.ndx -o water.tpr
echo Water | gmx trjconv     -s topol.tpr -f traj.trr -n index.ndx -o water.trr
```

with, for example, `nstxout = nstvout = 4` and `dt = 0.002` in the `.mdp` file.

**Volume per molecule.** The fluidicity (eq 13) and the hard-sphere entropy (eq 16) need the volume
per molecule, `V_PER_MOLECULE` near the top of the script. It is set to bulk TIP3P at 303 K
(0.02981 nm³). For a pure water box, use box volume / number of molecules. For a subset of water in
an inhomogeneous system the volume is not well defined; a ±20 % change shifts the entropy by about
0.2 J mol⁻¹ K⁻¹.

**Memory.** All frames are held in memory, about 0.7 kB per molecule and frame (2.3 GB for 1328
molecules × 2501 frames).

## Output

The script prints its progress and the result:

```
1328 water molecules, 2501 frames, dt = 8.0 fs, masses [15.99940014  1.00800002  1.00800002]
moments of inertia I = [0.01155114 0.00614571 0.01769685] amu nm^2
check eq 4: integral of S_trans = 3.000 N, of S_rot = 3.012 N   (both should be close to 3)
translation: Delta = 0.3567  f = 0.3491  S_trans = 56.252 J/(mol K)
rotation:    Delta = 0.0206  f = 0.0803  S_rot   = 13.825 J/(mol K)
total:       S = 70.077 J/(mol K) per water molecule
wrote water_entropy.dat and water_dos.dat
```

and writes two files (`entropy.dat` and `dos.dat` without `-o`):

| File | Content |
| --- | --- |
| `<basename>_entropy.dat` | one line: S_total, S_trans, S_rot (J mol⁻¹ K⁻¹ per molecule), f_trans, f_rot |
| `<basename>_dos.dat` | three columns: ν (ps⁻¹; × 33.356 for cm⁻¹), S_trans(ν)/N, S_rot(ν)/N (ps) |

The DoS is given per molecule, so each column integrates to 3. Both files start with a `#` header
that names the input files, the temperature and the number of molecules.

## Method

Equation numbers refer to Lin, Maiti and Goddard (2010).

1. **Velocities** (eqs 27, 28). For every molecule and frame: the centre-of-mass velocity, and the
   angular velocity about the three principal axes (H–O–H bisector, H–H direction, normal to the
   plane), obtained from the angular momentum and the principal moments of inertia. Bond vectors
   are taken in the minimum-image convention.
2. **Density of states** (eqs 1, 2, 20). S_trn(ν) and S_rot(ν) from the fast Fourier transform of
   the centre-of-mass and angular velocities, weighted with the molecular mass and the moments of
   inertia. Each integrates to 3N (eq 4), which the script checks.
3. **Fluidicity** (eqs 12, 13). The dimensionless diffusivity Δ from S(0), and the fluidicity f as
   the root of eq 12 (Brent's method), separately for translation and rotation.
4. **Gas–solid split** (eqs 10, 11, 21, 22). The gas-like DoS is a Lorentzian that starts at S(0) and
   holds 3fN modes; the solid-like DoS is the rest.
5. **Weighting functions** (eqs 8b, 15b–17, 25b). Quantum harmonic oscillator for the solid-like
   modes, hard-sphere gas (Carnahan–Starling) for gas-like translation, free rigid rotor (σ = 2)
   for gas-like rotation.
6. **Entropy** (eqs 23b, 24b). S = S_trn + S_rot, with S_vib = 0. The solid-like part is integrated
   with the trapezoidal rule up to the Nyquist frequency; the gas-like part analytically, as
   3fN·W_gas.

## Check

Bulk-like water (1328 molecules, one 20-ps segment at 303.15 K) compared with TIP3P in Lin et al.
(2010), Table 1, at 298 K:

| | T (K) | f_trans | f_rot | S_trans | S_rot | S (J mol⁻¹ K⁻¹) |
| --- | --- | --- | --- | --- | --- | --- |
| this script | 303.15 | 0.349 | 0.080 | 56.25 | 13.83 | 70.08 |
| Lin et al. (2010) | 298 | 0.34 | 0.08 | 55.59 | 12.90 | 68.49 |

## Reference

S.-T. Lin, P. K. Maiti, W. A. Goddard III, Two-phase thermodynamic model for efficient and accurate
absolute entropy of water from molecular dynamics simulations,
[J. Phys. Chem. B, Vol. 114, No. 24, 2010](https://doi.org/10.1021/jp103120q).

 
