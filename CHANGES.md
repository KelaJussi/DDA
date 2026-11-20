# Update: optimal vs single-mode excitation plot for real-frequency runs (November 2026)

`real_freqs_plots_from_file`: new plot `optimal`. For a unit-norm incident field E the
extinction is K Im(E^H A^-1 E) (the normalisation of `dda_extinction`). One eigenmode
gives K Im(1/lambda) (the `C.extsall` curves); the best of all fields gives K times the
largest eigenvalue of (A^-1 - A^-H)/2i. Because A is complex symmetric, not normal, the
modes are biorthogonal rather than orthogonal, and a mixture of modes can beat the best
single mode. Left tile: this run's extinction, the best single mode and the optimum vs
wavelength. Right tile: their ratio (1 = one mode is optimal) and the overlap
|u^H v_best| of the optimal field with the best mode. The plot rebuilds A with the
solver's own builders at `Input.plotting.optimal.n_points` evenly spaced wavelengths
(default 40), plus the peaks and the coef_map wavelengths (or `optimal.wls`); the
console lists ratio and overlap at those. Runs without saved eigenvalues (lean, or
the direct solver) compute them at the same points. Slow for large arrays (one dense
inverse and Hermitian eig per point): turn it off in the Plotting tab or set fewer
points. Example (HG01 chain, 1x120): ratio 1.04 and overlap 0.98 at the strongest
mode, up to 1.25 away from it. theory.pdf has a new section "Biorthogonal modes and
optimal excitation".

---

# Update: coefficient map (incident field on the modes) for real-frequency runs (November 2026)

`real_freqs_plots_from_file`: new plot `coef_map` (eig runs with mode tracking). |c| of
every labelled mode, c = V \ Einc (stored as `C.coefs_ordered`), at chosen wavelengths
(`Input.plotting.coef_map.wls`, nm; [] = extinction peak; several values = several
tiles): an (m_x, m_y) map for 2D arrays, |c| vs m for line arrays and cylinders.
Each tile names the strongest mode and its share |c_max|^2 / sum |c|^2 (100 % = no
other mode excited); the overall title gives the incident field (polarization,
angles, plane wave or Gaussian w0, z, HG order).
Polarizations/branches of a label are combined as sqrt(sum |c|^2). The labels are those
of the tracking (`Input.track.basis`, e.g. "sine" or "k"). Plotting tab: "Show
coefficient map" and its wavelengths. Re-plot a saved file with
`real_freqs_plots_from_file(file, struct('coef_map', struct('wls', [1060 1065])))`.

---

# Fix: half-period shifts split the array into equal parts (November 2026)

`DDA_codes/Make_generic_Array.m`: the shifts (`Input.array_shifts_x/y`, all lattice
types that support them) started at index ceil(n*l/(shifts+1)), so one shift of a
120-particle line gave parts of 59 and 61 and put the step one period off the centre.
One helper, `apply_halfshifts`, now splits the n positions into shifts+1 parts as
equal as possible (120: 60 + 60; 121: 61 + 60; two shifts of 120: 40 + 40 + 40).
The x shifts also used p(2)/2 instead of p(1)/2 (wrong for non-square periods).
Half-shifted geometries change; together with the new file names, old results are
not reused.

---

# Update: compact result names; one field profile, Gaussian, corrected (November 2026)

| File | |
|---|---|
| `Help_functions/construct_filename.m` | rewritten: compact names (below) |
| `DDA_codes/Make_3D_Field_wl_looped.m` | presets reduced to "none" and "Gaussian"; Gaussian corrected |
| `DDA_Master.m`, `Help_functions/DDA_parameter_GUI.m` | the old preset options removed |
| `Help_functions/dda_update_derived_inputs.m` | old preset names in saved defaults/favorites: gaussian_proper -> Gaussian, others -> none (printed) |
| `DDA_codes/Calculate_3D_DDA_generic.m` | no `Input.Gaussian.bool` |
| `README.md`, `README_results_index.md` | name format |

**Names**: fields separated by single spaces, no spaces or units inside them
(lengths nm, energies eV; the beam waist in um), e.g.
`1Dx 1x120 p(700,700) m[0000] PBCx700n10 HG10 r65 AuDM-BK7 L[0,0,90]Gau[80um,HG01] eig wl688.80-1239.84(N1500).mat`
(113 characters; before 210). Mirrors `m[1100]` (or `m[1,-1,0,0]`), periodic boundary
`PBCx700n10` (direction, cell, image cells per side), half-period shifts `HG10` (only
when there are any; `HG10b` for 1Dsquare_halfshifts2), particle-substrate `AuDM-BK7`
(AuDM, AuM, glM<n>, AuDMLWA; `n1.4665` for a constant host), field `L[theta,phi,psi]`
(L, RC, LC) plus `Gau[w0um,(z..um,)HGxy]`. Results saved under the old names are not
found under the new ones and are recalculated; `find_results` still finds them (it
reads the inputs in the files).

**Field profiles**: `Input.preset_field` is "none" (plane wave) or "Gaussian" (was
gaussian_proper, still accepted). sine_y, exp_y, gaussian_y, inverse_gaussian_y,
hermite_gaussian, hermite_JTK_tests and nullfield are removed, with their options
(Hermite_scale/offset/R, preset_sine_*). The PC mirror test is selected by its array
type only. Corrections to the Gaussian:
* centred on the array: the centre was half a period off (an HG01 node fell on
  particle 60 of 120 instead of between 60 and 61);
* evaluated at the actual particle positions (it used index x period, so half-period
  shifts, edited layouts and any other particle order were ignored);
* cylinders: a 1D beam in x with amplitude sqrt(w0/w) (Hermite_order_y not used);
* Hermite polynomials by recursion (the File Exchange helper is no longer needed).
Unchanged: waist w0, z_R = k w0^2/2 in the host, w(z), curvature exp(-i k r^2/2R)
(Gaussian_z > 0: waist beyond the array), no Gouy phase (constant).

---

# Fix: result files with long names could not be saved on Windows (November 2026)

`Unable to write to file ... because it appears to be corrupt` came from the path length:
folder + file name over 259 characters (e.g. a Gaussian/HG incident field adds ~45
characters to the name), which MATLAB's -v7.3 save cannot handle on Windows.
`construct_filename` now shortens names that would be too long, deterministically:
the beginning of the name, then ` ... #<fingerprint of the full name>.mat`. Re-runs
find the same file; names that fit are unchanged. `save_real_freq_results` prints
(instead of a suppressed warning) if a path is still too long.

---

# Update: k labels for real-frequency runs, wider GUI, no mode absorption ratio (November 2026)

| File | |
|---|---|
| `DDA_codes/dda_mode_basis.m` | modified: "k" labels (plane-wave coupling) shared by tracking and complex plots |
| `DDA_codes/Track_eigenmodes.m` | modified: `Input.track.basis = "k"`, `Input.track.k_step` |
| `DDA_codes/complex_freqs_plots_from_file.m` | modified: uses the shared k code; mode absorption ratio removed |
| `DDA_codes/Make_int_matrix_A_SI.m`, `DDA_codes/calc_alpha.m` | modified: `Input.proj_real_freq` removed (only the absorption ratio used it) |
| `DDA_codes/Calculate_mode_absorption_ratio.m` | deleted |
| `Help_functions/DDA_parameter_GUI.m` | modified: window 1.1x wider; "k" basis and k_step on the Tracking tab |
| `DDA_Master.m` | modified: `Input.track.k_step` |

**Real-frequency runs** can label modes like the complex plots' `labels = "k"`:
`Input.track.basis = "k"` (Tracking tab). The score of a label m is the strongest
plane-wave coupling |v.' e_k|^2 among the k values nearest to its calibration k_m,
so the best-scoring label is the k label of the complex plots (identical for line
arrays; 2D arrays use a coarser k grid, 4 instead of 10 points per pi/L). Labels
are assigned with the same matching as before ("basis", "overlap", "hybrid").
`Input.track.k_step = true` undoes the pi steps of half-shifted lattices first.
Along periodic directions with more than one cell, "k" falls back to "auto".

**Mode absorption ratio** at complex frequencies (tile titles of 2D-grid scalar modes)
removed with its helper and `Input.proj_real_freq`.

---

# Update: cleaned for the public repository (November 2026)

* No personal paths, host names or project numbers in the code: `Input.cluster.project`
  is "" (set yours in the GUI and press "Save defaults"; `DDA_saved_defaults.m` is
  git-ignored). `dda_results_folder` no longer falls back to a personal folder,
  `calc_alpha`/`Mie_polarizability_DDA` read `Johnson_gold.csv` from `DDA_codes`,
  `Make_DDA_save_folder_and_name(_TS)` save to `DDA/Data` instead of per-host folders.
* `Get_cpu_list_and_folder_hierarchy_TS.m` (machine names, local paths) is optional
  and git-ignored; without it the DDA folder is the root folder.
* `Make_checkerboard_array.m` (hexagonal arrays) and `Random_NP_sample_design_DDA.m`
  (square_random) moved from `DDA_codes/Archive` to `DDA_codes`; `Archive/` is
  git-ignored (only the legacy validation solvers need it).
* `.gitignore`: the above, `Data/`, and `theory.pdf` kept despite `*.pdf`.

---

# Update: plane-wave (k) mode labels and the 1D mode plot (November 2026)

| File | |
|---|---|
| `DDA_codes/complex_freqs_plots_from_file.m` | modified: `labels = "k"`, `k_step`, new 1D mode tiles |
| `Help_functions/DDA_parameter_GUI.m` | modified: "k" labels and the k_step option |
| `DDA_Master.m` | modified: `Input.plotting.complex.k_step` |
| `theory.pdf` | updated: section 4 (plane-wave labels, HG arrays), node-method recommendation |

**`labels = "k"`**: the mode number of the plane wave that couples most strongly to
the QNM, |v.' e_k|^2 summed over the signs of k, mapped to a mode number with the
standing waves sin((m+1) pi u) on the same positions as calibration. Confidence
1 - (second peak)/(first); below `min_confidence` the title shows the runner-up, (12|10).

**Half-shifted (HG01) chains**: each mode is a standing wave times the pi step, so a
plane wave couples to it at two |k| (k_m +- ~1.2 pi/L) and not at k_m itself. The
plane-wave labels give 1,0,3,2,5,4,7,6,9,8 going down from the centre of the branch,
then become ambiguous (equal lobes). `k_step = true` undoes the steps first (sign flip
at every half-period gap): one peak per mode, labels 0,1,2,3,... (the step-matched,
HG-shaped field). The node count measures neither and is fragile for these modes
(20-50 % standing wave): it gave pairs in one computation and an ascending sequence
in another.

**1D mode tiles** (line arrays, cylinders): |p| grey dashed and Re p solid on the left
axis, cumulative (unwrapped) phase/pi on the right axis. The mode is scaled to
max |p| = 1 and given the global phase 1/2 arg(v.'v), which makes it as real as
possible, so Re p and the phase no longer depend on the arbitrary eigenvector phase.
cyl_TE: p_z dotted (|p_z|) and dash-dotted (Re p_z).

---

# Update: sine-basis mode labels in the complex-frequency plots (November 2026)

| File | |
|---|---|
| `DDA_codes/dda_mode_basis.m` | new: the sine/Fourier position basis and basis scores, shared |
| `DDA_codes/Track_eigenmodes.m` | modified: uses dda_mode_basis (same labels as before) |
| `DDA_codes/complex_freqs_plots_from_file.m` | modified: sine and node-count labels, confidence, label map |
| `Help_functions/DDA_parameter_GUI.m` | modified: the new complex plotting options |
| `DDA_Master.m` | modified: the new complex plotting options |

The QNMs (local minima of |lambda|) are labelled by projecting their eigenvector on
the same basis the real-frequency mode tracking uses, so a mode has the same label
in both kinds of run. `Input.plotting.complex.*`:

* `labels = "sine"` (default): the basis function with the largest share of the
  mode's power; the share is shown as the confidence (tile titles), and labels below
  `min_confidence` (0.25) get a grey '?'. `"nodes"`: sign changes of the mode along
  x and y (made real by a global phase); confidence = how standing-wave-like it is.
  `"fft"`: the old FFT peak (`plot_fft`, `nufft` apply only to it).
* `basis = "auto"`: Fourier along periodic directions, sine otherwise.
* `label_map = true`: the label of the smallest-|lambda| eigenvector at every grid
  point, with the QNMs marked.

Tested on QNMs of a 1Dx 1 x 120 column (Au_Drude_Mie, 700 nm): uniform column, FFT,
sine and node count agree on 15/15 modes. Half-shifted column (array_shifts_y = 1):
the modes are not sine-like; FFT mislabels the odd modes (20/35 match the node
count), sine labels number the modes monotonically in energy but call every
symmetric mode nodes + 2 (20/35); use `labels = "nodes"` there.

---

# Update: mirror fix, line arrays and cylinders in the real-frequency plots (November 2026)

| File | |
|---|---|
| `DDA_codes/Make_int_matrix_A_SI.m` | modified: mirror fix; built in double precision; periodic cell size checked |
| `Help_functions/DDA_parameter_GUI.m` | modified: tip for "Periodic image cells" |
| `DDA_codes/Calculate_DDA_eigenmodes.m` | modified: `C.solver_version = 4` |
| `DDA_codes/Calculate_3D_DDA_generic.m` | modified: `C.solver_version = 4` |
| `DDA_Master.m` | modified: old mirror results are recalculated; note when a mirror has nothing to reflect |
| `DDA_codes/real_freqs_plots_from_file.m` | modified: line arrays and cylinders, band plot, mode plots, negative-extinction report |
| `Help_functions/calc_Q_varQ_DDA.m` | modified: own figure for the fits (the crash), Qprom |
| `DDA_codes/Lorentzian_polarizability_DDA.m` | modified: the passivity message is printed (warnings are off) |

**Mirrors (bug, all arrays with `Input.mirror_A`).** The planes go through the
outermost particles (x = 0, x = max(x), ...). The particles on a plane are their
own images, but their images were added again: their couplings were counted twice,
three times for a 1 x N column between the two x mirrors. The matrix was then not
passive: transmission > 1 and quasi-normal modes with Im(omega) > 0 (e.g. 1Dx,
1 x 120, mirror [1,1,0,0], Au_Drude_Mie: a mode at 1.171 eV growing at +1 meV).
Now these images are left out, and the dipole components that the symmetry makes
odd on the plane (s*M = -1) are decoupled. Checked against full symmetric arrays
(half array + mirror = full array to 1e-16 in 1Dx, 2D and 3D, s = +1 and -1); the
old code differed by about 0.5 % for every mirrored array. Mirror results without
`C.solver_version >= 4` are treated as outdated and recalculated.
A mirror through every particle (x mirrors on a 1 x N column) changes nothing;
DDA_Master prints a note.

**Interaction matrix in double precision.** `Array.r_i` is stored single, and
`Make_int_matrix_A_SI` built the whole matrix `'like'` it, so the eigen solver
and A\Einc ran in single precision (the cylinder builders were already double).
It now converts the positions to double. Eigen solves take somewhat longer.

**Periodic cell size is checked.** A cell given in nm as if in metres (700e-9
in the GUI's nm field = 7e-16 m) put every periodic image on top of a particle.
In single precision the enormous self-image terms overflowed to Inf and were
silently set to 1e-10, so the run gave the same non-passive result as the mirror
bug (transmission > 1). Now a cell not longer than the array in that direction
stops the run with an error that names the units. `periodic_Ncells` is the number
of image cells on EACH side: 0 = none, 1 = 3 cells in all, 10 = 21 cells.

**Line arrays** (1 x Ny, Nx x 1) and **cylinders** (cyl_TM, cyl_TE, always a line along x):
the driven-response and eigenmode figures have two panels, the transmission and
|p| vs wavelength and position along the line. The band plot shows only the
direction(s) with more than one particle. The incident-field plot uses line plots
and a 1D angular spectrum. In `Input.plotting.dipoles.modes`, (m,0) and (0,m) both
mean mode m along the line. 2D arrays are plotted as before.

**Crash** "The layout does not have sufficient space": the Q fit drew into
`figure(10)`, which could be a full 2 x 3 dipole figure. It now has its own figure
("Q fit"); `calc_Q_varQ_DDA(..., plotQ)` also takes a tiledlayout or axes handle.
A failing fit no longer stops the other plots. Qprom used a fixed 1060 nm and had
the wrong sign.

**Eigenmode |p| was zero**: `C.P_eigs` is single, and |p|^2 ~ 1e-62 underflows in
single precision. The mode plots now work in double.

**Transmission exp(-ext) > 1** (negative extinction) cannot happen with a passive
system, for any incident field. The plots now print what in the run breaks
passivity: particles above the radiation-reaction limit (e.g. a Lorentzian with
`lorentz_correction = "none"` and a large A0), `Input.Hermitian`, legacy cylinder
polarizability, gain, mirror results from before the fix, periodic boundaries.

---

# Merge: layered background / cyl_TE with the updates from the main chat (November 2026)

This package combines the layered-background and cyl_TE update with the changes
made in parallel: mirror reflection in 2D/3D, Input.P_method, one SI matrix for
the direct solvers, mldivide/linsolve, nm real axis for the complex plots, wall-time
handling and the array editor. Copy the files into your DDA folder, keeping the
folder structure.

After copying, delete the old duplicates
`Help_functions/complex_freqs_plots_from_file.m` and
`Help_functions/real_freqs_plots_from_file.m`: the current versions are in
`DDA_codes`, and MATLAB uses whichever copy comes first on the path.

# Update: measured time estimate with layers (November 2026)

| File | |
|---|---|
| `Help_functions/dda_layer_benchmark.m` | new: times the run's own matrix build with the layered background (cached) |
| `Help_functions/dda_time_estimate.m` | modified: uses the benchmark instead of the fixed per-cylinder guess |
| `Help_functions/DDA_parameter_GUI.m` | modified: shows the layer part of the estimate; option "Benchmark the layer cost" |
| `Help_functions/collect_dda_jobs.m` | modified: notes the layer part when comparing measured and estimated times |

The build is timed single-threaded at the corners of the frequency range
(the most negative imaginary part, where the integration path is longest),
for the run's own array and layer stack. For cluster runs it is scaled by
the speed of this computer relative to Roihu (an eig benchmark, once per
session). Options: `Input.layers.benchmark` (false: rough model),
`Input.layers.benchmark_budget` (s, default 5), `Input.layers.speed_factor`
(fixed Roihu/local factor).

---

# Update: layered background and cyl_TE (November 2026)

Copy the files into your DDA folder, keeping the folder structure. Modified
files replace yours; new files are added.

## What is new

**Layered background for the cylinder models** (`Input.layers`, see `dda_layers.m`).
The cylinders (along y, array along x) sit in a host layer; any number of layers
can be put below and above it, stacked along z:

    Input.n = 1.4665; Input.substrate = 'const';   % host layer (the film)
    Input.sample.type = "glass_Mie"; Input.Mie_glass_n = 1.33;   % the cylinders
    Input.layers.enable         = true;
    Input.layers.host_thickness = 1.18e-6;      % film thickness
    Input.layers.z_cyl          = 60e-9;        % cylinder axes above the film's lower interface
    Input.layers.below_n = [1.33 1.52];  Input.layers.below_d = 2e-6;   % 2 um of 1.33, then glass
    Input.layers.above_n = 1;            Input.layers.above_d = [];     % air

Indices are listed from the host layer outwards; the last one of each list is
the half-space. The same options are on the Physics tab of the parameter GUI.

**cyl_TE** (`Input.dimensionality = "cyl_TE"`): the same cylinders with the
field in the xz plane, two unknowns (p_x, p_z) per cylinder. In a film, cyl_TM
couples to the TE modes of the film and cyl_TE to its TM modes.

**Corrected cyl_TM polarizability.** The expression in `Mie_polarizability_DDA.m`
was too small: by the factor 1/(n_cyl + n_host) for thin cylinders (0.36 for
1.466 in 1.33), 0.51 at r = 101 nm. It is replaced by the exact Mie coefficient
(alpha -> pi a^2 eps0 (eps_cyl - eps_host) for thin cylinders, and
Im(1/alpha) = -k^2/(4 eps0 n^2) for lossless ones). Old cyl_TM result files are
treated as outdated and recalculated. `Input.particle.cyl_alpha_legacy = true`
gives the old values.

## Limits

- cyl_TE is solved by the eigen solver only (`inversion_tool = 'eig'`).
- The layered background does not work with `fft_gmres`.
- All cylinders must lie inside the host layer (0 < z_cyl + z < host_thickness);
  each is still a single line dipole, also next to an interface.
- Layer indices are constants (they may be complex). Very high indices make the
  Sommerfeld integral slow.
- With layers the incident plane wave comes from the top half-space; the angle
  is the vacuum angle, as before.
- Mirrors and periodic boundaries are not available for cylinders (as before).

## Tested

In Octave: `dda_layer_selftest` (run it once in MATLAB), an independent
numerical reference for the Green's function at a complex frequency, and
complex-frequency scans of a 184-cylinder chain in a film with both
polarizations. Not run at all (Octave cannot execute them): the edits to
`DDA_Master.m`, `Mie_polarizability_DDA.m`, the GUI, `construct_filename.m`,
`Make_3D_Field_wl_looped.m`, the eigen solver dispatch and the plot files.
They pass a syntax check only.

## Files

| File | |
|---|---|
| `DDA_Master.m` | modified: `Input.layers.*`, `cyl_TE`, `particle.cyl_alpha_legacy`, outdated cyl_TM results |
| `DDA_codes/dda_layers.m` | new: reads and checks `Input.layers` |
| `DDA_codes/dda_layer_green.m` | new: reflected part of the 2D Green's function of a layer stack |
| `DDA_codes/dda_layer_field.m` | new: plane-wave field inside the host layer |
| `DDA_codes/dda_layer_selftest.m` | new: checks of the three functions above |
| `DDA_codes/Make_cyl_TE_matrix_A_SI.m` | new: 2N x 2N matrix for cyl_TE |
| `DDA_codes/Make_cyl_TM_matrix_A_SI.m` | modified: layered background |
| `DDA_codes/Mie_polarizability_DDA.m` | modified: exact cylinder polarizabilities (TM corrected, TE new) |
| `DDA_codes/Calculate_DDA_eigenmodes.m` | modified: cyl_TE, `C.solver_version = 3` |
| `DDA_codes/Make_3D_Field_wl_looped.m` | modified: cyl_TE field, field inside the layers |
| `DDA_codes/dda_extinction.m` | modified: cyl_TE |
| `DDA_codes/Lorentzian_polarizability_DDA.m` | modified: error message covers cyl_TE |
| `DDA_codes/Calculate_3D_DDA_generic.m` | modified: clear errors for cyl_TE and for layers with fft_gmres |
| `DDA_codes/complex_freqs_plots_from_file.m` | modified: cyl_TE modes (|p_x| solid, |p_z| dashed) |
| `DDA_codes/real_freqs_plots_from_file.m` | modified: cyl_TE dipoles |
| `Help_functions/complex_freqs_plots_from_file.m` | modified: same change in the older copy |
| `Help_functions/real_freqs_plots_from_file.m` | modified: same change in the older copy |
| `Help_functions/DDA_parameter_GUI.m` | modified: cyl_TE and the layer options (Physics tab) |
| `Help_functions/construct_filename.m` | modified: cylinder materials and the layer stack in the file name |
| `Help_functions/dda_result_size.m` | modified: cyl_TE |
| `Help_functions/dda_time_estimate.m` | modified: cyl_TE, rough cost of the layers |

---

# What is in this update

Copy the contents of this `DDA` folder into your DDA folder (the one with
`DDA_Master.m`), keeping the folder structure. Modified files replace yours;
new files are added. Nothing needs to be deleted.

Not included on purpose: `DDA_saved_defaults.m` (written by the GUI's
"Save defaults" on your computer).

## Main folder

| File | |
|---|---|
| `DDA_Master.m` | modified: your latest version plus the GUI hook, saved defaults, size options |
| `README.md` | new: running locally, in the background and on the cluster; the GUI |
| `README_results_index.md` | new: finding and plotting saved results |

## dda_codes

| File | |
|---|---|
| `Calculate_DDA_eigenmodes.m` | new: one eigen solver for all dimensionalities (lean complex results) |
| `Track_eigenmodes.m` | new: mode labels and tracking |
| `compare_mode_tracking.m` | new: compare legacy and new mode labels |
| `real_freqs_plots_from_file.m` | new: real-frequency plots from a saved file (replaces JTK_plots) |
| `Lorentzian_polarizability_DDA.m` | new: MJH Lorentzian particles in SI units |
| `dda_solve_job.m` | new: everything one run computes, locally or in a job |
| `dda_extinction.m` | new: one extinction formula for all solvers |
| `dda_alpha_is_SI.m` | new: whether the polarizabilities are in SI units |
| `dda_alpha.m` | new: reads the polarizability (full or one per particle type) |
| `dda_compact_alpha.m` | new: one polarizability per particle type |
| `dda_expand_alpha.m` | new: full per-particle polarizability for the old solvers |
| `complex_freqs_plots_from_file.m` | modified: options from Input.plotting.complex, four bug fixes |
| `NLO_DDA_Make_3D_polarizability.m` | modified: Lorentzian option, compact polarizability |
| `Make_int_matrix_A_SI.m` | modified: matrix-inverse self term, compact polarizability |
| `Make_cyl_TM_matrix_A_SI.m` | modified: compact polarizability |
| `Make_1D_int_matrix_A_sub.m` | modified: SI check via dda_alpha_is_SI |
| `Make_3D_int_matrix_A_sub.m` | modified: SI check, missing n^2 |
| `Make_2D_int_matrix_A_sub_JTK.m` | modified: matrix-inverse self term, SI check |
| `Make_3D_int_matrix_A_sub_JTK.m` | modified: matrix-inverse self term, SI check |
| `Make_Afft_for_AxGain.m` | modified: SI check, missing n^2 |
| `Calculate_3D_DDA_generic.m` | modified: SI check in the extinction |

## Help_functions

| File | |
|---|---|
| `DDA_parameter_GUI.m` | new: the run-options window |
| `submit_dda_job.m` | new: submits background and cluster jobs |
| `collect_dda_jobs.m` | new: fetches and saves finished jobs |
| `dda_jobs.m` | new: lists, cancels and removes jobs |
| `dda_job_registry.m` | new: book-keeping of submitted jobs |
| `dda_job_futures.m` | new: keeps background jobs through `clear all` |
| `finalize_dda_result.m` | new: saves (and plots) a finished run |
| `save_real_freq_results.m` | new: saves real-frequency results |
| `dda_size_check.m` | new: stops runs that would not fit, before submitting |
| `dda_result_size.m` | new: size estimate used by the solver and the check |
| `plot_results.m` | new: plots rows found with find_results |
| `construct_filename.m` | new here (was missing from the repository): result file names |
| `dda_results_folder.m` | modified: real-frequency folder, root folder fix |
| `find_results.m` | modified: real-frequency results ('Kind') |
| `build_results_index.m` | modified: indexes run information of real-frequency results |
