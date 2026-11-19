# Finding and plotting saved results

Every DDA run is saved to a `.mat` file named after its inputs:

| Run | Folder | Plotted by |
|---|---|---|
| complex frequencies (`range.imag_parts` ≠ 0) | `DDA/Complex_data` | `complex_freqs_plots_from_file` |
| real frequencies | `DDA/Real_freq_results` | `real_freqs_plots_from_file` |

The file names are long and do not hold every input. To find results by any input
parameter, use `find_results`; to plot what you found, use `plot_results`.

---

## Quick start

```matlab
T = find_results('Kind', 'all', 'arrayN', [20 20 1], 'dimensionality', '3D');
T(:, {'kind','date','file','particle_r'})   % look at a few columns
plot_results(T, 1)                          % plot the newest match
```

```matlab
find_results('Kind', 'real')                % print a numbered list of all real-frequency results
```

---

## `find_results`

```matlab
T = find_results()                          % all complex-frequency results, newest first
T = find_results('Kind', 'real', ...)       % real-frequency results
T = find_results('Kind', 'all',  ...)       % both
T = find_results('Name', value, ...)        % only results that match every filter
find_results(...)                           % no output: print a numbered list
```

### Filters

A filter is an input name and a value. Nested inputs are written with `.` or `_`:
`'particle.r'` and `'particle_r'` are the same. Long arrays such as
`range.real_parts` are summarised as `range_real_parts_N`, `_min`, `_max`
(and `_imin`, `_imax` when complex).

| Value | Matches | Example |
|---|---|---|
| number or vector | equal within a relative tolerance (default 1e-6); vectors must have the same length | `'p', [700 700 100]*1e-9` |
| text | case-insensitive; `*` and `?` are wildcards | `'save_message', '*test*'` |
| cell array | any of the values | `'dimensionality', {'2D','3D'}` |
| function handle | your own test on the stored value | `'particle_r', @(r) r >= 60e-9` |

The columns `kind`, `file` and `date` can be filtered too, for example
`'date', @(d) d > datetime(2026,11,1)`.

### Options

| Option | Default | Meaning |
|---|---|---|
| `'Kind'` | `'complex'` | `'complex'`, `'real'` or `'all'` |
| `'Folder'` | the folder of `Kind` | another results folder (not with `'all'`) |
| `'Rebuild'` | `false` | re-read every file instead of using the cache |
| `'Tol'` | `1e-6` | relative tolerance for numbers |
| `'KeepEmpty'` | `false` | keep columns that are empty for all matches |

### The table

| Column | Content |
|---|---|
| `kind` | `'complex'` or `'real'` |
| `folder`, `file` | where the result is |
| `date`, `bytes` | file date and size |
| one column per input | e.g. `p`, `arrayN`, `particle_r`, `sample_type`, `cluster_num_workers`, `plotting_Q_enable` |
| `info_*` (real-frequency runs) | `info_created`, `info_solver`, and the run information: `info_run_solve_time_s`, `info_run_workers`, `info_run_threads`, `info_run_host`, `info_run_ran_s`, `info_run_queued_s`, ... |

The run information is stored for runs saved after this update; it comes from
`dda_solve_job` and, for cluster and background jobs, from `collect_dda_jobs`.

Examples:

```matlab
% Real-frequency 3D runs that took longer than 10 minutes to solve
T = find_results('Kind','real', 'dimensionality','3D', 'info_run_solve_time_s', @(t) t > 600);

% Complex-frequency runs of a 20 x 20 array with Lorentzian particles
T = find_results('arrayN',[20 20 1], 'sample_type','Lorentzian');

% All columns of the index
T = find_results('Kind','all', 'KeepEmpty', true);  T.Properties.VariableNames
```

If a name is not a column, the error message suggests similar column names.

---

## `plot_results`

```matlab
plot_results(T)              % row 1
plot_results(T, [1 3])       % rows 1 and 3
plot_results(T, 2, opts)     % with plot-option overrides
```

It picks the plot function from the `kind` column and opens the file from the
`folder` column. The plot options saved with each result are used unless you
override them:

* **Real frequencies:** `opts` is a plotting struct like `Input.plotting`, e.g.
  `struct('Q', struct('enable', false), 'emitted', struct('enable', true, 'wls', 1060))`.
* **Complex frequencies:** `opts` holds complex-plot options like
  `Input.plotting.complex`, e.g. `struct('quiver', false, 'full_range', false,
  're_lims', [1.12 1.18], 'im_lims', [25 30])`.

The same overrides work with the plot functions directly:

```matlab
real_freqs_plots_from_file(file, struct('dipoles', struct('surf_wls', [1053 1066])))
complex_freqs_plots_from_file(name, struct('modes', false))
```

### Complex-frequency plot options (`Input.plotting.complex`)

| Option | Default | Meaning |
|---|---|---|
| `real_axis` | `"eV"` | real axis in energy (`"eV"`) or vacuum wavelength (`"nm"`) |
| `full_range` | `true` | axes over the whole computed range |
| `re_lims`, `im_lims` | `[1.12 1.18]` eV, `[25 30]` meV | axes when `full_range` is false (`im_lims` is −Im) |
| `modes` | `true` | plot the modes at the minima of \|λ\| |
| `tiles_per_figure` | `20` | mode tiles per figure |
| `max_fft_order` | `1000` | skip modes with a higher spatial order |
| `quiver` | `true` | 2D/3D modes as arrows (needs `slanCM` on the path) |
| `imagesc_is_z` | `true` | 3D: the colour map shows p_z |
| `separate_pols` | `false` | 2D/3D: one tile per polarization |
| `surf_1D` | `true` | scalar modes as surfaces |
| `plot_fft` | `true` | show the spatial FFT used for the mode labels |
| `nufft` | `false` | non-uniform FFT for the labels |

They are set in `Make_predefined_Inputs` and in the parameter GUI (Plotting tab,
shown for complex-frequency runs), and saved with every complex-frequency result.
Files saved before these options existed use the defaults in the table.

---

## How the index works

* `build_results_index` reads only the `Input` variable of each file (and, in
  real-frequency files, the small `info` variable), never the results themselves,
  and flattens it to one row per file.
* The rows are cached in `results_index.mat` in each results folder. On later
  calls only new or changed files are read and deleted files are dropped, so the
  index always matches the folder. The first call after this update re-reads
  every file once (the index format changed).
* Inputs added after a file was saved are missing from its row. Give them the
  value old files effectively used in the `defaults` struct at the top of
  `build_results_index.m` (e.g. `periodic_boundary = [0 0]`), so filtering on them
  also finds older files.
* `'Rebuild', true` re-reads every file, e.g. if the cache is damaged.

---

## Files

| File | Folder | Role |
|---|---|---|
| `find_results.m` | `Help_functions` | search the results by their inputs |
| `plot_results.m` | `Help_functions` | plot rows of a `find_results` table |
| `build_results_index.m` | `Help_functions` | builds and caches the index (used by `find_results`) |
| `dda_results_folder.m` | `Help_functions` | the two results folders |
| `real_freqs_plots_from_file.m` | `dda_codes` | real-frequency plots from a file |
| `complex_freqs_plots_from_file.m` | `dda_codes` | complex-frequency plots from a file |
