# DDA codes: running simulations locally, in the background and on the cluster

`DDA_Master.m` runs a DDA simulation from the inputs in `Make_predefined_Inputs`
(inside `DDA_Master.m`). Every run is saved to a result file named after its
inputs, and running `DDA_Master` again with the same inputs plots that file
instead of recomputing.

This README covers **where the solve runs**: in your MATLAB session, in the
background on your own computer, or as a job on a cluster (CSC Roihu), and how
to collect results from jobs that ran without blocking MATLAB.

---

## 0. Choosing the run options in a window

With `Input.open_parameter_GUI = true` (top of `Make_predefined_Inputs`), `DDA_Master`
opens a window after the defaults are set. Tabs: **Array**, **Physics**, **Solver**,
**Tracking**, **Incident field**, **Plotting**.

* The window starts from the values in `Make_predefined_Inputs`; **Reset to defaults**
  returns to them. **Run** starts the run with the edited options and prints the ones
  that differ from the defaults; **Cancel** (or closing the window) stops without running.
* Options that only matter in a certain mode appear only when that mode is chosen:
  Lorentzian particle parameters, gmres parameters, thread/background/cluster settings,
  Hermite and Gaussian beam parameters, tracking settings (eigen solver only), and the
  settings of each plot once it is enabled.
* `inversion_tool` is on both the Physics and the Solver tab; changing one updates the other.
* **Figures as tabs:** with `Input.plotting.docked = true` (default) every plot figure opens
  as a tab in MATLAB's common figure window; `false` gives each figure its own window.
  (Since R2025a, giving a figure a size undocks it, so the DDA plots create sized figures
  through `dda_figure`, and docking also sets MATLAB's default window style for the session.)
* The Plotting tab shows the real-frequency plot options for real-frequency runs and the
  complex-frequency plot options (`Input.plotting.complex`) for complex-frequency runs.
* The Solver tab estimates the run time with `dda_time_estimate` (Help_functions):
  run time ≈ 6 s job start-up + 0.7 s per worker + frequencies per worker × t(Ndof),
  where t is the measured time of one frequency on one single-threaded worker on Roihu
  (1Dx 20 × 20: 0.09 s, 2D 20 × 20: 0.59 s, 3D 20 × 20: 2.8 s, 2D 30 × 30: 10.5 s,
  2D 40 × 40: 57 s), interpolated in the matrix size Ndof and growing as Ndof³ beyond
  the largest point (the measured slope there is 2.9). For cluster runs it is compared
  with the requested wall time. Queueing does not count against the wall time.
* `collect_dda_jobs` prints, for every collected job, the measured time per frequency
  per worker next to the model's value, and a ready-made row such as
  `[800 0.593]` for the table `dda_time_benchmarks` in `dda_time_estimate.m`.
* Vector fields take any MATLAB expression (`[700 700 100]`, `linspace(1,2,5)`). Units
  are given in the labels (nm, fs, um, eV, meV); the window converts to SI.
* Changing the range, `p` or `arrayN` recomputes what depends on them
  (`range.eV/wl/k/k_N`, `nmbofmodes`, `second_array_offset`, `Gaussian_sigma`).

Set `Input.open_parameter_GUI = false` for scripted runs.

### Editing the array layout

**Edit array layout...** (Array tab) opens a window that builds the array exactly as the run
will (`Make_generic_Array` with the current inputs) and shows every particle at its position
with its radius, coloured by z.

* **Select:** click a particle; Shift- or Ctrl-click adds or removes particles; drag on empty
  space to select everything inside a box (with Shift/Ctrl: add to the selection).
* **Move:** drag any selected particle to move the whole selection, or type values.
* **Values (nm):** x, y, z and radius apply to all selected particles; a field is empty when
  the selected particles differ. "Move selected by dx, dy, dz" shifts the selection.
* **Undo**, **Reset selected**, **Reset all particles**; overlapping particles are reported.
* **Done** stores the edits in `Input.array_edit`; **Cancel** leaves them unchanged. The Array
  tab then shows how many particles were moved and resized, with **Clear layout edits**.

How the edits are used:

* They are applied right after the array is built (`dda_apply_array_edit`), for local runs and
  cluster jobs alike, and the result file name gets `edited layout (n moved, m resized) #id`.
* They belong to the array they were made on. If `p`, `arrayN`, the lattice type or another
  array input changes afterwards, the GUI marks the edits as outdated and the run stops with
  a message instead of applying them to a different array.
* Edited radii (Mie particles) are computed once per distinct radius. Radius editing is not
  available for Lorentzian particles or for wall/Gaussian lattices, which set their own radii.
* FFT-based solvers (`fft_gmres`) need the regular lattice and do not accept edited layouts.
* **Save defaults** also saves the layout edits (in `DDA_saved_layout.mat`, next to
  `DDA_saved_defaults.m`), so the next run starts from the saved layout and the editor
  opens with it, ready to adjust. Saving without layout edits removes a saved layout.

### Favorites: saved input profiles

The **Favorites** tab of the GUI saves the current settings under a name (with an
optional note) and loads saved ones.

* A favorite stores every option of the window, in SI units, and the edited array layout,
  so loading it reproduces those settings exactly, even if the defaults in
  `Make_predefined_Inputs` change later.
* **Load** replaces the values in the window (check the tabs, then press **Run**);
  **Delete** removes the favorite. Saving under an existing name asks before replacing it.
* Favorites are `.mat` files in `DDA_favorites` next to `DDA_Master.m`, one per favorite.
  Decide whether that folder belongs in git (shared) or in `.gitignore` (personal).
* For scripted runs, set `Input.favorite = "name"` at the top of `Make_predefined_Inputs`:
  the favorite is applied after the saved defaults and before the GUI, if it opens.

Favorites and **Save defaults** are independent: saved defaults are applied on every run,
a favorite only when you load it or name it in `Input.favorite`.

### Saving the GUI settings as defaults

**Save defaults** writes every setting that differs from the values in
`Make_predefined_Inputs` to `DDA_saved_defaults.m`, next to `DDA_Master.m`, as plain
assignments in the GUI's units, for example:

```matlab
Input.p = [700 700 100]*1e-9;
Input.range.real_parts = linspace(1.1, 1.5, 40);
Input.parallel = "cluster";
```

Layout edits from the array editor are saved alongside, in `DDA_saved_layout.mat`, and the
generated file loads them.

`Make_predefined_Inputs` applies this file right after its own values whenever
`Input.use_saved_defaults = true`, also when the GUI is not opened, and prints which
settings it applied. The code in `DDA_Master.m` itself is never changed. Edit the file
by hand, or delete it to use only the values in `DDA_Master.m`.

* **Reset** goes back to the values the window opened with (including saved defaults).
* **Code values** goes back to the values written in `Make_predefined_Inputs`.
* Saving when nothing differs from the code values offers to delete the file.

If you keep the code in git, decide whether `DDA_saved_defaults.m` belongs in the
repository (shared defaults) or in `.gitignore` (personal ones).

---

## 1. Choosing where a run goes

In `Make_predefined_Inputs`, section **Computing options → Where the solve runs**:

```matlab
Input.parallel = "cluster";   % "none" | "threads" | "background" | "cluster"
```

| `Input.parallel` | Runs on | MATLAB while it runs | Good for |
|---|---|---|---|
| `"none"` | this MATLAB, one core | busy until done | small runs, debugging |
| `"threads"` | this MATLAB, local thread pool | busy until done | medium runs you want to see right away |
| `"background"` | local pool via `parfeval`, **one worker per run** | **free** | queuing several small/medium runs |
| `"cluster"` | a batch job (Roihu, or your own computer with `"Processes"`) | **free** (unless `Input.cluster.wait = true`) | large runs |

For one large run in the background on your own computer, using all cores, use
`"cluster"` with `Input.cluster.name = "Processes"`. The run gets its own pool,
so frequencies are computed in parallel, while MATLAB stays free.

`"none"` and `"threads"` save and plot as soon as the run is done. `"background"`
and `"cluster"` return immediately; results are fetched with `collect_dda_jobs`.

---

## 2. Everyday workflow for background and cluster runs

```matlab
DDA_Master                % submits the run and returns at once
% ... change inputs, submit more runs, plot old results, keep working ...
collect_dda_jobs          % fetches and saves every job that has finished
dda_jobs                  % lists submitted jobs and their state
```

* `collect_dda_jobs` takes no arguments. It checks all submitted jobs, saves the
  finished ones exactly as a blocking run would, deletes them from the cluster,
  and lists the ones still queued, running or failed. Run it whenever you like;
  nothing is collected automatically.
* For each collected job it prints how much was downloaded from the cluster and how
  long it took (the job files, which the MATLAB plugin copies when it first sees the job
  finished; if `dda_jobs` or `DDA_Master` saw it first, the copy happened then), how long
  reading the results into memory took, the size of the saved file, the time spent in the queue, the actual run
  time against the requested wall time, and what the run used (workers × threads,
  CPU or GPU, node, Slurm job id). For example:
  ```
  [1] finished 3D, 20 x 20, p=(700,700), ... (N=50) eV.mat
      fetching... downloaded 228.0 MB in 30.5 s (60 Mbit/s, during the state check), read 274.6 MB into memory in 1.4 s.
      saved: file 228.0 MB in 17.8 s.
      queued 2 min 10 s, ran 11 min 40 s of 30 min 00 s requested (39%), solve 10 min 55 s
      10 worker(s) x 4 thread(s), CPU on r1c05, Slurm job 4123456
  ```
  The same information is saved with the results in `C.run_info`. Running jobs
  show how long they have run against the limit, and failed jobs show their run
  time, which reveals jobs stopped by the wall-time limit. The run time is
  measured by MATLAB from start to finish of the job, so it includes starting
  the worker pool; Slurm's own accounting can differ by a few tens of seconds.
  The solvers use CPUs only; GPUs present on a node are reported as unused.
* After a job is collected, run `DDA_Master` with the same inputs to plot it, find it
  with `find_results` (see `README_results_index.md`), or plot the saved file directly:
  ```matlab
  real_freqs_plots_from_file(file)        % real-frequency runs
  complex_freqs_plots_from_file(name)     % complex-frequency runs
  ```
  `collect_dda_jobs` prints the saved file paths and also returns them.
* To plot automatically when a job is collected, set
  `Input.jobs.plot_on_collect = true` before submitting it.

### What happens when you run `DDA_Master` again with the same inputs

| Situation | What `DDA_Master` does |
|---|---|
| A valid result file exists | plots it (set `Input.recompute = true` to redo the run) |
| The run is queued, running or finished but not yet collected | nothing; tells you to run `collect_dda_jobs` |
| The previous job failed or was lost | submits the run again |

So a run is never submitted twice by accident, and a failed one can be retried
just by running `DDA_Master` again.

### Managing jobs

```matlab
dda_jobs                    % numbered list: state, where, when, result name
dda_jobs("cancel", 2)       % stop job 2 and forget it (also: [2 3] or "all")
dda_jobs("remove", 2)       % forget job 2, e.g. a failed one (deletes its cluster data)
```

States shown by `collect_dda_jobs` and `dda_jobs`:

| State | Meaning |
|---|---|
| `queued` | waiting for a worker or in the Slurm queue |
| `timeout` | ran past its wall time: Slurm stopped it (MATLAB still says running) |
| `running` | being computed |
| `finished` | done; `collect_dda_jobs` fetches and saves it |
| `FAILED` | stopped with an error; the error and the end of the job's output are printed. The job is kept until you remove it or submit the run again |
| `lost` | a background job from an earlier MATLAB session, or its pool was deleted; removed from the list |
| `missing` | the cluster job was deleted outside these tools; removed from the list |

---

## 3. Background runs (`Input.parallel = "background"`)

```matlab
Input.parallel               = "background";
Input.background.pool        = "Processes";  % "Processes" (recommended) or "Threads"
Input.background.num_workers = [];           % [] = MATLAB's default pool size
```

* Each run uses **one worker** and goes through its frequencies one after another.
  Several submitted runs are computed side by side. This is the right mode for
  queuing many small or medium runs, not for one very large run (see the
  `"Processes"` cluster option above).
* Background jobs exist **only in the current MATLAB session**. If MATLAB is
  restarted or the pool is deleted (`delete(gcp)`), unfinished jobs are reported
  as `lost` and have to be submitted again. Collect results before closing MATLAB.
* MATLAB allows only one pool at a time. If a pool is already open (for example a
  thread pool kept by a `"threads"` run), the background job uses that pool and
  prints a note. To get the pool type set in `Input.background.pool`, run
  `delete(gcp)` before submitting.
* A blocking `"threads"` run uses the same pool, so it has to wait for workers that
  background jobs occupy. Use `"none"` for quick local runs while background jobs
  are going.
* The legacy validation solvers (`Input.use_legacy_eig_solvers = true`) run only
  with `"none"` or `"threads"`.

---

## 4. Cluster runs (`Input.parallel = "cluster"`)

### One-time setup

1. Set up the MATLAB cluster profile for Roihu following CSC's MATLAB instructions
   (MATLAB: *Parallel → Create and Manage Clusters*). The profile name used here
   is `"Roihu R2026a"`.
2. In `Make_predefined_Inputs`, set your CSC project. **This is required on Roihu**;
   without it submission fails with "Required field AccountName is empty":
   ```matlab
   Input.cluster.project = "project_2001234";   % your CSC project (Slurm AccountName)
   ```
   `csc-projects` on Roihu lists your projects.
3. Check that the property names match your profile:
   ```matlab
   c = parcluster("Roihu R2026a"); c.AdditionalProperties
   ```
   The defaults are the names of the MathWorks Slurm plugin used on Roihu
   (`AccountName`, `WallTime`, `Partition`, `MemPerCPU`, `EmailAddress`,
   `AdditionalSubmitArgs`). If one differs, map it without editing code, e.g.
   `Input.cluster.property_names.project = 'ComputingProject';`

### Options (`Input.cluster.*`)

| Option | Default in `DDA_Master` | Meaning |
|---|---|---|
| `name` | `"Roihu R2026a"` | cluster profile; `"Processes"` = this computer |
| `project` | `""` | CSC project → Slurm account (**required on Roihu**) |
| `wall_time` | `"0:03:00"` | Slurm time limit (h:mm:ss); the job is killed after this. **Raise it for real runs** |
| `num_workers` | `10` | pool workers computing frequencies in parallel |
| `num_threads` | `4` | threads per worker (used by `eig`) |
| `partition` | `""` | Slurm partition; `""` = profile default |
| `mem_per_cpu` | `""` | e.g. `"4g"`; `""` = profile default |
| `email` | `""` | Slurm e-mail notifications |
| `submit_args` | `""` | any extra `sbatch` arguments |
| `additional_properties` | `struct()` | any other plugin property, copied as is |
| `property_names` | `struct()` | remap option names to your plugin's names |
| `wait` | `false` | `true` = block until done, then save and plot (old behaviour) |
| `delete_job` | `true` | delete the cluster job after collecting it (failed jobs are kept) |
| `auto_attach_files` | `true` | send the needed `.m` files with the job |
| `additional_paths` | `{}` | code folders that already exist on the cluster |
| `attached_files` | `{}` | extra files to send |
| `current_folder` | `"."` | working folder of the workers |
| `auto_add_client_path` | `false` | add your MATLAB path on the workers (always on for `"Processes"`) |
| `save_profile` | `false` | store these settings in the profile |

**Cores reserved** = (`num_workers` + 1) × `num_threads`. `batch` uses one extra
worker to run the job itself.

**A full Roihu CPU node** has 384 cores. Use the `medium` partition (it gives each job a
whole node) and fill the node: (`num_workers` + 1) × `num_threads` = 384, e.g. 95 workers ×
4 threads, 191 × 2 or 383 × 1, with `mem_per_cpu = ""` so the node's memory is available.
The GUI's **Use a full node** button sets this for the chosen threads per worker. Many
single-threaded workers suit small matrices; for large ones (Ndof above about 1000) use
4-8 threads per worker. Each worker needs a MATLAB Parallel Server licence, so very
large worker counts may exceed what is available.

**Memory.** The Solver tab shows the estimated size of the results, the memory of each
worker and of the job's controlling worker, and for cluster runs compares them with
`mem_per_cpu` × `num_threads` and the 2 GB `.mat` limit; `dda_size_check` applies the same
estimates before a job is submitted.

**Code version.** With `auto_attach_files = true`, the code is copied when the job
is submitted, so editing files while a job waits in the queue does not change it.
With `additional_paths`, the job runs whatever code is in those folders when it
starts.

**Choosing the wall time.** A run solves one eigenproblem of size
Ndof × Ndof per frequency (Ndof = 3 × particles in 3D). As a rough guide, a
20 × 20 3D array (Ndof = 1200) on a 50 × 50 complex frequency grid is 2,500
eigenproblems of a second or two each: roughly 5–15 minutes with 10 workers,
about an hour on a single worker. When unsure, give generous time: a job that
hits the limit is killed and its work is lost.

---

## 5. Where things are stored

| What | Where |
|---|---|
| Complex-frequency results | `DDA/Complex_data` |
| Real-frequency results | `DDA/Real_freq_results` (created on first use) |
| List of submitted jobs | `%LOCALAPPDATA%\DDA_jobs` (Windows) or `~/.local/share/DDA_jobs`: one small file per job, outside OneDrive so it is never synced |
| Cluster job data | the profile's job storage folders, locally and on the cluster, until the job is deleted |

MATLAB never deletes cluster job data by itself; `collect_dda_jobs` deletes a job
after saving it. On Roihu, files in `/scratch` that have not been accessed for
180 days are removed automatically, and scratch has a per-project quota (250 GiB
by default). Collect jobs well before that, and check where your profile's remote
job storage location points. These are CSC's policies at the time of writing;
see CSC's Roihu disk documentation for current values.

### Result size

The largest part of a real-frequency result is the per-mode dipole moments
`C.P_eigs`, about 8 · Ndof² · Nk bytes:

| 3D array | Ndof | Nk = 600 |
|---|---|---|
| 1 × 20 | 60 | ≈ 17 MB |
| 10 × 10 | 300 | ≈ 0.43 GB |
| 20 × 20 | 1200 | ≈ 6.9 GB |

When a job is collected, its output is downloaded, held in memory and saved again,
so your computer briefly needs about twice the result size in RAM. This is fine up
to about 1–2 GB per job. `Input.save.max_file_GB` (default 4) drops `C.P_eigs`
from very large result files when saving (with a warning); the mode dipole plots
are then skipped.

### Physics options worth knowing

* **Driven response P** (`Input.P_method`, eigen solver): `"direct"` (default) solves
  P = A\E, exact and robust near degenerate modes; `"eigenmodes"` builds the same P as a sum
  over eigenmodes, V·(c/λ), which shows how the modes make up the response but loses accuracy
  when the eigenvectors are nearly dependent.
* **Mirrors** (`Input.mirror_A`, planes at x = 0, x = max x, y = 0, y = max y): in 2D/3D,
  +1 makes the field mirror-symmetric (dipole components along the mirror even, the normal
  component odd; a magnetic-wall/PMC plane) and −1 antisymmetric (PEC-like). In 1D modes the
  one component satisfies p(−x) = ±p(x). The far mirrors lie on the outermost particles.
* **One interaction matrix for all solvers.** The eigen solver and the dense direct solvers
  (`mldivide`, `linsolve`, `gmres`, `gmres_gpu`, `inv`) build the same SI matrix, with mirrors and
  periodic boundaries, for 1D, 2D, 3D and cyl_TM. The FFT solvers (`fft_gmres…`) support neither
  mirrors nor periodic boundaries and stop if they are set.
* **Periodic boundaries** are summed directly over ±`periodic_Ncells` cells. This converges
  slowly at real frequencies (check by changing `periodic_Ncells`) and is not well defined at
  complex frequencies, where the far terms grow with distance.

### Size limits and how the code stays within them

* **Check before submitting.** `submit_dda_job` runs `dda_size_check` first. It measures
  the inputs, estimates the results, and compares them with the 2 GB per-variable limit
  of `.mat` files (used for job inputs, job outputs and complex-frequency result files)
  and, for cluster runs, with `mem_per_cpu` × `num_threads`. If something does not fit it
  stops with the reason and suggests how many frequencies fit in one run; otherwise it
  prints one line, e.g. `Sizes: inputs 29 MB, results ~389 MB (lean complex results), memory per worker 4.00 GB.`
* **One polarizability per particle type.** Identical particles share one tensor:
  `Array.alpha_types` (types × 3 × 3 × frequencies) with `Array.alpha_type` (the type of each
  particle) replaces the full `Array.alpha`. For 400 particles and 120,000 frequencies this
  is 17 MB instead of 6.9 GB. Read it with `dda_alpha`; the old direct and legacy solvers
  get the full array back from `dda_expand_alpha`.
* **Lean complex-frequency results.** Complex-frequency runs keep only what the complex
  plots use: the smallest eigenvalue, its eigenvector (single precision) and the
  degeneracy, about 8 × smaller. `Input.save.complex_all = true` keeps all eigenvalues,
  coefficients, `P` and the incident field. Runs with mode tracking keep everything.

---

## 6. Troubleshooting

| Message or symptom | Cause and fix |
|---|---|
| `Required field AccountName is empty in AdditionalProperties` | set `Input.cluster.project` (section 4) |
| `Failed to start mirror ... Unable to send channel-open request` | MATLAB's cached SSH connection died (CSC certificates expire after 24 h). Renew the certificate; the next submission then reconnects and retries automatically. To fix it by hand without restarting MATLAB: `dda_cluster_reconnect` |
| Warning: `"X" is not an existing AdditionalProperties entry` | your plugin uses another name; map it with `Input.cluster.property_names` |
| Cluster job `FAILED` with no error message | usually wall time or memory exceeded: raise `wall_time` or `mem_per_cpu` |
| Job shown as `TIMEOUT` | it ran past its wall time and Slurm stopped it. MATLAB itself keeps such a job "running" (the killed worker cannot report), so the job list compares the run time with the wall time (+2 min margin). Raise `wall_time` and run `DDA_Master` again; the job is replaced |
| Requested 10 s, got 1 min | Slurm counts whole minutes and rounds up (`0:00:10` → 1 min, `0:05:40` → 6 min); the GUI and the job list use the same rounding |
| Job `lost` | background job from an earlier MATLAB session or a deleted pool; run `DDA_Master` again to resubmit |
| A blocking `"threads"` run seems stuck | its pool is busy with background jobs; wait, or use `"none"` |
| `The legacy eigen solvers (validation) run only with ...` | legacy validation runs only with `"none"` or `"threads"` |
| `Complex frequencies need Input.inversion_tool = 'eig'` | complex-frequency runs use the eigen solver only |
| `This run does not fit (...)` | the run exceeds a size limit; split the frequency grid as suggested, or see "Size limits" above |
| `Unable to save data because data is too large` (from `batch`) | a job input over 2 GB; with the size check this is caught before submitting |
| `DDA_Master` says "Already submitted" | the run is pending; run `collect_dda_jobs`, or `dda_jobs("cancel", n)` to drop it |

To look at a failed cluster job in detail, list it with `dda_jobs`, then:

```matlab
c   = parcluster("Roihu R2026a");
job = findJob(c, 'ID', <job number shown by dda_jobs>);
getDebugLog(c, job)
```

---

## 7. Files involved

| File | Folder | Role |
|---|---|---|
| `DDA_Master.m` | main | inputs; `DDA_calc_looped` decides plot / submit / solve |
| `DDA_parameter_GUI.m` | `Help_functions` | the run-options window (`Input.open_parameter_GUI`) |
| `dda_favorites.m` | `Help_functions` | saved input profiles (Favorites tab, `Input.favorite`) |
| `DDA_array_editor.m` | `Help_functions` | the array layout editor (opened from the GUI) |
| `dda_apply_array_edit.m` | `Help_functions` | applies layout edits to the constructed array |
| `find_results.m`, `plot_results.m` | `Help_functions` | find saved results by their inputs and plot them (see `README_results_index.md`) |
| `submit_dda_job.m` | `Help_functions` | submits a background or cluster job, cluster options |
| `collect_dda_jobs.m` | `Help_functions` | fetches and saves finished jobs (no arguments) |
| `dda_jobs.m` | `Help_functions` | lists, cancels and removes jobs |
| `dda_job_registry.m` | `Help_functions` | book-keeping of submitted jobs (used internally) |
| `dda_cluster_reconnect.m` | `Help_functions` | renews a dead cluster connection without restarting MATLAB |
| `finalize_dda_result.m` | `Help_functions` | saves a finished run, optionally plots it |
| `dda_solve_job.m` | `dda_codes` | everything one run computes, wherever it runs |
| `Calculate_DDA_eigenmodes.m` | `dda_codes` | the eigen solver (`"none"` / `"threads"` inside each run) |
