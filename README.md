# utrALS

**A developmentally programmed RNA regulatory code integrates 3′UTR structure, m6A and RNA-binding proteins in motor neuron maturation and ALS**


The repository holds the paper analyses and figures.

---

## Reproducing the results

The processed data and results live outside this repository. Their locations are set in one local file, `paths.env`
(git-ignored), created from the template [`paths.env.example`](paths.env.example):
```bash
cp paths.env.example paths.env      # then set DATA_ROOT in paths.env
```
`DATA_ROOT` (the root of the data storage) is its only absolute path, and every other location (`PROJECT_DIR`, the
project folder holding the results; `ANNOT_DIR`; `CLIP_DIR`; ...) is relative to it. Adapt them in `paths.env` if your
folders are organised differently, or export any of these variables to override it (e.g. `PROJECT_DIR=/my/scratch`). The notebooks and
scripts read the locations through `utrals.paths`, the pipeline through `pipeline/config.sh`.

Set up the software once (details in [`envs/`](envs/README.md) and Software below): create the conda environments,
install the lab packages of [AI-for-RNA-Biology](https://github.com/AI-for-RNA-Biology) in `utrals`, and install this
repository's `utrals` package (`src/utrals/`, editable, so it follows the repository) in `utrals` and `bioinfo`:
```bash
for e in utrals bioinfo coverage_env sra_download; do conda env create -f envs/$e.yml; done
conda run -n utrals pip install --no-deps -r envs/lab_packages.txt
conda run -n utrals pip install --no-deps --no-build-isolation -e .
conda run -n bioinfo pip install --no-deps --no-build-isolation -e .
```
(The pipeline also finds the `utrals` package without installing it: it puts `src/` on `PYTHONPATH`.)

Two steps:

1. **Results** – [`pipeline/`](pipeline/README.md) regenerates them from the raw FASTQ files, locally or on SLURM, for the
   datasets used in the figures (MN_diff, Astrocytes, VCP-KI-Ml240, SFPQ_TDP_FUS_KD, ENCODE-KD_K562, ENCODE-KD_HepG2):
   ```bash
   bash pipeline/setup_external_repos.sh                       # apatools and bulkanalysis scripts, at v0.2.0
   PROJECT_DIR=/path/to/output bash pipeline/datasets/Astrocytes/run_all.sh --jobs 4
   ```
2. **Figures** – run the notebooks in [`notebooks/`](notebooks/).
   They write their figures to `Figures/<notebook>/` and their tables to `Tables/<notebook>/`.

### Software

Conda environments: [`envs/`](envs/README.md) (`utrals` for the pipeline's Python steps and the notebooks, `bioinfo`
for featureCounts and the m6A-switch scripts, `coverage_env` and `sra_download` for samtools and sra-tools).

Besides public packages, the code needs these packages of the
[AI-for-RNA-Biology](https://github.com/AI-for-RNA-Biology) GitHub organisation, all at their release `v0.2.0`, listed in
[`envs/lab_packages.txt`](envs/lab_packages.txt) and installed in the `utrals` environment with
`pip install --no-deps -r envs/lab_packages.txt` 

## Repository structure

Code is organised by role. The library is the `utrals` package, the command-line scripts run by the
pipeline are in `scripts/`, grouped by analysis stage, and what differs between datasets (samples, settings, configs)
is in `pipeline/datasets/<dataset>/`.

| Folder | Content |
|---|---|
| [`notebooks/`](notebooks/) | the 10 paper-figure notebooks, named by topic (see [`docs/notebooks_and_scripts.md`](docs/notebooks_and_scripts.md) §1) |
| [`src/utrals/`](src/utrals/) | the `utrals` package, imported by the notebooks and the scripts (`from utrals.figures import ...`): `paths` (data locations, `paths.env`), `figures` (`save_figure` / `save_table`, which keep `Figures/README.md` and `Tables/README.md` up to date; `fix_random_seeds`), `rbp_names`, `stats`, `string_api` and `string_analysis` (STRING queries, networks, enrichment), `enrichr`, `expression` and `go_expression` (expression over the time course), `welch` (Welch profile figures), `m6a_switch` (m6A structural-switch analysis) |
| [`scripts/`](scripts/) | command-line scripts run by the pipeline: `counts/` (summing sequencing runs), `clip/` (`clip_analysis.py`: Welch profiles of every dataset, layout chosen with `--layout`; Welch table merging; lncRNA cross-link density), `m6a_switch/` (m6A structural-switch pipeline) |
| [`pipeline/`](pipeline/README.md) | shell pipeline regenerating the results, one folder per dataset |
| `paths.env.example` | template of `paths.env` (local, git-ignored): locations of the data outside the repository (see above) |
| `pyproject.toml` | installs `src/utrals` (editable) |
| [`envs/`](envs/) | conda environments |
| [`cache/`](cache/) | archived STRING v12.5 responses, Enrichr library, QuickGO download and one reference figure, read by the notebooks so that a rerun gives the same figures without contacting these servers |
| `Figures/` | figures saved by the notebooks; `Figures/README.md` gives the legend of each one (git-ignored, rebuilt by the notebooks) |
| `Tables/` | tables (.csv / .xlsx) saved by the notebooks, one folder per notebook; `Tables/README.md` describes their columns (versioned, unlike `Figures/`) |
| [`docs/`](docs/) | [`data_provenance.md`](docs/data_provenance.md): how every result read by the notebooks was produced; [`notebooks_and_scripts.md`](docs/notebooks_and_scripts.md): which notebook makes which figure, what it reads, and what every module and script does |

### Paper-figure notebooks

- `Figure_3UTR_dynamics_over_time.ipynb` — nuclear vs cytoplasmic 3′UTR isoform usage over MN differentiation, CTRL vs VCP; RNA export
- `Figure_3UTR_CTRL_vs_VCP.ipynb` — 3′UTR shifts in VCP vs CTRL MNs and in independent datasets (VCP-KI ± ML240, SFPQ / TDP-43 / FUS KD, ENCODE KD, astrocytes); pathways
- `Figure_global_regulators_NuclearWT_D35.ipynb` — regulators of 3′UTR isoform choice selected at one time point (nuclear CTRL, D35)
- `Figure_regulators_lengthening.ipynb` — regulators of nuclear 3′UTR lengthening (D22 → D35, CTRL); writes the candidate RBP lists used by the m6A notebooks and the pipeline
- `Figure_RBP_SVD.ipynb` — SVD of per-patient RBP activity over differentiation
- `Figure_m6A_loading.ipynb`, `Figure_m6A_nucleotide_content.ipynb`, `Figure_m6A_RBP_correlations.ipynb`, `Figure_m6A_switch_regions_summary.ipynb`, `Figure_m6A_lncRNA_switch.ipynb` — m6A density, sequence content, m6A / RBP profile correlations and m6A structural switches
