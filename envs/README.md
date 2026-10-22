# Conda environments

Environments used by the pipeline (`pipeline/config.sh`: `ENV_PY`, `ENV_BIO`, `ENV_COV`, `ENV_SRA`) and
by the notebooks:

| Environment | File | Used for |
|---|---|---|
| `utrals` | `utrals.yml` + `lab_packages.txt` | Python analyses: genes_preprocessing, bulkapa, CLIP integration, MN_diff scripts; kernel of every notebook in `notebooks/` except `Figure_m6A_lncRNA_switch` |
| `bioinfo` | `bioinfo.yml` | featureCounts, RSeQC `infer_experiment.py`; MN_diff m6A-switch scripts (`scripts/m6a_switch`); kernel of `Figure_m6A_lncRNA_switch` (needs pysam and ViennaRNA) |
| `coverage_env` | `coverage_env.yml` | samtools (BAM check after STAR) |
| `sra_download` | `sra_download.yml` | sra-tools `fasterq-dump` |

Each file lists only the packages the code of this repository (and of the lab packages it uses) imports or runs, at
the versions the results were produced with; conda and pip resolve their dependencies.

## Setup

From the repository root:
```bash
for e in utrals bioinfo coverage_env sra_download; do conda env create -f envs/$e.yml; done

# Packages of the AI-for-RNA-Biology GitHub organisation (apatools, bulkanalysis, cliptools, genomemanager,
# kallistomanager, omicstools, scanalysis), at their release v0.2.0
conda run -n utrals pip install --no-deps -r envs/lab_packages.txt

# This repository's utrals package (editable, without touching the environment's packages)
conda run -n utrals pip install --no-deps --no-build-isolation -e .
conda run -n bioinfo pip install --no-deps --no-build-isolation -e .
```

The lab packages are installed with `--no-deps`: `utrals.yml` provides their dependencies, and their `pyproject.toml`
files require each other from GitHub without a version, which pip cannot reconcile with the tags of
`lab_packages.txt`. Their command-line scripts (`apatools/scripts/bulkapa.py`, ...) are not installed: the pipeline
runs them from clones (`bash pipeline/setup_external_repos.sh`, see `pipeline/README.md`).

To run the notebooks in Jupyter, register the two kernels:
```bash
conda run -n utrals python -m ipykernel install --user --name utrals
conda run -n bioinfo python -m ipykernel install --user --name bioinfo
```

`Figure_3UTR_CTRL_vs_VCP` and the pipeline's genes_preprocessing / bulkapa steps call R through rpy2 (`limma`, in
`utrals`): when R is not found (e.g. a notebook run headless), set `R_HOME` to the `lib/R` folder of the `utrals`
environment.

Other environment names can be used for the pipeline by setting `ENV_PY`, `ENV_BIO`, `ENV_COV` and `ENV_SRA` (e.g.
`ENV_PY=my_env bash pipeline/datasets/...`).

The figures are byte-identical to the published ones when made with these package versions. Packages that conda or
pip resolve as dependencies (not pinned here) can differ between installations; this does not change the results, but
can change details of the rendered figures (e.g. font rasterisation in PNG files).
