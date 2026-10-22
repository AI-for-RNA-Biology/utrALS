# Pipeline: regenerating the results read by the figures

Plain shell scripts that regenerate, for each dataset used in the paper figures, the `data/` and `results/` folders
of the project, from the raw FASTQ files to the CLIP and RNA-modification analyses. Produced: the results read by the
notebooks in `notebooks/`, and the Welch profiles of every CLIP dataset and RNA modification used in the paper
(`clip_analysis`, `rnamod_analysis` and `m6a_analysis` steps, written to `.../analysis*/welch_profiles*/`). They run on any Linux
machine; SLURM is optional.

What each output contains, and the few results the pipeline does not regenerate, are described in
[`docs/data_provenance.md`](../docs/data_provenance.md).

| Dataset | Folder | Read by | Notes |
|---|---|---|---|
| MN_diff | `datasets/MN_diff/` | all notebooks | 31 steps. The `m6a_regions` step reads the candidate RBP lists written by `notebooks/Figure_regulators_lengthening.ipynb` |
| Astrocytes | `datasets/Astrocytes/` | `Figure_3UTR_CTRL_vs_VCP` | |
| VCP-KI-Ml240 | `datasets/VCP-KI-Ml240/` | `Figure_3UTR_CTRL_vs_VCP`, `Figure_global_regulators_NuclearWT_D35`, `Figure_regulators_lengthening` | |
| SFPQ_TDP_FUS_KD | `datasets/SFPQ_TDP_FUS_KD/` | `Figure_3UTR_CTRL_vs_VCP`, `Figure_regulators_lengthening` | kallisto file pairing corrected (see the provenance doc) |
| ENCODE-KD_K562, ENCODE-KD_HepG2 | `datasets/ENCODE-KD_*/` + `datasets/encode_kd.sh` | `Figure_3UTR_CTRL_vs_VCP` | raw processing kept in `ENCODE_DATA_DIR` |

## Layout

```
pipeline/
  config.sh                  tools and conda environments; data locations from ../paths.env (all overridable by
                             environment variables)
  lib.sh                     shared helpers (conda, config templating, strandedness, CLIP file selection, provenance);
                             puts ../src on PYTHONPATH, so the scripts import utrals without installing it
  runner.sh                  generic step runner used by every datasets/<dataset>/run_all.sh
  external_repos.tsv         library repositories whose scripts are run (apatools, bulkanalysis), pinned commits
  setup_external_repos.sh    clones them at those commits
  iclip_excluded.txt         the 3 iCLIP datasets not used
  steps/                     one script per unit of work (one sample, one CLIP file, ...), no scheduler code:
                             download_sra, star_align, kallisto_quant, featurecounts_utrs, aggregate_featurecounts,
                             genes_preprocessing, bulkapa, clip_integration, clip_analysis (Welch profiles with
                             scripts/clip/clip_analysis.py, layout CLIP_ANALYSIS_LAYOUT of the dataset's params.sh),
                             run_script (any script of ../scripts/ in a given environment)
  datasets/<dataset>/
    samples.tsv              sample sheet
    params.sh                dataset-specific settings (indices, strandedness, scripts, ...)
    configs/*.yaml           templates of the analysis configs (${RESULTS_DIR}-style placeholders)
    run_all.sh               runs the steps in order, locally
  datasets/encode_kd.sh      steps shared by the two ENCODE datasets (configs written per RBP)
  slurm/
    slurm.conf               account and resources per step
    submit.sh                submits the same tasks as SLURM job arrays, chained with dependencies
```

## Steps

| Dataset | Steps (`run_all.sh --list-steps`) |
|---|---|
| MN_diff | `download star kallisto featurecounts aggregate sum_utr_counts sum_gene_counts genes_preprocessing bulkapa clip_integration clip_stop_codon_integration clip_analysis clip_merge_welch clip_per_patient_integration clip_per_patient_merge_welch paraspeckles_prepare paraspeckles_run paraspeckles_merge rnamod_integration rnamod_stop_codon_integration rnamod_analysis rnamod_merge_welch m6a_regions m6a_regions_clip_list m6a_regions_clip_run m6a_regions_clip_merge m6a_regions_compare lncrna_prepare lncrna_features lncrna_run lncrna_merge` |
| Astrocytes | `download star kallisto featurecounts aggregate genes_preprocessing bulkapa` |
| VCP-KI-Ml240 | `download star kallisto featurecounts aggregate genes_preprocessing bulkapa clip_integration clip_stop_codon_integration clip_analysis m6a_integration m6a_stop_codon_integration m6a_analysis` |
| SFPQ_TDP_FUS_KD | `star kallisto featurecounts aggregate genes_preprocessing bulkapa clip_integration clip_stop_codon_integration clip_analysis m6a_integration m6a_stop_codon_integration m6a_analysis` (in-house data, no download) |
| ENCODE-KD_K562, ENCODE-KD_HepG2 | `download star kallisto featurecounts genes_preprocessing bulkapa` |

## Requirements

1. **Conda environments**, created from `envs/*.yml` (`conda env create -f envs/<name>.yml`):

   | Variable | Default | Used for |
   |---|---|---|
   | `ENV_PY` | `utrals` | Python analyses (bulkapa, genes_preprocessing, CLIP integration, MN_diff scripts) |
   | `ENV_BIO` | `bioinfo` | featureCounts, RSeQC; MN_diff m6A-switch scripts |
   | `ENV_COV` | `coverage_env` | samtools (BAM check after STAR) |
   | `ENV_SRA` | `sra_download` | fasterq-dump |

   Other environment names can be used by setting these variables; `""` uses tools already on your `PATH`.
   The `utrals` environment also needs the AI-for-RNA-Biology packages:
   `conda run -n utrals pip install --no-deps -r envs/lab_packages.txt` (see `envs/README.md`).
2. **STAR 2.7.10a and kallisto 0.46.1.** Taken from `PATH` if present, otherwise from the lab's `group/genomics/bin`
   (`STAR_BIN`, `KALLISTO_BIN`).
3. **Library repositories** `apatools` and `bulkanalysis`, whose command-line scripts the pipeline runs, at their
   `v0.2.0` commits (`external_repos.tsv`): `bash pipeline/setup_external_repos.sh` clones the missing ones into
   `REPOS_DIR`. Existing clones are not modified; it only reports those at another code version or with uncommitted
   changes (later commits that only update `poetry.lock` are accepted).
4. **Data locations:** results (`PROJECT_DIR`), raw FASTQ files (`DATABASE_DIR`), genome indices and annotations
   (`ANNOT_DIR`), CLIP (`CLIP_DIR`) and RMBase (`RMBASE_DIR`) files, library repositories (`REPOS_DIR`) and the ENCODE
   raw processing (`EXTERNAL_DATA_DIR`) are defined in `paths.env` at the repository root (local, git-ignored; create it
   with `cp paths.env.example paths.env`), relative to `DATA_ROOT` (the data storage). Edit `DATA_ROOT` there, or
   export any of these variables to override it.
5. **MN_diff only:** the m6A-switch steps read the candidate RBP lists
   `Tables/Figure_regulators_lengthening/regulators_of_nuclear_3UTR_lengthening_D22_to_D35_CTRL_<region>_APA_vs_nuclear_retention_class.csv`
   (written by `notebooks/Figure_regulators_lengthening.ipynb`).

## Running

```bash
# Where to write the results (default: PROJECT_DIR of paths.env)
export PROJECT_DIR=/path/to/ALS_Comp

# Everything, locally, 4 tasks in parallel
bash pipeline/datasets/Astrocytes/run_all.sh --jobs 4

# Part of it
bash pipeline/datasets/MN_diff/run_all.sh --list-steps
bash pipeline/datasets/MN_diff/run_all.sh --from genes_preprocessing --to bulkapa
bash pipeline/datasets/MN_diff/run_all.sh --only clip_integration,clip_merge_welch
ENCODE_RBPS="TARDBP FUS" bash pipeline/datasets/ENCODE-KD_K562/run_all.sh --from genes_preprocessing

# See what would run, without running it
bash pipeline/datasets/MN_diff/run_all.sh --dry-run

# On SLURM (from a login node): one job array per step, each waiting for the previous one
bash pipeline/slurm/submit.sh MN_diff --from star --dry-run   # print the sbatch commands
bash pipeline/slurm/submit.sh MN_diff --from star
```

- **Skipping finished work:** raw-data steps (download, STAR, kallisto, featureCounts) skip a sample whose output
  already exists. Set `FORCE=1` to recompute. Analysis steps always rerun.
- **Logs:** each run writes one log per task, plus `provenance.txt` (commits of this repo and of the library repos,
  paths and environments), into `$PROJECT_DIR/<dataset>/pipeline_logs/<date_time>/`.
- **Configs:** the analysis configs are rendered from `datasets/<dataset>/configs/*.yaml` and written next to their
  results (e.g. `results/3UTR_analysis/bulkapa.yaml`), as in the original runs.
- **SLURM chaining:** a step starts only when every task of the previous step succeeded. If one task fails
  (e.g. one CLIP file), fix it, then resubmit from that step.
- **Running one step by hand:** every script in `steps/` also works on its own, e.g.
  `bash pipeline/steps/featurecounts_utrs.sh SAMPLE BAM OUTDIR`. Each script's header gives its usage.

## Adding a dataset

Copy the closest dataset folder (e.g. `datasets/Astrocytes/`) and adapt:
- `samples.tsv`
- `params.sh` (paths and settings that differ)
- `configs/*.yaml` (copy the configs that produced the existing results and replace absolute paths with
  `${ANNOT_DIR}`, `${DATA_DIR}`, `${RESULTS_DIR}`, ...)
- the `STEPS` list and the `tasks_<step>` functions of `run_all.sh`

A step is a function `tasks_<step>` that prints one `task <script> <args>` line per independent unit of work; the
runner and the SLURM submitter both execute those lines.

## Verification

Run on 2026-10-07 into scratch `PROJECT_DIR`s and compared with the existing results in the project folder.
Not re-run anywhere: download, STAR and kallisto on all samples (hours per sample); their commands were checked
against the original scripts and the commands recorded in `STAR/*Log.out` and `kallisto/*/run_info.json`.
featureCounts, aggregation, bulkapa and the CLIP / m6A integration were also checked on a dataset not used in the
figures (TDP-43 / UPF1 knockdown): counts, shift lists and Welch statistics identical.

| Dataset | Step | Result |
|---|---|---|
| all | dry run | expected task counts; ENCODE: all 78 generated configs identical to the existing ones |
| SFPQ_TDP_FUS_KD | featureCounts (1 sample, strandedness inferred) | counts, summary and RSeQC report identical |
| SFPQ_TDP_FUS_KD | aggregate, genes_preprocessing, bulkapa (with the original kallisto output) | 212 files identical (incl. all 180 `shifts/`, 30 `genes_utrs/`); genes in a different row order; `table_of_pairs*.csv` differ by ≤ 9e-16 (rounding) |
| SFPQ_TDP_FUS_KD | kallisto with the corrected file pairing (TYZ4001A1, SLURM) | 86.4 % of reads pseudoaligned (original run: 7.5 %; STAR: 85 % uniquely mapped) |
| MN_diff | sum_utr_counts, genes_preprocessing (SLURM) | `df_counts_summed.csv` identical; genes in a different row order |
| MN_diff | bulkapa (SLURM) | 337 files identical (all current `shifts/` and `genes_utrs/` files); genes in a different row order; `table_of_pairs*.csv` differ by ≤ 5e-15 (rounding) |
| MN_diff | clip_merge_welch, clip_per_patient_merge_welch (SLURM, on the existing pickles) | `CLIP/analysis/welch_*.csv` and `CLIP_per_patient/analysis/welch_*.csv` identical |
| MN_diff | clip_per_patient_integration (1 eCLIP file, SLURM) | pickle identical |
| MN_diff | rnamod_merge_welch (2026-10-22, on the existing pickles) | `welch_*.csv` identical for m6A, m1A, m7G, Nm, Pseudo and RNA-editing. m5C differs: the existing `RNA_modifs/m5C/welch/analysis/welch_*.csv` (2026-08-24) predate the current m5C pickles (2026-09-18) and are labelled `hg38.m6A.gff_RMBase` (see the provenance doc) |
| MN_diff | m6a_regions_compare | all 15 output tables identical |
| VCP-KI-Ml240 | m6a_analysis (2026-10-22, on the existing m6A pickles, PAS and stop codon) | rendered config identical; `welch_profiles/welch_profile_hg38.m6A.gff_RMBase.pdf` renders identically to the existing one |
| MN_diff, VCP-KI-Ml240, SFPQ_TDP_FUS_KD | stop-codon integration and CLIP / RNA-modification analysis steps, dry run (2026-10-22) | task commands identical to the original pipeline's (the Pierre Klein m6A files are not included) |
| MN_diff, VCP-KI-Ml240, SFPQ_TDP_FUS_KD | `scripts/clip/clip_analysis.py` (one script, `--layout`) vs the three former per-dataset scripts (2026-10-22, on the existing pickles: a 10-CLIP subset with stop codon, and m6A / pseudouridine in full; MN_diff with and without `--welch-only`; the 3 knockdowns) | same output files; all 316 figures render identically (Welch profiles, by-region profiles, Fisher heatmaps and profiles) |
| MN_diff | scripts moved to `scripts/` and `utrals` (2026-10-22, on the existing results): `merge_welch` (m1A), `sum_utr_counts`, `paraspeckles_merge` (MALAT1), `m6a_regions_compare`, `m6a_regions_clip_merge`, `lncrna_merge`, switch group statistics | all outputs byte-identical to the existing ones |
| Astrocytes | genes_preprocessing, bulkapa (SLURM) | all 32 `shifts/`/`genes_utrs/` files and normalised counts identical; genes in a different row order; `table_of_pairs*.csv` differ by ≤ 5e-15 (rounding) |
| ENCODE K562 TARDBP | genes_preprocessing, bulkapa (SLURM) | `df_counts.csv` and all `shifts/`/`genes_utrs/` files identical; genes in a different row order; `table_of_pairs*.csv` differ by ≤ 5e-15 (rounding) |
| VCP-KI-Ml240 | genes_preprocessing, bulkapa (SLURM) | all 64 `shifts/`/`genes_utrs/` files and normalised counts identical; genes in a different row order; `table_of_pairs*.csv` differ by ≤ 9e-16 (rounding) |
