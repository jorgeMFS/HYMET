# CAMI Benchmark Harness (`HYMET/bench`)

This directory houses the automation used to benchmark HYMET and a panel of external profilers on CAMI-style data. It covers database preparation, tool execution, evaluation against CAMI truth, aggregation, and figure generation. Everything runs with standard POSIX shell and Python utilities so the workflow is reproducible on bare metal or minimal containers, provided sufficient CPU, RAM, and disk are available.

---

## 1. Directory Overview

```
bench/
├── aggregate_metrics.py        # Merge evaluation outputs into TSVs
├── cami_manifest.tsv           # Sample manifest with local paths + optional URLs
├── convert/                    # Tool output → CAMI profile adapters
├── data/                       # Lightweight CAMI subsets (derived from sample_0)
├── db/                         # Database builders + cached indices
├── environment.yml             # Conda/micromamba environment specification
├── environment.lock.yml        # Frozen environment capture (optional)
├── fetch_cami.sh               # Helper to fetch manifest URLs
├── lib/                        # Shared shell helpers (common.sh, measure.sh, run_eval.sh, ...)
├── out/                        # Benchmark results (per sample/tool)
├── plot/                       # Matplotlib summary figure generator
├── refsets/                    # Shared reference FASTA subsets
└── run_*.sh                    # Tool-specific wrappers
```

Key support scripts:

| Path | Purpose |
|------|---------|
| `lib/common.sh` | Logging, path resolution, root discovery. Sourced by every wrapper. |
| `lib/measure.sh` | Wraps commands with `/usr/bin/time -v` and appends to `out/runtime_memory.tsv`. |
| `lib/run_eval.sh` | Invokes `HYMET/tools/eval_cami.py`; deletes contig reports when no usable pairs exist. |
| `convert/*.py` | Converts raw tool outputs into CAMI-compliant profiles. |
| `aggregate_metrics.py` | Builds `summary_per_tool_per_sample.tsv`, `leaderboard_by_rank.tsv`, and `contig_accuracy_per_tool.tsv` after pruning empty contig rows. |

---

## 2. Prerequisites

| Requirement | Notes |
|-------------|-------|
| OS | Linux/x86_64 with Bash ≥ 4 and Python ≥ 3.9. |
| Package manager | `micromamba` (preferred) or `conda` to materialise the environment. |
| CPU/RAM | Designed for ~16 hardware threads and ≥ 32 GB RAM. MetaPhlAn4 stays < 20 GB when using `--split_reads` and `METAPHLAN_THREADS=4`. |
| Disk | Expect ~160 GB total (MetaPhlAn DB ≈ 34 GB, HYMET references ≈ 45 GB, tool outputs ≈ 2 GB, plus working space). Removing HYMET minimap2 indices after runs keeps usage < 200 GB. |
| Taxonomy dump | `HYMET/taxonomy_files/` must contain `names.dmp`, `nodes.dmp`, etc. |
| MetaPhlAn database | Install with `metaphlan --install -x mpa_vJun23_CHOCOPhlAnSGB_202307 --db_dir bench/db/metaphlan`. |

---

## 3. Input Data

### 3.1 Manifest (`cami_manifest.tsv`)

Tab-separated file with one entry per sample:

```
sample_id  contigs_fa  truth_contigs_tsv  truth_profile_tsv  contigs_url  truth_contigs_url  truth_profile_url
cami_sample_0  /data/cami/sample_0.fna  /data/cami/sample_0/.../gsa_mapping_new.tsv  /data/cami/sample_0/taxonomic_profile_0.txt
cami_i_lc      data/cami_i_lc/contigs.fna      data/cami_i_lc/truth_contigs.tsv      data/cami_i_lc/truth_profile.tsv
...
```

Columns:
- `*_fa` / `*_tsv` – local paths. Relative paths are resolved against `bench/`.
- URL columns (optional) are used by `fetch_cami.sh` to download missing entries without disturbing existing files.

### 3.2 Derived CAMI subsets

`bench/data/` stores lightweight subsets generated from `cami_sample_0`. Recreate them with:

```bash
python ../tools/generate_cami_subsets.py \
  --fasta /data/cami/sample_0.fna \
  --mapping /data/cami/sample_0/2017.12.29_11.37.26_sample_0/contigs/gsa_mapping_new.tsv \
  --outdir $(pwd)/data
```

Each subset (`cami_i_lc`, `cami_ii_mousegut`, …) contains `contigs.fna`, `truth_contigs.tsv`, and `truth_profile.tsv`.

### 3.3 Shared reference FASTA

Database builders expect a unified FASTA such as `refsets/combined_subset.fasta`. You can generate one via `lib/subset_fasta.py` (see §5) or point `REF_FASTA` to a custom file. Keeping all builders aligned on the same FASTA is critical for fair comparisons.

---

## 4. Environment Setup

```bash
cd HYMET/bench
micromamba env create -f environment.yml -n hymet-benchmark
micromamba activate hymet-benchmark
```

The environment includes HYMET dependencies, Kraken2, Centrifuge, ganon, sourmash, MetaPhlAn4, taxonkit, minimap2, pandas, matplotlib, etc.

After any local edits, capture an exact lockfile:

```bash
micromamba env export -n hymet-benchmark > environment.lock.yml
```

---

## 5. Database Builders (`db/`)

All builders are idempotent: a `.build.stamp` prevents redundant work unless inputs change.

| Script | Description | Important env/flags |
|--------|-------------|---------------------|
| `build_kraken2.sh` | Builds Kraken2 DB against `REF_FASTA`. | `REF_FASTA`, `THREADS`, `KRAKEN2_DB_DIR`. |
| `build_centrifuge.sh` | Produces Centrifuge indices. | `REF_FASTA`, `THREADS`, `CENTRIFUGE_BMAX`, `CENTRIFUGE_DCV`. |
| `build_ganon2.sh` | Creates ganon2 HIBF using a generated manifest. | `REF_FASTA`, `THREADS`, `GANON_FILTER_SIZE`. |
| `build_sourmash.sh` | Generates sourmash sketch (`reference.k<k>.scaled<scaled>.sig`) + SBT. | `REF_FASTA`, `SOURMASH_KSIZE`, `SOURMASH_SCALED`. |

All builders require `HYMET/data/detailed_taxonomy.tsv`, produced by running `tools/build_id_map.py`.

To avoid memory pressure, consider trimming the shared FASTA:

```bash
python lib/subset_fasta.py \
  --input ../data/downloaded_genomes/combined_genomes.fasta \
  --output refsets/combined_subset.fasta \
  --max-seqs 1000 \
  --max-bases 500000000
export REF_FASTA=$(pwd)/refsets/combined_subset.fasta
```

---

## 6. Running the Benchmark

### 6.1 One-button driver

```bash
THREADS=16 METAPHLAN_THREADS=4 METAPHLAN_OPTS="--split_reads" ./run_all_cami.sh \
  --tools hymet,kraken2,centrifuge,ganon2,sourmash_gather,metaphlan4
```

Relevant flags:

| Flag | Meaning |
|------|---------|
| `--tools A,B,...` | Comma-separated list (or `all`). |
| `--no-build` | Skip database builders. |
| `--threads N` | Override thread count passed to runners. |
| `--max-samples N` | Process only the first N manifest rows. |
| `--resume` | Do not truncate `out/runtime_memory.tsv`. |

Environment variables respected by runners:

| Variable | Default | Description |
|----------|---------|-------------|
| `THREADS` | 8 | Baseline thread count (overridden by `--threads`). |
| `KEEP_HYMET_WORK` | 0 | When 0, HYMET deletes the ~52 GB minimap2 index after copying outputs. |
| `KEEP_GANON_RAW` / `KEEP_KRAKEN2_RAW` / `KEEP_CENTRIFUGE_RAW` | 0 | Preserve heavy raw files if set to 1. |
| `GANON_REL_CUTOFF` | 0 | ganon2 relative cutoff (0 disables). |
| `GANON_REL_FILTER` | 1 | ganon2 relative filter (1 keeps all matches). |
| `SOURMASH_TOP_HITS` | 500 | Limit gather results; set to `all` for no limit. |
| `METAPHLAN_DB_DIR` | `bench/db/metaphlan` | Directory containing MetaPhlAn bowtie2 indices. |
| `METAPHLAN_INDEX` | `mpa_vJun23_CHOCOPhlAnSGB_202307` | Index basename (passed as `-x`). |
| `METAPHLAN_THREADS` | `THREADS` | Bowtie2 threads (lower to reduce RAM). |
| `METAPHLAN_OPTS` | empty | Extra MetaPhlAn CLI options. |

Execution order:
1. Optionally run builders (skipped when `--no-build`).
2. Iterate manifest rows sequentially.
3. For each selected tool:
   - Run wrapper under `lib/measure.sh` (records runtime + memory).
   - Execute `lib/run_eval.sh` to compare predictions vs. truth.
4. After the manifest completes, call `aggregate_metrics.py` and optional figures.

### 6.2 Running individual tools

You can invoke any `run_<tool>.sh` directly, e.g.:

```bash
./run_ganon2.sh --sample cami_i_lc --contigs data/cami_i_lc/contigs.fna --threads 16
```

Outputs land in `out/<sample>/<tool>/`. Rebuild aggregate tables afterwards:

```bash
python aggregate_metrics.py --bench-root . --outdir out
python plot/make_figures.py  --bench-root . --outdir out
```

### 6.3 Partial reruns

- Rerun only ganon2 with relaxed cutoffs:
  ```bash
  THREADS=16 GANON_REL_CUTOFF=0 GANON_REL_FILTER=1 \
    ./run_all_cami.sh --tools ganon2 --no-build
  ```
- Rerun evaluation only:
  ```bash
  ./lib/run_eval.sh --sample cami_sample_0 --tool ganon2 \
    --pred-profile out/cami_sample_0/ganon2/profile.cami.tsv \
    --truth-profile /data/cami/sample_0/taxonomic_profile_0.txt \
    --pred-contigs out/cami_sample_0/ganon2/classified_sequences.tsv \
    --truth-contigs /data/cami/sample_0/.../gsa_mapping_new.tsv \
    --pred-fasta /data/cami/sample_0.fna --threads 16
  python aggregate_metrics.py --bench-root . --outdir out
  ```

---

## 7. Outputs

### 7.1 Per-sample directories

`out/<sample>/<tool>/` typically contains:

| Item | Description |
|------|-------------|
| `profile.cami.tsv` | CAMI-compatible profile (bioboxes format). |
| `classified_sequences.tsv` | Per-contig assignments (when provided by the tool). |
| `metadata.json` | Provenance (sample, tool, paths to key artefacts). |
| `resultados.paf` | HYMET-specific alignment output (symlinked to run directory). |
| `eval/profile_summary.tsv` | Rank-wise abundance metrics. |
| `eval/contigs_exact.tsv` and `eval/contigs_per_rank.tsv` | Contig metrics only emitted when usable pairs exist (`usable_pairs > 0`). |
| `eval/_debug_info.txt` | Diagnostic dump of the evaluation inputs. |

Tools that do **not** produce contig assignments (MetaPhlAn4, sourmash gather) will only have `profile_summary.tsv`; contig TSVs are deleted to avoid misleading zero rows.

### 7.2 Aggregates (`out/`)

| File | Description |
|------|-------------|
| `summary_per_tool_per_sample.tsv` | Profile metrics for every sample/tool/rank. |
| `leaderboard_by_rank.tsv` | Average metrics per tool/rank across samples. |
| `contig_accuracy_per_tool.tsv` | Contig accuracy; rows with `n <= 0` are filtered. |
| `runtime_memory.tsv` | Time, CPU, RSS, and I/O for each stage (produced by `measure.sh`). |
| `fig_*.png` | Visual summaries (F1, accuracy, L1/Bray–Curtis, per-sample stacks). |

---

## 8. Resource Management & Troubleshooting

### 8.1 Disk usage

- HYMET generates a ~52 GB minimap2 index in `out/<sample>/hymet/run/work/reference.mmi`. It is removed post-run unless `KEEP_HYMET_WORK=1`.
- MetaPhlAn database download leaves a 12 GB `.tar`; remove it after extraction to save space.
- Clearing `out/<sample>/<tool>/` before reruns avoids stale metrics and frees disk.

### 8.2 ganon2 coverage

- Relaxed defaults (`GANON_REL_CUTOFF=0`, `GANON_REL_FILTER=1`) retain more long-contig hits. Tighten them if you see false positives.
- Low coverage on CAMI subsets is expected because the shared reference FASTA is intentionally compact. Rebuild against a richer FASTA for higher recall.

### 8.3 MetaPhlAn4 memory

- Use `METAPHLAN_THREADS` to cap bowtie2 threads. `METAPHLAN_THREADS=4` kept memory < 20 GB during testing.
- Always set `METAPHLAN_DB_DIR` and ensure the index (`mpa_vJun23_CHOCOPhlAnSGB_202307.*`) exists; otherwise the runner will exit early with a detailed error.

### 8.4 Tools without contig output

- MetaPhlAn4 and sourmash gather do not produce contig-level predictions. The evaluation script now removes `contigs_*` TSVs when `usable_pairs == 0`, and the aggregator ignores them. This prevents meaningless zero rows in the contig accuracy table.

### 8.5 Logs & debugging

- Every tool invocation logs to STDOUT/STDERR via `measure.sh`—inspect direct output or `runtime_memory.tsv` for command lines and resource usage.
- `_debug_info.txt` in each `eval/` folder lists the exact files fed into `eval_cami.py`.

---

## 9. Suggested Workflow

1. **Activate environment**
   ```bash
   micromamba activate hymet-benchmark
   ```
2. **Validate manifest paths** (`cami_manifest.tsv`).
3. **Build databases**
   ```bash
   ./db/build_kraken2.sh
   ./db/build_centrifuge.sh
   ./db/build_ganon2.sh
   ./db/build_sourmash.sh
   # MetaPhlAn database must already exist in bench/db/metaphlan
   ```
4. **Run benchmark**
   ```bash
   THREADS=16 METAPHLAN_THREADS=4 METAPHLAN_OPTS="--split_reads" \
     ./run_all_cami.sh --threads 16 --tools all
   ```
5. **Regenerate aggregates/figures (optional)**
   ```bash
   python aggregate_metrics.py --bench-root . --outdir out
   python plot/make_figures.py  --bench-root . --outdir out
   ```
6. **Review results**
   - `out/runtime_memory.tsv`
   - `out/summary_per_tool_per_sample.tsv`
   - `out/contig_accuracy_per_tool.tsv`
   - `out/fig_*.png`

---

## 10. Extending the Harness

- **Add a new profiler**: create `run_<tool>.sh`, add a converter in `convert/`, register both in `run_all_cami.sh` (`TOOL_SCRIPTS` and optionally `TOOL_BUILDERS`), and document configuration knobs here.
- **Custom plots**: extend `plot/make_figures.py` or create additional scripts invoked after `aggregate_metrics.py`.
- **Alternate datasets**: edit `cami_manifest.tsv` to point at new samples or generate new subsets with `tools/generate_cami_subsets.py`.
- **Environment tweaks**: record additional packages in `environment.yml` and regenerate `environment.lock.yml` for reproducibility.

---

## 11. Checklist Before Sharing Results

- [ ] Environment created from `environment.yml` (or `environment.lock.yml`).
- [ ] MetaPhlAn database installed and referenced by `METAPHLAN_DB_DIR`.
- [ ] `cami_manifest.tsv` paths verified.
- [ ] Database builders completed successfully (check `.build.stamp` files).
- [ ] `run_all_cami.sh` run with appropriate `THREADS` and tool list.
- [ ] `aggregate_metrics.py` re-run after any manual tool execution.
- [ ] Disk usage monitored (`df -h`) to prevent quota issues.
- [ ] `out/` directory archived along with `runtime_memory.tsv`, aggregate TSVs, and plots.

Following the steps above yields a complete, reproducible CAMI benchmark for HYMET and companion profilers, ready for regression tracking or publication.
