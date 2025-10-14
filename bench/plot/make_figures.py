#!/usr/bin/env python3
"""Create benchmark plots from aggregated CAMI metrics."""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path


def load_table(path: Path, key_cols):
    if not path.is_file():
        return {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        data = {}
        for row in reader:
            key = tuple(row[col] for col in key_cols)
            data[key] = row
    return data


def load_rows(path: Path):
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def ensure_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: F401
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("matplotlib is required to generate plots. Install it or skip plotting.") from exc


def plot_f1_by_rank(summary_rows, out_path: Path):
    import matplotlib.pyplot as plt

    data = defaultdict(lambda: defaultdict(list))
    for row in summary_rows:
        tool = row["tool"]
        rank = row["rank"]
        data[rank][tool].append(safe_float(row.get("F1_%")))

    ranks = sorted(data.keys())
    tools = sorted({tool for rank_data in data.values() for tool in rank_data})

    fig, ax = plt.subplots(figsize=(10, 5))
    width = 0.8 / max(1, len(tools))
    x = range(len(ranks))
    for idx, tool in enumerate(tools):
        offsets = [xi + idx * width for xi in x]
        means = [sum(data[rank].get(tool, [])) / max(1, len(data[rank].get(tool, []))) for rank in ranks]
        ax.bar(offsets, means, width=width, label=tool)

    ax.set_xticks([xi + width * (len(tools) - 1) / 2 for xi in x])
    ax.set_xticklabels(ranks, rotation=20)
    ax.set_ylabel("F1 (%)")
    ax.set_title("Mean F1 by Rank")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_l1_bray(summary_rows, out_path: Path):
    import matplotlib.pyplot as plt

    tool_metrics = defaultdict(lambda: {"L1": defaultdict(list), "Bray": defaultdict(list)})
    for row in summary_rows:
        tool = row["tool"]
        rank = row["rank"]
        tool_metrics[tool]["L1"][rank].append(safe_float(row.get("L1_total_variation_pctpts")))
        tool_metrics[tool]["Bray"][rank].append(safe_float(row.get("BrayCurtis_pct")))

    tools = sorted(tool_metrics.keys())
    ranks = sorted({rank for tm in tool_metrics.values() for rank in tm["L1"].keys()})

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, metric_key, label in zip(axes, ["L1", "Bray"], ["L1 Total Variation (pct-pts)", "Bray-Curtis (%)"]):
        for tool in tools:
            means = [sum(tool_metrics[tool][metric_key].get(rank, [])) / max(1, len(tool_metrics[tool][metric_key].get(rank, []))) for rank in ranks]
            ax.plot(ranks, means, marker="o", label=tool)
        ax.set_title(label)
        ax.set_xticks(range(len(ranks)))
        ax.set_xticklabels(ranks, rotation=25)
    axes[0].set_ylabel("Mean Value")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_accuracy(contig_rows, out_path: Path):
    import matplotlib.pyplot as plt

    data = defaultdict(lambda: defaultdict(list))
    for row in contig_rows:
        if safe_float(row.get("n")) <= 0:
            continue
        tool = row["tool"]
        rank = row["rank"]
        data[rank][tool].append(safe_float(row.get("accuracy_percent")))

    ranks = sorted(data.keys())
    tools = sorted({tool for rank_data in data.values() for tool in rank_data})

    fig, ax = plt.subplots(figsize=(10, 5))
    width = 0.8 / max(1, len(tools))
    x = range(len(ranks))
    for idx, tool in enumerate(tools):
        offsets = [xi + idx * width for xi in x]
        means = [sum(data[rank].get(tool, [])) / max(1, len(data[rank].get(tool, []))) for rank in ranks]
        ax.bar(offsets, means, width=width, label=tool)

    ax.set_xticks([xi + width * (len(tools) - 1) / 2 for xi in x])
    ax.set_xticklabels(ranks, rotation=20)
    ax.set_ylabel("Contig Accuracy (%)")
    ax.set_title("Mean Contig Accuracy by Rank")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_per_sample_stack(summary_rows, out_path: Path):
    import matplotlib.pyplot as plt

    samples = sorted({row["sample"] for row in summary_rows})
    tools = sorted({row["tool"] for row in summary_rows})

    data = {sample: defaultdict(list) for sample in samples}
    for row in summary_rows:
        data[row["sample"]][row["tool"]].append(safe_float(row.get("F1_%")))

    fig, ax = plt.subplots(figsize=(12, 5))
    bottoms = [0.0] * len(samples)
    x = range(len(samples))
    for tool in tools:
        heights = []
        for sample in samples:
            values = data[sample].get(tool, [])
            heights.append(sum(values) / max(1, len(values)))
        ax.bar(x, heights, bottom=bottoms, label=tool)
        bottoms = [b + h for b, h in zip(bottoms, heights)]

    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=25)
    ax.set_ylabel("Average F1 (%)")
    ax.set_title("Per-sample stacked F1 scores")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate plots from CAMI aggregate metrics.")
    ap.add_argument("--bench-root", default=str(Path(__file__).resolve().parent.parent), help="Bench directory root.")
    ap.add_argument("--outdir", default="out", help="Relative output directory for figures.")
    args = ap.parse_args()

    bench_root = Path(args.bench_root)
    out_root = bench_root / args.outdir
    out_root.mkdir(parents=True, exist_ok=True)

    summary_path = out_root / "summary_per_tool_per_sample.tsv"
    contig_path = out_root / "contig_accuracy_per_tool.tsv"

    summary_rows = load_rows(summary_path)
    contig_rows = load_rows(contig_path)
    if not summary_rows:
        print("[plot] No summary data available; skipping figure generation.")
        return

    ensure_matplotlib()

    plot_f1_by_rank(summary_rows, out_root / "fig_f1_by_rank.png")
    plot_l1_bray(summary_rows, out_root / "fig_l1_braycurtis.png")
    if contig_rows:
        plot_accuracy(contig_rows, out_root / "fig_accuracy_by_rank.png")
    plot_per_sample_stack(summary_rows, out_root / "fig_per_sample_f1_stack.png")


if __name__ == "__main__":
    main()
