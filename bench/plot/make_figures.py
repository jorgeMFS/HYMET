#!/usr/bin/env python3
"""Create benchmark plots from aggregated CAMI metrics."""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path

RANK_ORDER = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]


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


def load_runtime_rows(path: Path):
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
        try:
            plt.style.use("seaborn-v0_8-colorblind")
        except Exception:
            plt.style.use("ggplot")
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("matplotlib is required to generate plots. Install it or skip plotting.") from exc


def order_ranks(ranks):
    ranks = list(ranks)
    ordered = [rank for rank in RANK_ORDER if rank in ranks]
    extras = [rank for rank in ranks if rank not in RANK_ORDER]
    return ordered + sorted(extras)


def get_tool_colors(tools):
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("tab10")
    return {tool: cmap(i % cmap.N) for i, tool in enumerate(tools)}


def mean(values):
    values = [v for v in values if v is not None]
    return sum(values) / len(values) if values else 0.0


def clean_axis(ax):
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def plot_f1_by_rank(summary_rows, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    data = defaultdict(lambda: defaultdict(list))
    for row in summary_rows:
        tool = row["tool"]
        rank = row["rank"]
        data[rank][tool].append(safe_float(row.get("F1_%")))

    ranks = order_ranks(data.keys())
    tools = [tool for tool in tool_colors if any(tool in rank_data for rank_data in data.values())]

    fig, ax = plt.subplots(figsize=(11, 5.5))
    width = 0.8 / max(1, len(tools))
    x = list(range(len(ranks)))
    for idx, tool in enumerate(tools):
        offsets = [xi + idx * width for xi in x]
        means = [mean(data[rank].get(tool, [])) for rank in ranks]
        ax.bar(offsets, means, width=width, label=tool, color=tool_colors.get(tool))

    ax.set_xticks([xi + width * (len(tools) - 1) / 2 for xi in x])
    ax.set_xticklabels(ranks, rotation=20)
    ax.set_ylabel("F1 (%)")
    ax.set_title("Mean F1 by Rank")
    ax.legend(frameon=False, fontsize=9, ncol=min(len(tools), 4))
    ax.set_ylim(0, 100)
    clean_axis(ax)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_l1_bray(summary_rows, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    tool_metrics = defaultdict(lambda: {"L1": defaultdict(list), "Bray": defaultdict(list)})
    for row in summary_rows:
        tool = row["tool"]
        rank = row["rank"]
        tool_metrics[tool]["L1"][rank].append(safe_float(row.get("L1_total_variation_pctpts")))
        tool_metrics[tool]["Bray"][rank].append(safe_float(row.get("BrayCurtis_pct")))

    tools = [tool for tool in tool_colors if tool in tool_metrics]
    ranks = order_ranks({rank for tm in tool_metrics.values() for rank in tm["L1"].keys()})

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)
    for ax, metric_key, label in zip(axes, ["L1", "Bray"], ["L1 Total Variation (pct-pts)", "Bray-Curtis (%)"]):
        for tool in tools:
            means = [mean(tool_metrics[tool][metric_key].get(rank, [])) for rank in ranks]
            ax.plot(ranks, means, marker="o", linewidth=2, label=tool, color=tool_colors.get(tool))
        ax.set_title(label)
        ax.set_xticks(range(len(ranks)))
        ax.set_xticklabels(ranks, rotation=25)
        clean_axis(ax)
    axes[0].set_ylabel("Mean Value")
    axes[1].legend(frameon=False, fontsize=9, ncol=min(len(tools), 3))
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_accuracy(contig_rows, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    data = defaultdict(lambda: defaultdict(list))
    for row in contig_rows:
        if safe_float(row.get("n")) <= 0:
            continue
        tool = row["tool"]
        rank = row["rank"]
        data[rank][tool].append(safe_float(row.get("accuracy_percent")))

    if not data:
        return

    ranks = order_ranks(data.keys())
    tools = [tool for tool in tool_colors if any(tool in rank_data for rank_data in data.values())]

    fig, ax = plt.subplots(figsize=(11, 5.5))
    width = 0.8 / max(1, len(tools))
    x = list(range(len(ranks)))
    for idx, tool in enumerate(tools):
        offsets = [xi + idx * width for xi in x]
        means = [mean(data[rank].get(tool, [])) for rank in ranks]
        ax.bar(offsets, means, width=width, label=tool, color=tool_colors.get(tool))

    ax.set_xticks([xi + width * (len(tools) - 1) / 2 for xi in x])
    ax.set_xticklabels(ranks, rotation=20)
    ax.set_ylabel("Contig Accuracy (%)")
    ax.set_title("Mean Contig Accuracy by Rank")
    ax.legend(frameon=False, fontsize=9, ncol=min(len(tools), 4))
    ax.set_ylim(0, 100)
    clean_axis(ax)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_per_sample_stack(summary_rows, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    samples = sorted({row["sample"] for row in summary_rows})
    tools = [tool for tool in tool_colors if any(row["tool"] == tool for row in summary_rows)]

    data = {sample: defaultdict(list) for sample in samples}
    for row in summary_rows:
        data[row["sample"]][row["tool"]].append(safe_float(row.get("F1_%")))

    fig, ax = plt.subplots(figsize=(12, 5.5))
    bottoms = [0.0] * len(samples)
    x = list(range(len(samples)))
    for tool in tools:
        heights = []
        for sample in samples:
            values = data[sample].get(tool, [])
            heights.append(mean(values))
        ax.bar(x, heights, bottom=bottoms, label=tool, color=tool_colors.get(tool))
        bottoms = [b + h for b, h in zip(bottoms, heights)]

    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=25)
    ax.set_ylabel("Average F1 (%)")
    ax.set_title("Per-sample stacked F1 scores")
    ax.legend(frameon=False, fontsize=9, ncol=min(len(tools), 4))
    clean_axis(ax)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def summarise_runtime(runtime_rows):
    data = defaultdict(lambda: {"cpu_sec": [], "max_gb": []})
    for row in runtime_rows:
        if row.get("stage") != "run":
            continue
        tool = row.get("tool")
        if not tool:
            continue
        try:
            user = float(row.get("user_seconds", 0.0))
        except ValueError:
            user = 0.0
        try:
            sys_time = float(row.get("sys_seconds", 0.0))
        except ValueError:
            sys_time = 0.0
        try:
            rss = float(row.get("max_rss_gb", 0.0))
        except ValueError:
            rss = 0.0
        data[tool]["cpu_sec"].append(user + sys_time)
        data[tool]["max_gb"].append(rss)
    summary = {}
    for tool, metrics in data.items():
        cpu_vals = [v for v in metrics["cpu_sec"] if v >= 0]
        mem_vals = [v for v in metrics["max_gb"] if v >= 0]
        summary[tool] = {
            "cpu_min": mean(cpu_vals) / 60.0 if cpu_vals else 0.0,
            "max_gb": max(mem_vals) if mem_vals else 0.0,
        }
    return summary


def plot_runtime(summary, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    if not summary:
        return

    tools = sorted(summary.keys())
    values = [summary[t]["cpu_min"] for t in tools]
    x = list(range(len(tools)))

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(x, values, color=[tool_colors.get(t) for t in tools])
    ax.set_ylabel("CPU time (minutes)")
    ax.set_title("Mean CPU time per tool (run stage)")
    ax.set_xticks(x)
    ax.set_xticklabels(tools, rotation=25)
    clean_axis(ax)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{val:.1f}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_memory(summary, out_path: Path, tool_colors):
    import matplotlib.pyplot as plt

    if not summary:
        return

    tools = sorted(summary.keys())
    values = [summary[t]["max_gb"] for t in tools]
    x = list(range(len(tools)))

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(x, values, color=[tool_colors.get(t) for t in tools])
    ax.set_ylabel("Peak RSS (GB)")
    ax.set_title("Peak memory per tool (run stage)")
    ax.set_xticks(x)
    ax.set_xticklabels(tools, rotation=25)
    clean_axis(ax)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{val:.1f}", ha="center", va="bottom", fontsize=8)
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
    runtime_rows = load_runtime_rows(out_root / "runtime_memory.tsv")
    if not summary_rows:
        print("[plot] No summary data available; skipping figure generation.")
        return

    ensure_matplotlib()

    tools = sorted({row["tool"] for row in summary_rows})
    if contig_rows:
        tools = sorted(set(tools).union({row["tool"] for row in contig_rows}))
    tool_colors = get_tool_colors(tools)

    plot_f1_by_rank(summary_rows, out_root / "fig_f1_by_rank.png", tool_colors)
    plot_l1_bray(summary_rows, out_root / "fig_l1_braycurtis.png", tool_colors)
    if contig_rows:
        plot_accuracy(contig_rows, out_root / "fig_accuracy_by_rank.png", tool_colors)
    plot_per_sample_stack(summary_rows, out_root / "fig_per_sample_f1_stack.png", tool_colors)

    runtime_summary = summarise_runtime(runtime_rows)
    if runtime_summary:
        plot_runtime(runtime_summary, out_root / "fig_cpu_time_by_tool.png", tool_colors)
        plot_memory(runtime_summary, out_root / "fig_peak_memory_by_tool.png", tool_colors)


if __name__ == "__main__":
    main()
