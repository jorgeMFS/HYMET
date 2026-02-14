#!/usr/bin/env python3
"""Generate publication-quality figures for ZymoGut D6331 ground-truth validation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np


ZYMOGUT_ROOT = Path(__file__).resolve().parent
DEFAULT_RESULTS_DIR = ZYMOGUT_ROOT / "work" / "results"
DEFAULT_FIGURES_DIR = DEFAULT_RESULTS_DIR / "figures"

# ── Palette & style (matches bench/plot/make_figures.py) ──────────────────
COLOR_TRUTH = "#0f4c81"      # deep sapphire
COLOR_PREDICTED = "#ff6b6b"  # coral rose
COLOR_OVER = "#4d908e"       # muted teal  (overestimation)
COLOR_UNDER = "#f28482"      # salmon blush (underestimation)
COLOR_CONNECTOR = "#d0d5dd"  # light grey for connecting lines
COLOR_GENUS = "#6a4c93"      # muted purple for genus accent
BG_AXES = "#f4f6fb"
BG_FIGURE = "white"


def _setup_style():
    """Apply the shared HYMET plotting style."""
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except Exception:
        plt.style.use("ggplot")
    plt.rcParams.update({
        "axes.facecolor": BG_AXES,
        "figure.facecolor": BG_FIGURE,
        "axes.edgecolor": "#d0d5dd",
        "grid.alpha": 0.3,
        "grid.linestyle": "--",
        "axes.titleweight": "semibold",
        "axes.labelweight": "semibold",
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
    })


def _legend_style(legend):
    """White rounded legend box matching benchmark figures."""
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(0.65)


# ── Data loading ──────────────────────────────────────────────────────────

def load_comparison_table(path: Path) -> List[Dict]:
    """Load comparison table (species or genus)."""
    data = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Works for both species and genus tables
            name_key = "species" if "species" in row else "genus"
            data.append({
                "name": row[name_key],
                "truth_pct": float(row["truth_pct"]),
                "predicted_count": int(row["predicted_count"]),
                "predicted_pct": float(row["predicted_pct"]),
                "diff_pct": float(row["diff_pct"]),
            })
    return data


def load_metrics(path: Path) -> Dict[str, float]:
    """Load profile metrics."""
    metrics = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            metrics[row["metric"]] = float(row["value"])
    return metrics


def format_name(name: str) -> str:
    """Format taxon name for display (title case)."""
    name = name.replace("_", " ")
    if name and name[0].islower():
        name = name.title()
    return name


# ── Plot: Cleveland dot plot ─────────────────────────────────────────────

def plot_abundance_comparison(
    data: List[Dict],
    metrics: Dict[str, float],
    fig_path: Path,
    level: str = "species",
) -> None:
    """Paired dot plot comparing ground truth vs HYMET predicted abundances."""
    # Include taxa with truth > 0.1% OR predicted > 0.1%
    filtered = [d for d in data
                if d["truth_pct"] > 0.1 or d["predicted_pct"] > 0.1]
    filtered.sort(key=lambda x: x["truth_pct"])

    if not filtered:
        print(f"  Warning: no {level} with abundance > 0.1%")
        return

    names = [format_name(d["name"]) for d in filtered]
    truth = [d["truth_pct"] for d in filtered]
    pred = [d["predicted_pct"] for d in filtered]

    n = len(names)
    fig_h = max(5, 0.45 * n + 1.8)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    y = np.arange(n)

    for i in range(n):
        ax.plot([truth[i], pred[i]], [y[i], y[i]],
                color=COLOR_CONNECTOR, linewidth=1.8, zorder=2)

    ax.scatter(truth, y, s=90, color=COLOR_TRUTH, edgecolors="white",
               linewidths=0.8, zorder=3, label="Ground truth")
    ax.scatter(pred, y, s=90, color=COLOR_PREDICTED, edgecolors="white",
               linewidths=0.8, zorder=3, label="HYMET predicted")

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=10)
    ax.set_xlabel("Relative abundance (%)", fontsize=12)

    level_label = level.capitalize()
    ax.set_title(f"ZymoGut D6331: {level_label}-Level Abundance Comparison",
                 fontsize=14, weight="semibold", pad=12)

    pfx = "genus_" if level == "genus" else ""
    l1 = metrics.get(f"{pfx}l1_distance", 0)
    bc = metrics.get(f"{pfx}bray_curtis", 0)
    corr = metrics.get(f"{pfx}correlation", 0)
    ax.text(
        0.98, 0.03,
        f"L1 = {l1:.2f}   BC = {bc:.2f}   r = {corr:.3f}",
        transform=ax.transAxes, fontsize=9.5,
        ha="right", va="bottom",
        bbox=dict(facecolor="white", edgecolor="#d0d5dd",
                  boxstyle="round,pad=0.4", alpha=0.85),
    )

    leg = ax.legend(fontsize=10, loc="lower right",
                    bbox_to_anchor=(0.97, 0.08),
                    handlelength=1.4, columnspacing=0.8)
    _legend_style(leg)

    ax.set_xlim(left=-0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path.name}")


# ── Plot: Correlation scatter ────────────────────────────────────────────

def plot_correlation_scatter(
    data: List[Dict],
    metrics: Dict[str, float],
    fig_path: Path,
    level: str = "species",
) -> None:
    """Scatter plot of ground truth vs predicted abundances."""
    filtered = [d for d in data
                if d["truth_pct"] > 0.05 or d["predicted_pct"] > 0.05]
    if not filtered:
        print(f"  Warning: no {level} with abundance > 0.05%")
        return

    truth = np.array([d["truth_pct"] for d in filtered])
    pred = np.array([d["predicted_pct"] for d in filtered])
    names = [format_name(d["name"]) for d in filtered]

    pfx = "genus_" if level == "genus" else ""
    corr = metrics.get(f"{pfx}correlation", 0)

    fig, ax = plt.subplots(figsize=(7.5, 7.5))

    dot_color = COLOR_GENUS if level == "genus" else COLOR_TRUTH
    ax.scatter(truth, pred, s=85, color=dot_color, edgecolors="white",
               linewidths=0.8, alpha=0.85, zorder=3)

    axis_max = max(truth.max(), pred.max()) * 1.15
    ax.plot([0, axis_max], [0, axis_max], color="#adb5bd", ls="--", lw=1.4,
            zorder=2, label="y = x")

    # Label points that deviate significantly from the diagonal
    thresh = axis_max * 0.08
    for t, p, sp in zip(truth, pred, names):
        if abs(p - t) < thresh and not (t > axis_max * 0.15):
            continue
        if p > t:
            ha, xoff = "right", -8
        else:
            ha, xoff = "left", 8
        yoff = 5 if p >= t else -5
        ax.annotate(
            sp, (t, p), xytext=(xoff, yoff), textcoords="offset points",
            fontsize=9, color="#374151", ha=ha, va="center",
            bbox=dict(facecolor="white", edgecolor="none",
                      boxstyle="round,pad=0.15", alpha=0.75),
            arrowprops=dict(arrowstyle="-", color="#9ca3af", lw=0.8,
                            alpha=0.5),
        )

    level_label = level.capitalize()
    ax.set_xlabel("Ground truth (%)", fontsize=12)
    ax.set_ylabel("HYMET predicted (%)", fontsize=12)
    ax.set_title(f"ZymoGut D6331: {level_label}-Level Correlation (r = {corr:.3f})",
                 fontsize=13, weight="semibold", pad=10)
    ax.set_xlim(-0.5, axis_max)
    ax.set_ylim(-0.5, axis_max)
    ax.set_aspect("equal")

    leg = ax.legend(fontsize=10, loc="upper left", handlelength=1.4)
    _legend_style(leg)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path.name}")


# ── Plot: Diverging error bars ───────────────────────────────────────────

def plot_error_distribution(
    data: List[Dict],
    fig_path: Path,
    level: str = "species",
) -> None:
    """Horizontal diverging bar chart of prediction errors."""
    filtered = [d for d in data
                if d["truth_pct"] > 0.1 or d["predicted_pct"] > 0.1]
    filtered.sort(key=lambda x: x["diff_pct"])

    if not filtered:
        return

    names = [format_name(d["name"]) for d in filtered]
    diffs = [d["diff_pct"] for d in filtered]

    n = len(names)
    fig_h = max(4.5, 0.42 * n + 1.6)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    y = np.arange(n)
    colors = [COLOR_OVER if d >= 0 else COLOR_UNDER for d in diffs]

    bars = ax.barh(y, diffs, color=colors, edgecolor="white", linewidth=0.6,
                   height=0.65, zorder=3)

    for i, (bar, val) in enumerate(zip(bars, diffs)):
        if abs(val) < 0.3:
            continue
        sign = "+" if val > 0 else ""
        xpos = val
        ha = "left" if val >= 0 else "right"
        offset = 0.4 if val >= 0 else -0.4
        ax.text(
            xpos + offset, i, f"{sign}{val:.1f}",
            va="center", ha=ha, fontsize=8.5, color="#374151",
            bbox=dict(facecolor="white", edgecolor="none",
                      boxstyle="round,pad=0.12", alpha=0.8),
        )

    ax.axvline(x=0, color="#6b7280", linewidth=0.9, zorder=2)

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=10)
    ax.set_xlabel("Prediction error (predicted % \u2212 truth %)", fontsize=12)

    level_label = level.capitalize()
    ax.set_title(f"ZymoGut D6331: {level_label}-Level Prediction Error",
                 fontsize=14, weight="semibold", pad=12)

    leg = ax.legend(
        handles=[
            Patch(facecolor=COLOR_OVER, edgecolor="white", label="Overestimated"),
            Patch(facecolor=COLOR_UNDER, edgecolor="white", label="Underestimated"),
        ],
        fontsize=10, loc="lower right", handlelength=1.4,
    )
    _legend_style(leg)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path.name}")


# ── Plot: Combined 2-panel correlation (species vs genus) ────────────────

def plot_dual_correlation(
    species_data: List[Dict],
    genus_data: List[Dict],
    metrics: Dict[str, float],
    fig_path: Path,
) -> None:
    """Side-by-side scatter: species-level (left) vs genus-level (right)."""
    sp_f = [d for d in species_data if d["truth_pct"] > 0.05 or d["predicted_pct"] > 0.05]
    gn_f = [d for d in genus_data if d["truth_pct"] > 0.05 or d["predicted_pct"] > 0.05]

    if not sp_f or not gn_f:
        return

    fig, (ax_sp, ax_gn) = plt.subplots(1, 2, figsize=(14, 6.5))

    for ax, filt, pfx, color, label in [
        (ax_sp, sp_f, "", COLOR_TRUTH, "Species"),
        (ax_gn, gn_f, "genus_", COLOR_GENUS, "Genus"),
    ]:
        truth = np.array([d["truth_pct"] for d in filt])
        pred = np.array([d["predicted_pct"] for d in filt])
        names = [format_name(d["name"]) for d in filt]
        corr = metrics.get(f"{pfx}correlation", 0)

        ax.scatter(truth, pred, s=80, color=color, edgecolors="white",
                   linewidths=0.8, alpha=0.85, zorder=3)

        axis_max = max(truth.max(), pred.max()) * 1.15
        ax.plot([0, axis_max], [0, axis_max], color="#adb5bd", ls="--",
                lw=1.4, zorder=2)

        # Label major outliers
        thresh = axis_max * 0.10
        for t, p, sp in zip(truth, pred, names):
            if abs(p - t) < thresh and t < axis_max * 0.2:
                continue
            ha = "right" if p > t else "left"
            xoff = -8 if p > t else 8
            yoff = 5 if p >= t else -5
            ax.annotate(
                sp, (t, p), xytext=(xoff, yoff),
                textcoords="offset points", fontsize=8, color="#374151",
                ha=ha, va="center",
                bbox=dict(facecolor="white", edgecolor="none",
                          boxstyle="round,pad=0.12", alpha=0.7),
                arrowprops=dict(arrowstyle="-", color="#9ca3af", lw=0.7,
                                alpha=0.5),
            )

        ax.set_xlabel("Ground truth (%)", fontsize=11)
        ax.set_ylabel("HYMET predicted (%)", fontsize=11)
        ax.set_title(f"{label} level  (r = {corr:.3f})",
                     fontsize=12, weight="semibold", pad=8)
        ax.set_xlim(-0.5, axis_max)
        ax.set_ylim(-0.5, axis_max)
        ax.set_aspect("equal")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle("ZymoGut D6331: Abundance Correlation",
                 fontsize=14, weight="semibold", y=1.01)
    fig.tight_layout()
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path.name}")


# ── Plot: Metrics summary card ───────────────────────────────────────────

def plot_metrics_summary(metrics: Dict[str, float], fig_path: Path) -> None:
    """Compact table-style card showing key metrics at both levels."""
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.axis("off")

    rows = [
        ("L1 distance", f"{metrics.get('l1_distance', 0):.3f}",
         f"{metrics.get('genus_l1_distance', 0):.3f}"),
        ("Bray-Curtis", f"{metrics.get('bray_curtis', 0):.3f}",
         f"{metrics.get('genus_bray_curtis', 0):.3f}"),
        ("Correlation (r)", f"{metrics.get('correlation', 0):.3f}",
         f"{metrics.get('genus_correlation', 0):.3f}"),
        ("Genus recall", "\u2014",
         f"{metrics.get('genus_recall', 0) * 100:.1f}%"),
    ]

    col_labels = ["Metric", "Species", "Genus"]
    cell_text = [[r[0], r[1], r[2]] for r in rows]

    table = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.0, 1.8)

    # Style header
    for j in range(3):
        cell = table[0, j]
        cell.set_facecolor(COLOR_TRUTH)
        cell.set_text_props(color="white", weight="bold")
        cell.set_edgecolor("white")

    # Style data rows
    for i in range(1, len(rows) + 1):
        for j in range(3):
            cell = table[i, j]
            cell.set_facecolor("#f8f9fc" if i % 2 == 0 else "white")
            cell.set_edgecolor("#e5e7eb")
            if j == 0:
                cell.set_text_props(weight="semibold", ha="left")

    # Highlight the genus column with subtle accent
    for i in range(1, len(rows) + 1):
        cell = table[i, 2]
        cell.set_text_props(weight="bold", color=COLOR_GENUS)

    n_contigs = int(metrics.get("total_contigs", 0))
    ax.set_title(
        f"ZymoGut D6331: Classification Metrics ({n_contigs} contigs)",
        fontsize=13, weight="semibold", pad=16,
    )

    fig_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fig_path.name}")


# ── Metadata helper ───────────────────────────────────────────────────────

def save_metadata(
    species_data: List[Dict],
    genus_data: List[Dict],
    metrics: Dict[str, float],
    fig_dir: Path,
) -> None:
    payload = {
        "sample": "zymogut_d6331",
        "dataset": "ZymoBIOMICS Gut Microbiome Standard (D6331)",
        "metrics": metrics,
        "species_count": len([d for d in species_data if d["truth_pct"] > 0]),
        "genus_count": len([d for d in genus_data if d["truth_pct"] > 0]),
        "figures": [
            "fig_zymogut_genus_abundance.png",
            "fig_zymogut_genus_correlation.png",
            "fig_zymogut_dual_correlation.png",
            "fig_zymogut_metrics_summary.png",
            "fig_zymogut_abundance_comparison.png",
            "fig_zymogut_species_correlation.png",
            "fig_zymogut_species_error.png",
            "fig_zymogut_genus_error.png",
        ],
    }
    out = fig_dir / "zymogut_figures_metadata.json"
    out.write_text(json.dumps(payload, indent=2))
    print(f"  Saved: {out.name}")


# ── Main driver ───────────────────────────────────────────────────────────

def generate_figures(args: argparse.Namespace) -> None:
    _setup_style()

    results_dir = Path(args.results_dir).resolve()
    figures_dir = Path(args.figures_dir).resolve()

    species_path = results_dir / "comparison_table.tsv"
    genus_path = results_dir / "genus_comparison_table.tsv"
    metrics_path = results_dir / "profile_metrics.tsv"

    if not species_path.exists():
        raise SystemExit(f"Species comparison not found: {species_path}")
    if not metrics_path.exists():
        raise SystemExit(f"Metrics file not found: {metrics_path}")

    print("Loading data...")
    species_data = load_comparison_table(species_path)
    genus_data = load_comparison_table(genus_path) if genus_path.exists() else []
    metrics = load_metrics(metrics_path)

    print(f"  Species in comparison: {len(species_data)}")
    if genus_data:
        print(f"  Genera in comparison:  {len(genus_data)}")
    print(f"  Species r: {metrics.get('correlation', 0):.4f}")
    print(f"  Genus r:   {metrics.get('genus_correlation', 0):.4f}")

    print("\nGenerating figures...")
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Primary figures: genus-level (the strong result)
    if genus_data:
        plot_abundance_comparison(
            genus_data, metrics,
            figures_dir / "fig_zymogut_genus_abundance.png",
            level="genus",
        )
        plot_correlation_scatter(
            genus_data, metrics,
            figures_dir / "fig_zymogut_genus_correlation.png",
            level="genus",
        )
        plot_error_distribution(
            genus_data,
            figures_dir / "fig_zymogut_genus_error.png",
            level="genus",
        )

    # Combined dual-panel (genus vs species side-by-side)
    if genus_data:
        plot_dual_correlation(
            species_data, genus_data, metrics,
            figures_dir / "fig_zymogut_dual_correlation.png",
        )

    # Metrics summary card
    plot_metrics_summary(metrics, figures_dir / "fig_zymogut_metrics_summary.png")

    # Species-level figures (supplementary)
    plot_abundance_comparison(
        species_data, metrics,
        figures_dir / "fig_zymogut_abundance_comparison.png",
        level="species",
    )
    plot_correlation_scatter(
        species_data, metrics,
        figures_dir / "fig_zymogut_species_correlation.png",
        level="species",
    )
    plot_error_distribution(
        species_data,
        figures_dir / "fig_zymogut_species_error.png",
        level="species",
    )

    save_metadata(species_data, genus_data, metrics, figures_dir)
    print(f"\nGenerated figures in {figures_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        default=str(DEFAULT_RESULTS_DIR),
        help="Directory containing analysis results (default: work/results)",
    )
    parser.add_argument(
        "--figures-dir",
        help="Directory to store generated figures (default: <results-dir>/figures)",
    )
    args = parser.parse_args()
    if args.figures_dir is None:
        args.figures_dir = str(Path(args.results_dir).resolve() / "figures")
    return args


def main() -> None:
    args = parse_args()
    generate_figures(args)


if __name__ == "__main__":
    main()
