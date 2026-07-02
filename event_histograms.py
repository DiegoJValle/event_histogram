#!/usr/bin/env python3
"""Create monthly DAS event histograms from annotation CSV files.

The program scans annotation CSV files, extracts the repeated event fields
(`event1_*`, `event2_*`, ...), groups events by month, and writes one stacked
histogram figure per month plus a text metrics log.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

TOTAL_LENGTH_KM = 50.0
LABELS = ("other", "noise", "vessel", "seismic")
LABEL_COLORS = {
    "other": "#7f7f7f",
    "noise": "#ff7f0e",
    "vessel": "#1f77b4",
    "seismic": "#2ca02c",
}
BAR_ALPHA = 0.55
CURVE_POINTS_PER_BIN = 12


@dataclass(frozen=True)
class Event:
    month: str
    km: float
    label: str
    source_file: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate monthly fiber-length histograms from DAS event annotation CSV files."
    )
    parser.add_argument("--input-dir", default="event_data", help="Folder containing annotation CSV files.")
    parser.add_argument("--output-dir", default="figures", help="Folder where figures are written.")
    parser.add_argument("--log-file", default="event_metrics.log", help="Path to the metrics log file.")
    parser.add_argument(
        "--bin-km",
        type=float,
        default=5.0,
        help="Histogram spacing/chunk size in kilometers along the 50 km fiber.",
    )
    parser.add_argument(
        "--total-length-km",
        type=float,
        default=TOTAL_LENGTH_KM,
        help="Total interrogated fiber length in kilometers.",
    )
    parser.add_argument("--dpi", type=int, default=150, help="Output figure DPI.")
    parser.add_argument(
        "--format",
        default="svg",
        choices=("png", "pdf", "svg"),
        help="Figure file format. SVG works without third-party packages; PNG/PDF require matplotlib.",
    )
    return parser.parse_args()


def normalize_label(raw_label: str | None) -> str:
    label = (raw_label or "").strip().lower()
    if label in {"noise", "vessel", "seismic"}:
        return label
    return "other"


def parse_datetime(value: str | None) -> datetime | None:
    if not value:
        return None
    value = value.strip()
    if not value:
        return None
    try:
        return datetime.fromisoformat(value)
    except ValueError:
        return None


def month_from_filename(filename: str) -> str | None:
    match = re.search(r"(20\d{2})(\d{2})\d{2}", filename)
    if not match:
        return None
    return f"{match.group(1)}-{match.group(2)}"


def iter_event_indexes(fieldnames: Iterable[str]) -> list[str]:
    indexes = set()
    for name in fieldnames:
        match = re.fullmatch(r"event(\d+)_label", name)
        if match:
            indexes.add(match.group(1))
    return sorted(indexes, key=int)


def read_events(input_dir: Path) -> list[Event]:
    csv_paths = sorted(input_dir.glob("*.csv"))
    if not csv_paths:
        raise FileNotFoundError(f"No CSV files found in {input_dir}")

    events: list[Event] = []
    for csv_path in csv_paths:
        with csv_path.open(newline="", encoding="utf-8-sig") as csv_file:
            reader = csv.DictReader(csv_file)
            indexes = iter_event_indexes(reader.fieldnames or [])
            for row in reader:
                fallback_dt = parse_datetime(row.get("init_dt"))
                fallback_month = (
                    fallback_dt.strftime("%Y-%m") if fallback_dt else month_from_filename(row.get("filename", ""))
                )
                for index in indexes:
                    label_value = row.get(f"event{index}_label")
                    km_value = row.get(f"event{index}_km")
                    if not (label_value and label_value.strip() and km_value and km_value.strip()):
                        continue
                    try:
                        km = float(km_value)
                    except ValueError:
                        continue
                    event_dt = parse_datetime(row.get(f"event{index}_dt"))
                    month = event_dt.strftime("%Y-%m") if event_dt else fallback_month
                    if not month:
                        continue
                    events.append(
                        Event(
                            month=month,
                            km=km,
                            label=normalize_label(label_value),
                            source_file=csv_path.name,
                        )
                    )
    return events


def bin_edges(total_length_km: float, bin_km: float) -> list[float]:
    if bin_km <= 0:
        raise ValueError("--bin-km must be greater than zero")
    if total_length_km <= 0:
        raise ValueError("--total-length-km must be greater than zero")
    bins = [round(i * bin_km, 10) for i in range(math.ceil(total_length_km / bin_km) + 1)]
    if bins[-1] < total_length_km:
        bins.append(total_length_km)
    bins[-1] = total_length_km
    return bins


def event_bin(km: float, edges: list[float]) -> int | None:
    if km < 0 or km > edges[-1]:
        return None
    if km == edges[-1]:
        return len(edges) - 2
    return max(0, min(len(edges) - 2, int(km / (edges[1] - edges[0]))))


def aggregate(events: list[Event], edges: list[float]) -> dict[str, list[Counter[str]]]:
    monthly = defaultdict(lambda: [Counter() for _ in range(len(edges) - 1)])
    for event in events:
        idx = event_bin(event.km, edges)
        if idx is not None:
            monthly[event.month][idx][event.label] += 1
    return dict(monthly)


def bin_centers(edges: list[float]) -> list[float]:
    return [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]


def smooth_distribution(counts: list[int], edges: list[float]) -> list[tuple[float, float]]:
    """Build a smooth count-scaled distribution curve for the histogram bars."""
    centers = bin_centers(edges)
    if not centers:
        return []

    bin_width = edges[1] - edges[0]
    sigma = max(bin_width * 0.75, 0.001)
    start = edges[0]
    stop = edges[-1]
    point_count = max(2, (len(edges) - 1) * CURVE_POINTS_PER_BIN + 1)
    step = (stop - start) / (point_count - 1)
    curve: list[tuple[float, float]] = []

    for point_index in range(point_count):
        x = start + point_index * step
        weighted_sum = 0.0
        weight_total = 0.0
        for center, count in zip(centers, counts):
            weight = math.exp(-0.5 * ((x - center) / sigma) ** 2)
            weighted_sum += count * weight
            weight_total += weight
        y = weighted_sum / weight_total if weight_total else 0.0
        curve.append((x, y))
    return curve


def create_figures(monthly_counts: dict[str, list[Counter[str]]], edges: list[float], output_dir: Path, dpi: int, fmt: str) -> list[Path]:
    if fmt == "svg":
        return create_svg_figures(monthly_counts, edges, output_dir)

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "matplotlib is required for PNG/PDF figures. Install it with "
            "'pip install matplotlib' or use '--format svg'."
        ) from exc

    output_dir.mkdir(parents=True, exist_ok=True)
    figure_paths: list[Path] = []
    width = (edges[1] - edges[0]) * 0.82
    centers = bin_centers(edges)

    for month in sorted(monthly_counts):
        fig, ax = plt.subplots(figsize=(12, 6))
        for label in LABELS:
            heights = [monthly_counts[month][i][label] for i in range(len(edges) - 1)]
            ax.bar(
                centers,
                heights,
                width=width,
                align="center",
                label=label,
                color=LABEL_COLORS[label],
                alpha=BAR_ALPHA,
                edgecolor="black",
                linewidth=0.7,
            )
            curve = smooth_distribution(heights, edges)
            if curve:
                ax.plot(
                    [point[0] for point in curve],
                    [point[1] for point in curve],
                    color=LABEL_COLORS[label],
                    linewidth=2.2,
                )
        ax.set_title(f"Overlapped DAS event histogram with distribution curves - {month}")
        ax.set_xlabel("Length along fiber (km)")
        ax.set_ylabel("Number of events recorded")
        ax.set_xlim(edges[0], edges[-1])
        ax.set_xticks(edges)
        ax.legend(title="Event class")
        ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()
        path = output_dir / f"event_histogram_{month}.{fmt}"
        fig.savefig(path, dpi=dpi)
        plt.close(fig)
        figure_paths.append(path)
    return figure_paths


def create_svg_figures(monthly_counts: dict[str, list[Counter[str]]], edges: list[float], output_dir: Path) -> list[Path]:
    """Write dependency-free overlapped histogram SVG files with distribution curves."""
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_paths: list[Path] = []
    width, height = 1200, 650
    margin_left, margin_right, margin_top, margin_bottom = 90, 220, 70, 90
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    centers = bin_centers(edges)
    bar_width_km = (edges[1] - edges[0]) * 0.82

    for month in sorted(monthly_counts):
        max_total = max((counter[label] for counter in monthly_counts[month] for label in LABELS), default=0) or 1
        x_scale = plot_width / (edges[-1] - edges[0])
        y_scale = plot_height / max_total
        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            "<style>text{font-family:Arial,Helvetica,sans-serif}.axis{stroke:#222;stroke-width:1.5}.grid{stroke:#ddd;stroke-width:1}.bar{stroke:#222;stroke-width:.8}.curve{fill:none;stroke-width:3;stroke-linejoin:round;stroke-linecap:round}</style>",
            f'<rect width="{width}" height="{height}" fill="white"/>',
            f'<text x="{width / 2}" y="35" text-anchor="middle" font-size="24">Overlapped DAS event histogram with distribution curves - {month}</text>',
        ]
        # Axes and grid.
        x0, y0 = margin_left, margin_top + plot_height
        parts.append(f'<line class="axis" x1="{x0}" y1="{y0}" x2="{x0 + plot_width}" y2="{y0}"/>')
        parts.append(f'<line class="axis" x1="{x0}" y1="{margin_top}" x2="{x0}" y2="{y0}"/>')
        for tick in range(0, max_total + 1, max(1, math.ceil(max_total / 5))):
            y = y0 - tick * y_scale
            parts.append(f'<line class="grid" x1="{x0}" y1="{y:.2f}" x2="{x0 + plot_width}" y2="{y:.2f}"/>')
            parts.append(f'<text x="{x0 - 10}" y="{y + 5:.2f}" text-anchor="end" font-size="14">{tick}</text>')
        for edge in edges:
            x = x0 + (edge - edges[0]) * x_scale
            parts.append(f'<line class="axis" x1="{x:.2f}" y1="{y0}" x2="{x:.2f}" y2="{y0 + 6}"/>')
            parts.append(f'<text x="{x:.2f}" y="{y0 + 25}" text-anchor="middle" font-size="12">{edge:g}</text>')

        # Overlapped bars, one transparent histogram per event class.
        bar_width = bar_width_km * x_scale
        for label in LABELS:
            for idx, counter in enumerate(monthly_counts[month]):
                value = counter[label]
                bar_height = value * y_scale
                if value:
                    x = x0 + (centers[idx] - edges[0]) * x_scale - bar_width / 2
                    y = y0 - bar_height
                    parts.append(
                        f'<rect class="bar" x="{x:.2f}" y="{y:.2f}" width="{bar_width:.2f}" '
                        f'height="{bar_height:.2f}" fill="{LABEL_COLORS[label]}" fill-opacity="{BAR_ALPHA}">'
                        f'<title>{label}: {value}</title></rect>'
                    )

            # Count-scaled smooth distribution curve over the matching bars.
            curve = smooth_distribution([counter[label] for counter in monthly_counts[month]], edges)
            if curve:
                path_points = []
                for point_index, (x_value, y_value) in enumerate(curve):
                    x = x0 + (x_value - edges[0]) * x_scale
                    y = y0 - y_value * y_scale
                    command = "M" if point_index == 0 else "L"
                    path_points.append(f"{command}{x:.2f},{y:.2f}")
                parts.append(
                    f'<path class="curve" d="{" ".join(path_points)}" stroke="{LABEL_COLORS[label]}">'
                    f'<title>{label} smoothed distribution</title></path>'
                )

        parts.append(f'<text x="{x0 + plot_width / 2}" y="{height - 25}" text-anchor="middle" font-size="16">Length along fiber (km)</text>')
        parts.append(f'<text x="25" y="{margin_top + plot_height / 2}" text-anchor="middle" font-size="16" transform="rotate(-90 25 {margin_top + plot_height / 2})">Number of events recorded</text>')
        legend_x, legend_y = width - margin_right + 45, margin_top + 20
        parts.append(f'<text x="{legend_x}" y="{legend_y - 15}" font-size="16" font-weight="bold">Event class</text>')
        for offset, label in enumerate(LABELS):
            y = legend_y + offset * 28
            parts.append(f'<rect x="{legend_x}" y="{y}" width="18" height="18" fill="{LABEL_COLORS[label]}" fill-opacity="{BAR_ALPHA}" stroke="#222"/>')
            parts.append(f'<text x="{legend_x + 28}" y="{y + 14}" font-size="15">{label}</text>')
        parts.append("</svg>")

        path = output_dir / f"event_histogram_{month}.svg"
        path.write_text("\n".join(parts), encoding="utf-8")
        figure_paths.append(path)
    return figure_paths


def write_log(events: list[Event], monthly_counts: dict[str, list[Counter[str]]], edges: list[float], log_file: Path) -> None:
    month_totals = {month: sum(sum(counter.values()) for counter in counters) for month, counters in monthly_counts.items()}
    month_labels = {month: Counter() for month in monthly_counts}
    chunk_totals: Counter[str] = Counter()

    for month, counters in monthly_counts.items():
        for idx, counter in enumerate(counters):
            label = f"{edges[idx]:g}-{edges[idx + 1]:g} km"
            chunk_total = sum(counter.values())
            chunk_totals[label] += chunk_total
            month_labels[month].update(counter)

    log_file.parent.mkdir(parents=True, exist_ok=True) if log_file.parent != Path("") else None
    with log_file.open("w", encoding="utf-8") as handle:
        handle.write("DAS Event Histogram Metrics\n")
        handle.write("===========================\n\n")
        handle.write(f"Total events counted: {len(events)}\n")
        handle.write(f"Months processed: {len(monthly_counts)}\n")
        handle.write(f"Fiber range: 0-{edges[-1]:g} km\n")
        handle.write(f"Chunk spacing: {edges[1] - edges[0]:g} km\n\n")

        handle.write("Months ordered by vessel events\n")
        handle.write("-------------------------------\n")
        for month, labels in sorted(month_labels.items(), key=lambda item: item[1]["vessel"], reverse=True):
            handle.write(f"{month}: vessels={labels['vessel']}, total={month_totals[month]}, noise={labels['noise']}, seismic={labels['seismic']}, other={labels['other']}\n")

        handle.write("\nFiber chunks with most events recorded\n")
        handle.write("--------------------------------------\n")
        for chunk, count in chunk_totals.most_common():
            handle.write(f"{chunk}: {count}\n")

        handle.write("\nMonthly label distribution\n")
        handle.write("--------------------------\n")
        for month in sorted(month_labels):
            labels = month_labels[month]
            handle.write(f"{month}: " + ", ".join(f"{label}={labels[label]}" for label in LABELS) + "\n")

        handle.write("\nTop chunks per month\n")
        handle.write("--------------------\n")
        for month in sorted(monthly_counts):
            ranked = []
            for idx, counter in enumerate(monthly_counts[month]):
                ranked.append((sum(counter.values()), f"{edges[idx]:g}-{edges[idx + 1]:g} km", counter))
            handle.write(f"{month}:\n")
            for total, chunk, counter in sorted(ranked, reverse=True)[:5]:
                handle.write(f"  {chunk}: total={total}, " + ", ".join(f"{label}={counter[label]}" for label in LABELS) + "\n")


def main() -> int:
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    log_file = Path(args.log_file)
    edges = bin_edges(args.total_length_km, args.bin_km)
    events = read_events(input_dir)
    monthly_counts = aggregate(events, edges)
    if not monthly_counts:
        print("No events found in the input CSV files.", file=sys.stderr)
        return 1
    figures = create_figures(monthly_counts, edges, output_dir, args.dpi, args.format)
    write_log(events, monthly_counts, edges, log_file)
    print(f"Created {len(figures)} figure(s) in {output_dir}")
    print(f"Wrote metrics log to {log_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
