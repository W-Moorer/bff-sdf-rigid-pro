#!/usr/bin/env python3
"""Generate lightweight SVG figures from benchmark CSV files.

The script intentionally uses only the Python standard library so the paper
artifact path remains reproducible on a clean vcpkg/CMake setup.
"""

from __future__ import annotations

import csv
import math
import pathlib
from collections import defaultdict


ROOT = pathlib.Path(__file__).resolve().parents[1]
BENCH = ROOT / "results" / "benchmarks"
FIG = ROOT / "results" / "figures"


def read_csv(name: str) -> list[dict[str, str]]:
    path = BENCH / name
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fnum(row: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, default))
    except ValueError:
        return default


def esc(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_svg(path: pathlib.Path, width: int, height: int, body: str) -> None:
    path.write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n'
        '<rect width="100%" height="100%" fill="#ffffff"/>\n'
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#202124}.title{font-size:18px;font-weight:700}.label{font-size:11px}.small{font-size:10px}.axis{stroke:#4b5563;stroke-width:1}.grid{stroke:#e5e7eb;stroke-width:1}.bff{fill:#2563eb}.base{fill:#dc6429}.accent{fill:#0f766e}.muted{fill:#6b7280}</style>\n'
        + body
        + "\n</svg>\n",
        encoding="utf-8",
    )


def bar_chart(path: pathlib.Path, title: str, rows: list[tuple[str, float, str]], ylabel: str) -> None:
    width, height = 980, 520
    left, right, top, bottom = 90, 30, 56, 130
    plot_w = width - left - right
    plot_h = height - top - bottom
    vmax = max([v for _, v, _ in rows] + [1e-12])
    vmax *= 1.12
    parts = [f'<text class="title" x="{left}" y="30">{esc(title)}</text>']
    for i in range(6):
        y = top + plot_h - plot_h * i / 5
        val = vmax * i / 5
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}"/>')
        parts.append(f'<text class="small" x="{left-8}" y="{y+4:.1f}" text-anchor="end">{val:.3g}</text>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}"/>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top+plot_h}" x2="{width-right}" y2="{top+plot_h}"/>')
    n = max(1, len(rows))
    slot = plot_w / n
    bw = max(12, slot * 0.62)
    for i, (label, value, cls) in enumerate(rows):
        x = left + i * slot + (slot - bw) / 2
        h = plot_h * value / vmax
        y = top + plot_h - h
        parts.append(f'<rect class="{cls}" x="{x:.1f}" y="{y:.1f}" width="{bw:.1f}" height="{h:.1f}" rx="2"/>')
        parts.append(f'<text class="small" x="{x + bw/2:.1f}" y="{y-4:.1f}" text-anchor="middle">{value:.3g}</text>')
        parts.append(
            f'<text class="small" transform="translate({x + bw/2:.1f},{top+plot_h+12}) rotate(45)" text-anchor="start">{esc(label)}</text>'
        )
    parts.append(f'<text class="label" transform="translate(22,{top + plot_h/2:.1f}) rotate(-90)" text-anchor="middle">{esc(ylabel)}</text>')
    write_svg(path, width, height, "\n".join(parts))


def line_chart(path: pathlib.Path, title: str, series: dict[str, list[tuple[float, float]]], ylabel: str) -> None:
    width, height = 820, 480
    left, right, top, bottom = 86, 145, 56, 70
    plot_w = width - left - right
    plot_h = height - top - bottom
    xs = [x for vals in series.values() for x, _ in vals]
    ys = [y for vals in series.values() for _, y in vals]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if ymin <= 0:
        ymin = min(v for v in ys if v > 0) * 0.7
    ymax *= 1.25

    def sx(x: float) -> float:
        return left + (x - xmin) / max(1e-12, xmax - xmin) * plot_w

    def sy(y: float) -> float:
        ly = math.log10(max(y, 1e-16))
        lmin = math.log10(ymin)
        lmax = math.log10(ymax)
        return top + plot_h - (ly - lmin) / max(1e-12, lmax - lmin) * plot_h

    colors = ["#2563eb", "#dc6429", "#0f766e", "#7c3aed"]
    parts = [f'<text class="title" x="{left}" y="30">{esc(title)}</text>']
    for i in range(6):
        y = top + plot_h * i / 5
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{left+plot_w}" y2="{y:.1f}"/>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}"/>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}"/>')
    for idx, (name, vals) in enumerate(series.items()):
        vals = sorted(vals)
        color = colors[idx % len(colors)]
        points = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in vals)
        parts.append(f'<polyline fill="none" stroke="{color}" stroke-width="2.5" points="{points}"/>')
        for x, y in vals:
            parts.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3.5" fill="{color}"/>')
        parts.append(f'<rect x="{left+plot_w+22}" y="{top+idx*22}" width="12" height="12" fill="{color}"/>')
        parts.append(f'<text class="label" x="{left+plot_w+40}" y="{top+idx*22+11}">{esc(name)}</text>')
    parts.append(f'<text class="label" x="{left + plot_w/2:.1f}" y="{height-24}" text-anchor="middle">resolution / mesh resolution</text>')
    parts.append(f'<text class="label" transform="translate(24,{top+plot_h/2:.1f}) rotate(-90)" text-anchor="middle">{esc(ylabel)} log scale</text>')
    write_svg(path, width, height, "\n".join(parts))


def runtime_breakdown(path: pathlib.Path, rows: list[dict[str, str]]) -> None:
    width, height = 1040, 560
    left, right, top, bottom = 92, 130, 58, 145
    plot_w = width - left - right
    plot_h = height - top - bottom
    rows = sorted(rows, key=lambda r: fnum(r, "runtime_ms"), reverse=True)[:10]
    vmax = max([fnum(r, "runtime_ms") for r in rows] + [1e-12]) * 1.12
    colors = {
        "integration_ms": "#2563eb",
        "projection_ms": "#dc6429",
        "other_ms": "#6b7280",
    }
    labels = {
        "integration_ms": "integration",
        "projection_ms": "projection",
        "other_ms": "other",
    }
    parts = [f'<text class="title" x="{left}" y="30">Runtime breakdown by benchmark row</text>']
    for i in range(6):
        y = top + plot_h - plot_h * i / 5
        val = vmax * i / 5
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{left+plot_w}" y2="{y:.1f}"/>')
        parts.append(f'<text class="small" x="{left-8}" y="{y+4:.1f}" text-anchor="end">{val:.3g}</text>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}"/>')
    parts.append(f'<line class="axis" x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}"/>')
    n = max(1, len(rows))
    slot = plot_w / n
    bw = max(14, slot * 0.58)
    for i, row in enumerate(rows):
        x = left + i * slot + (slot - bw) / 2
        ybase = top + plot_h
        for key in ("integration_ms", "projection_ms", "other_ms"):
            val = fnum(row, key)
            h = plot_h * val / vmax
            ybase -= h
            parts.append(f'<rect x="{x:.1f}" y="{ybase:.1f}" width="{bw:.1f}" height="{h:.1f}" fill="{colors[key]}" rx="1"/>')
        label = row["scene"].replace("_", " ")
        parts.append(f'<text class="small" transform="translate({x + bw/2:.1f},{top+plot_h+12}) rotate(45)" text-anchor="start">{esc(label[:34])}</text>')
    lx = left + plot_w + 24
    for j, key in enumerate(("integration_ms", "projection_ms", "other_ms")):
        parts.append(f'<rect x="{lx}" y="{top + j*24}" width="13" height="13" fill="{colors[key]}"/>')
        parts.append(f'<text class="label" x="{lx+20}" y="{top + j*24 + 12}">{labels[key]}</text>')
    parts.append(f'<text class="label" transform="translate(24,{top+plot_h/2:.1f}) rotate(-90)" text-anchor="middle">runtime (ms)</text>')
    write_svg(path, width, height, "\n".join(parts))


def pipeline() -> None:
    width, height = 1050, 270
    labels = [
        "disk-like patches",
        "official BFF charts",
        "2D UV acceleration",
        "pulled-back SDF g(u)",
        "marching triangles",
        "chart projection",
        "CSV/OBJ/report",
    ]
    parts = ['<text class="title" x="40" y="34">BFF-SDF Contact Atlas pipeline</text>']
    x0, y, w, h, gap = 40, 92, 126, 68, 18
    for i, label in enumerate(labels):
        x = x0 + i * (w + gap)
        cls = "bff" if i in (1, 3, 5) else "accent"
        parts.append(f'<rect class="{cls}" x="{x}" y="{y}" width="{w}" height="{h}" rx="6" opacity="0.92"/>')
        words = label.split()
        parts.append(f'<text class="label" x="{x+w/2}" y="{y+28}" text-anchor="middle" fill="#fff">{esc(" ".join(words[:2]))}</text>')
        parts.append(f'<text class="label" x="{x+w/2}" y="{y+46}" text-anchor="middle" fill="#fff">{esc(" ".join(words[2:]))}</text>')
        if i + 1 < len(labels):
            ax = x + w + 3
            parts.append(f'<line x1="{ax}" y1="{y+h/2}" x2="{ax+gap-6}" y2="{y+h/2}" stroke="#374151" stroke-width="2"/>')
            parts.append(f'<polygon points="{ax+gap-6},{y+h/2-5} {ax+gap},{y+h/2} {ax+gap-6},{y+h/2+5}" fill="#374151"/>')
    parts.append('<text class="label" x="40" y="220">BFF is preprocessing-only; runtime updates rigid poses and evaluates the fixed chart field.</text>')
    write_svg(FIG / "pipeline.svg", width, height, "\n".join(parts))


def main() -> int:
    FIG.mkdir(parents=True, exist_ok=True)
    pipeline()

    sdf_rows = read_csv("sdf_resolution_study.csv")
    series: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in sdf_rows:
        series[row["method"]].append((fnum(row, "sdf_resolution"), fnum(row, "gap_error")))
    if series:
        line_chart(FIG / "grid_sdf_refinement.svg", "GridSDF error and projection refinement floor", dict(series), "gap error")

    mesh_rows = read_csv("mesh_resolution_study.csv")
    mesh_series = {
        "sphere-plane": [(fnum(r, "mesh_faces"), fnum(r, "area_error")) for r in mesh_rows if fnum(r, "area_error") > 0]
    }
    if mesh_series["sphere-plane"]:
        line_chart(FIG / "analytic_convergence.svg", "Analytic sphere-plane convergence", mesh_series, "relative area error")

    quality = [r for r in read_csv("chart_quality.csv") if r["scene"].startswith(("spherical_cap", "ellipsoid_cap"))]
    bars = []
    for row in quality:
        label = row["scene"].replace("spherical_", "sph_").replace("_curved", "") + "\n" + row["chart_method"].replace("fallback-", "")
        cls = "bff" if row["chart_method"] == "official-bff" else "base"
        bars.append((label, fnum(row, "mean_angle_distortion_deg"), cls))
    if bars:
        bar_chart(FIG / "bff_vs_planar_quality.svg", "Curved patch chart angle distortion", bars, "mean angle distortion (deg)")

    acc = read_csv("accuracy.csv")
    real = [r for r in acc if r.get("scene_class") == "real_mesh"]
    real_bars = [(r["scene"].replace("_real", ""), fnum(r, "contact_area"), "accent") for r in real]
    if real_bars:
        bar_chart(FIG / "real_asset_contact_visualization.svg", "Real-asset contact areas with reference-backed rows", real_bars, "contact area")

    if acc:
        runtime_breakdown(FIG / "runtime_breakdown.svg", acc)

    print(f"figures written to {FIG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
