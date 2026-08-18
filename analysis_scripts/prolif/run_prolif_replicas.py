#!/usr/bin/env python3
"""Run the run1-compatible ProLIF workflow for Z944/mZ944 replicas 2 and 3."""

import argparse
import os
import shutil
from pathlib import Path

import MDAnalysis as mda
import matplotlib.pyplot as plt
import pandas as pd
import prolif as plf
from rdkit import Chem
from prolif.plotting.network import LigNetwork


HERE = Path(__file__).resolve().parent
CONFIGS = {
    "z944_run2": ("Z944", 2, Path("/work1/ted/June2023/WithLigand"), Path("/work1/ted/June2023/WithLigand/md_run2"), 6, 30, 0.020),
    "z944_run3": ("Z944", 3, Path("/work1/ted/June2023/WithLigand"), Path("/work1/ted/June2023/WithLigand/md_run3"), 6, 30, 0.020),
    "mz944_run1": ("mZ944", 1, Path("/mnt/passport/work1/ted/June2023/WithmodLigand"), Path("/work1/ted/June2023/WithmodLigand/md_run1"), 35, 64, 0.010),
    "mz944_run2": ("mZ944", 2, Path("/mnt/passport/work1/ted/June2023/WithmodLigand"), Path("/mnt/passport/work1/ted/June2023/WithmodLigand/md_run2"), 6, 30, 0.020),
    "mz944_run3": ("mZ944", 3, Path("/mnt/passport/work1/ted/June2023/WithmodLigand"), Path("/mnt/passport/work1/ted/June2023/WithmodLigand/md_run3"), 6, 30, 0.020),
}
FRAME_STRIDE = 300


def guess_element_symbols(atom_names):
    """Match the element-guessing behavior of the previous run1 analysis."""
    periodic_table = Chem.GetPeriodicTable()
    valid = {"H", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I", "Na", "Mg", "Ca", "K", "Zn", "Fe"}
    symbols = []
    for name in atom_names:
        raw = str(name).strip()
        if not raw:
            symbols.append("C")
            continue
        upper = raw.upper()
        special = {"CL": "Cl", "BR": "Br", "MG": "Mg", "CA": "Ca", "NA": "Na", "ZN": "Zn", "FE": "Fe"}
        symbol = next((value for prefix, value in special.items() if upper.startswith(prefix)), upper[0])
        if symbol not in valid:
            symbol = "C"
        try:
            periodic_table.GetAtomicNumber(symbol)
        except Exception:
            symbol = "C"
        symbols.append(symbol)
    return symbols


def run_one(key):
    system, replica, root, trajectory_dir, start_dcd, end_dcd, dt_ns = CONFIGS[key]
    out = HERE / "prolif" / system / f"run{replica}"
    out.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib"))

    dcds = [trajectory_dir / f"md-{i:02d}.dcd" for i in range(start_dcd, end_dcd + 1)]
    missing = [str(path) for path in dcds if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing DCD files for {key}: {missing}")

    print(f"[{key}] Loading {len(dcds)} DCD files...", flush=True)
    u = mda.Universe(str(root / "ionized.psf"), str(root / "ionized.pdb"), *map(str, dcds))
    u.add_TopologyAttr("elements", guess_element_symbols(u.atoms.names))
    protein = u.select_atoms("protein")
    ligand = u.select_atoms("resname DZR")
    if len(ligand.residues) != 1:
        raise RuntimeError(f"Expected one DZR residue, found {len(ligand.residues)}")

    sampled_indices = list(range(0, len(u.trajectory), FRAME_STRIDE))
    print(f"[{key}] {len(u.trajectory)} total frames; {len(sampled_indices)} sampled frames", flush=True)
    fp = plf.Fingerprint()
    fp.run(u.trajectory[::FRAME_STRIDE], ligand, protein, n_jobs=1, progress=True)
    df = fp.to_dataframe()
    df.to_csv(out / "interaction_fingerprint.csv")
    frequency = pd.DataFrame({
        "Interaction": [" | ".join(map(str, col)) for col in df.columns],
        "Occupancy": df.mean(axis=0).to_numpy(),
    }).sort_values("Occupancy", ascending=False)
    frequency.to_csv(out / "interaction_frequency.csv", index=False)

    pd.DataFrame({
        "sample_index": range(len(sampled_indices)),
        "trajectory_frame": sampled_indices,
        "time_ns": [i * dt_ns for i in sampled_indices],
    }).to_csv(out / "sampled_frames.csv", index=False)

    top = frequency.head(20).sort_values("Occupancy")
    fig, ax = plt.subplots(figsize=(11, 8))
    bars = ax.barh(top["Interaction"], 100 * top["Occupancy"], color="#2878B5")
    ax.bar_label(bars, fmt="%.1f%%", padding=3, fontsize=8)
    ax.set(xlabel="Occupancy (%)", ylabel="Ligand | Protein residue | Interaction type",
           title=f"Cav3.1–{system} run{replica} interaction occupancy")
    ax.set_xlim(0, 108)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "interaction_occupancy.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    lig = plf.Molecule.from_mda(ligand, force=True, inferrer=None)
    network = LigNetwork.from_fingerprint(
        fp, lig, kind="aggregate", threshold=0.20, use_coordinates=False, kekulize=True
    )
    network.save(out / "interaction_network.html", width="1200px", height="850px",
                 fontsize=16, show_interaction_data=True)
    js_source = Path("/work1/ted/June2023/WithLigand/md_run1/analyze/ProLif/output/vis-network.min.js")
    if js_source.is_file():
        shutil.copy2(js_source, out / js_source.name)
        html_path = out / "interaction_network.html"
        html = html_path.read_text()
        html = html.replace("https://unpkg.com/vis-network@9.0.4/dist/vis-network.min.js", "vis-network.min.js")
        html_path.write_text(html)

    with open(out / "summary.txt", "w") as handle:
        handle.write(f"Cav3.1-{system} run{replica} ProLIF Analysis\n")
        handle.write("=" * 44 + "\n\n")
        handle.write(f"DCD files      : md-{start_dcd:02d}.dcd through md-{end_dcd:02d}.dcd\n")
        handle.write(f"Total frames   : {len(u.trajectory)}\n")
        handle.write(f"Frames analyzed: {len(df)}\n")
        handle.write(f"Frame stride   : {FRAME_STRIDE}\n")
        handle.write(f"Nominal interval: 0-{(len(u.trajectory)-1)*dt_ns:.2f} ns\n\n")
        handle.write("Top interactions\n\n")
        for _, row in frequency.head(20).iterrows():
            handle.write(f"{row['Interaction']} {100 * row['Occupancy']:.1f}%\n")
    print(f"[{key}] Finished: {out}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cases", nargs="*", metavar="CASE")
    args = parser.parse_args()
    cases = args.cases or list(CONFIGS)
    unknown = sorted(set(cases) - set(CONFIGS))
    if unknown:
        parser.error(f"unknown case(s): {', '.join(unknown)}")
    for key in cases:
        run_one(key)


if __name__ == "__main__":
    main()
