#!/usr/bin/env python3
"""Match genome filenames to pESI assembly filenames without reading sequences."""

import argparse
import csv
from pathlib import Path


def find_assemblies(directory, plasmids=False):
    """Return {genome_id: filename} for FASTA files in one directory."""
    assemblies = {}
    for path in sorted(directory.iterdir()):
        if not path.is_file() or path.name.startswith("."):
            continue
        name = path
        if name.suffix.lower() in {".gz", ".bz2", ".xz"}:
            name = name.with_suffix("")
        if name.suffix.lower() not in {".fna", ".fasta", ".fa"}:
            continue

        genome_id = name.stem
        if plasmids:
            if not genome_id.endswith("_pESI"):
                raise ValueError(f"Expected <genome_id>_pESI filename: {path}")
            genome_id = genome_id[:-len("_pESI")]
        if not genome_id:
            raise ValueError(f"Missing genome ID in filename: {path}")
        if genome_id in assemblies:
            raise ValueError(f"Duplicate genome ID {genome_id} in {directory}")
        assemblies[genome_id] = path.name

    if not assemblies:
        raise ValueError(f"No FASTA assembly files found in {directory}")
    return assemblies


def main():
    project = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genomes-dir", type=Path,
                        default=project / "data/genome_assemblies")
    parser.add_argument("--plasmids-dir", type=Path,
                        default=project / "data/pESI_assemblies")
    parser.add_argument("--output", type=Path,
                        default=project / "metadata/genome_plasmids.tsv",
                        help="New TSV file to create; existing files are never overwritten")
    args = parser.parse_args()

    try:
        genomes = find_assemblies(args.genomes_dir)
        plasmids = find_assemblies(args.plasmids_dir, plasmids=True)
        unmatched = sorted(set(plasmids) - set(genomes))
        if unmatched:
            raise ValueError("Plasmids without matching genomes: " + ", ".join(unmatched))

        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["genome_id", "genome_file", "has_pESI_assembly", "pESI_file"])
            for genome_id in sorted(genomes):
                plasmid_file = plasmids.get(genome_id, "")
                writer.writerow([genome_id, genomes[genome_id],
                                 "yes" if plasmid_file else "no", plasmid_file])
    except (OSError, ValueError) as error:
        parser.error(str(error))

    print(f"Wrote {len(genomes)} genomes to {args.output}: "
          f"{len(plasmids)} with pESI assemblies, "
          f"{len(genomes) - len(plasmids)} without a matching pESI assembly.")


if __name__ == "__main__":
    main()
