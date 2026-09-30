# Genome plasmid metadata

`genome_plasmids.tsv` has one row per genome assembly. It matches the exact
genome ID in `data/genome_assemblies/<genome_id>.fna` to
`data/pESI_assemblies/<genome_id>_pESI.fasta`. For example, `10133415.fna`
matches `10133415_pESI.fasta`.

Columns:

- `genome_id`: filename ID, kept as text.
- `genome_file`: genome assembly filename.
- `has_pESI_assembly`: `yes` if a matching pESI assembly file exists; otherwise `no`.
- `pESI_file`: matching plasmid filename, blank when there is no match.

This records the availability of pESI assembly files. A `no` does not establish
that a genome lacks pESI or other plasmids. The script only reads filenames,
not sequence contents.

The initial table was generated on 2026-09-30 from the server project at
`salmonella:~/master/EnvironmentalSalmonella2026`: 750 genomes, 479 matching
pESI assemblies, 271 genomes without a matching file, and no unmatched plasmids.

Run from the project directory on the server, using Python 3 and no extra packages:

```sh
python3 bin/create_plasmid_metadata.py
```

The script creates `metadata/` when needed and refuses to overwrite an existing
output. To generate a new snapshot once the default output exists:

```sh
python3 bin/create_plasmid_metadata.py --output metadata/genome_plasmids_updated.tsv
```

Use `--genomes-dir` and `--plasmids-dir` to select different input directories.
The script accepts `.fna`, `.fasta`, and `.fa`, optionally compressed with
`.gz`, `.bz2`, or `.xz`. It rejects duplicate IDs, empty input collections,
unexpected plasmid FASTA names, and plasmids without matching genome files.
