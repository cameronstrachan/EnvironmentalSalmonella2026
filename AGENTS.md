# EnvironmentalSalmonella2026

This repository controls bioinformatic analyses for an environmental Salmonella project.

## Local vs remote

The local Mac repository contains code, workflows, configuration, documentation,
and small metadata tables. Large biological datasets and computational results
live on the remote server; the local data directories may contain only `.gitkeep`
placeholders.

Use the configured SSH alias to access the server:

```sh
ssh salmonella
cd ~/master/EnvironmentalSalmonella2026
```

The remote project is `/home/strachan/master/EnvironmentalSalmonella2026`.
The previously documented `~/EnvironmentalSalmonella2026` path does not exist.

For a single remote inspection command:

```sh
ssh salmonella 'cd ~/master/EnvironmentalSalmonella2026 && ls -la data'
```

Do not assume that local code changes are already present on the server. Remote
non-interactive shells may have a different `PATH`; check `command -v python3`
or `command -v git` when needed.

## Directory structure

```text
bin/                          Small analysis and metadata scripts
  create_plasmid_metadata.py   Map genome filenames to pESI assembly filenames
workflows/                    Workflow definitions and multi-step pipelines
metadata/                     Small metadata tables and their documentation
  README.md                   Column definitions, provenance, and usage
  genome_plasmids.tsv          One row per genome, with pESI assembly availability
data/                         Biological inputs, stored on the server
  genome_assemblies/           Genome FASTA files: <genome_id>.fna
  pESI_assemblies/             Plasmid FASTA files: <genome_id>_pESI.fasta
  alignments/                 Genome and plasmid alignments (.aln)
```

Use `results/` for generated analysis outputs, `logs/` for logs, and `tmp/` for
temporary files, creating them when needed. These directories and biological
data files are ignored by Git. Keep small metadata intended for version control
in `metadata/`, outside the ignored `data/` tree.

## Scripts and metadata

Prefer simple, readable Python 3 scripts in `bin/`. Use the standard library
when practical, clear variable names, and minimal dependencies. Document input
and output paths and how to run each script.

Inspect the actual filenames before mapping samples across datasets. Preserve
genome IDs as text and match complete IDs, rather than partial names. For the
current inputs, `10133415.fna` matches `10133415_pESI.fasta` by removing the exact
`_pESI` suffix from the plasmid filename stem.

`metadata/genome_plasmids.tsv` records whether a matching pESI assembly file
exists. A `no` does not establish the absence of pESI or other plasmids. See
`metadata/README.md` for the schema and provenance.

To create a new metadata snapshot on the server, from the project directory:

```sh
python3 bin/create_plasmid_metadata.py --output metadata/genome_plasmids_updated.tsv
```

The script reads filenames only and refuses to overwrite an existing output.
Choose a new output filename if the example above already exists.

## Safety

Do not delete, move, overwrite, or modify remote biological data unless explicitly asked.

Do not submit computationally expensive jobs unless explicitly asked.

Before destructive operations such as `rm`, `mv`, overwriting files, or deleting directories, ask for confirmation.

Prefer inspecting files and reporting findings before making changes.

## Git

Code changes should be made locally first unless there is a specific reason to edit remotely.

Do not commit or push changes unless explicitly asked.
