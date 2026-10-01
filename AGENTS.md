# EnvironmentalSalmonella2026

This repository controls bioinformatic analyses for an environmental Salmonella project.

## Local vs remote

The local Mac repository contains code, workflows, documentation, and small
metadata tables. Biological datasets and computational results live on the server.

Use the configured SSH alias to access the server:

```sh
ssh salmonella
cd ~/master/EnvironmentalSalmonella2026
```

Use `~/miniconda3/bin/git` for Git commands on the server.

## Conda environment

The project environment is `environmental-salmonella-2026`. Its dependencies
are recorded in `environment.yml`. Activate it in each server session:

```sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate environmental-salmonella-2026
```

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

Keep small metadata intended for version control in `metadata/`. Biological
data under `data/` is ignored by Git.

## Scripts and metadata

Prefer simple, readable Python 3 scripts in `bin/`. Use the standard library
when practical, clear variable names, and minimal dependencies. Document input
and output paths and how to run each script.

Inspect filenames before mapping samples across datasets. Preserve genome IDs
as text and match complete IDs: `<genome_id>.fna` pairs with
`<genome_id>_pESI.fasta`.

`metadata/genome_plasmids.tsv` records whether a matching pESI assembly file
exists. A `no` does not establish the absence of pESI or other plasmids. See
`metadata/README.md` for the schema and provenance.

## Safety

Do not delete, move, overwrite, or modify remote biological data unless explicitly asked.

Do not submit computationally expensive jobs unless explicitly asked.

Before destructive operations such as `rm`, `mv`, overwriting files, or deleting directories, ask for confirmation.

Prefer inspecting files and reporting findings before making changes.

## Git

Code changes should be made locally first unless there is a specific reason to edit remotely.

Do not commit or push changes unless explicitly asked.
