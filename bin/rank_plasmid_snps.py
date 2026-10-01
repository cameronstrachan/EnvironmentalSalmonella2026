#!/usr/bin/env python3
"""Rank SNPs in the core alignment by how specific they are to plasmid-carrying genomes.

Groups come from metadata/genome_plasmids.tsv: genomes with has_pESI_assembly
"yes" are the plasmid group, "no" the non-plasmid group. Alignment sequences
that are not in the metadata (the Reference) are left out of the statistics
and only used to report the reference allele.

For every variable column, the allele most enriched in the plasmid group is
chosen and scored:
  sensitivity  share of plasmid genomes carrying the allele
  specificity  share of non-plasmid genomes NOT carrying the allele
  youden_j     sensitivity + specificity - 1 (1.0 = allele carried by every
               plasmid genome and by no non-plasmid genome)
  fisher_p     two-sided Fisher's exact test on the 2x2 table
               (allele present/absent x plasmid/non-plasmid genome)
  fisher_q     Benjamini-Hochberg adjusted p-value across all variable sites
Sites are ranked by youden_j (highest first), then fisher_p, then position.

The core alignment has no gap columns, so aln_col is also the position on the
Reference sequence.

Requires: numpy, scipy

Usage:
  python3 bin/rank_plasmid_snps.py ALIGNMENT.aln METADATA.tsv OUTPUT.tsv

Example:
  python3 bin/rank_plasmid_snps.py \
      data/alignments/infantis_core_alignment_noambiguous.aln \
      metadata/genome_plasmids.tsv \
      results/ranked_snps/ranked_snps.tsv
"""

import csv
import sys

import numpy as np
from scipy.stats import fisher_exact

BASES = b"ACGT"
CHUNK = 200_000


def read_alignment(path):
    """Return (ids, uint8 matrix of sequences) in file order."""
    ids, rows = [], []
    current = None
    with open(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if current is not None:
                    rows.append(np.frombuffer(b"".join(current), dtype=np.uint8))
                ids.append(line[1:].strip().decode())
                current = []
            else:
                current.append(line.strip())
    rows.append(np.frombuffer(b"".join(current), dtype=np.uint8))
    return ids, np.vstack(rows)


def benjamini_hochberg(p_values):
    order = np.argsort(p_values)
    ranked = p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(p_values)
    q[order] = np.minimum(ranked, 1.0)
    return q


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    aln_path, meta_path, out_path = sys.argv[1:]

    meta = {row["genome_id"]: row["has_pESI_assembly"]
            for row in csv.DictReader(open(meta_path), delimiter="\t")}
    ids, matrix = read_alignment(aln_path)
    n_columns = matrix.shape[1]
    print(f"alignment: {len(ids)} sequences x {n_columns} columns", flush=True)

    with_plasmid = np.array([meta.get(i) == "yes" for i in ids])
    without_plasmid = np.array([meta.get(i) == "no" for i in ids])
    reference_rows = [k for k, i in enumerate(ids) if i not in meta]
    print(f"plasmid group: {with_plasmid.sum()}, non-plasmid group: "
          f"{without_plasmid.sum()}, not in metadata (excluded): "
          f"{[ids[k] for k in reference_rows]}", flush=True)
    reference = matrix[reference_rows[0]] if reference_rows else None

    n_with, n_without = int(with_plasmid.sum()), int(without_plasmid.sum())
    sub_with, sub_without = matrix[with_plasmid], matrix[without_plasmid]

    records = []  # one tuple of arrays per chunk
    for start in range(0, n_columns, CHUNK):
        end = min(start + CHUNK, n_columns)
        count_with = np.stack([(sub_with[:, start:end] == b).sum(axis=0) for b in BASES])
        count_without = np.stack([(sub_without[:, start:end] == b).sum(axis=0) for b in BASES])
        n_alleles = ((count_with + count_without) > 0).sum(axis=0)
        variable = np.flatnonzero(n_alleles > 1)
        if variable.size == 0:
            continue
        cw, cwo = count_with[:, variable], count_without[:, variable]
        delta = cw / n_with - cwo / n_without
        best = delta.argmax(axis=0)
        col = np.arange(variable.size)
        records.append((variable + start + 1, best, cw[best, col], cwo[best, col],
                        n_alleles[variable], cw.T, cwo.T))
        print(f"  columns {start + 1}-{end}: {variable.size} variable", flush=True)

    position = np.concatenate([r[0] for r in records])
    best = np.concatenate([r[1] for r in records])
    a = np.concatenate([r[2] for r in records])      # plasmid genomes with allele
    c = np.concatenate([r[3] for r in records])      # non-plasmid genomes with allele
    n_alleles = np.concatenate([r[4] for r in records])
    counts_with = np.vstack([r[5] for r in records])
    counts_without = np.vstack([r[6] for r in records])
    print(f"variable sites: {position.size}", flush=True)

    sensitivity = a / n_with
    specificity = 1 - c / n_without
    youden = sensitivity + specificity - 1

    # Fisher's exact test, computed once per distinct (a, c) pair.
    key = a.astype(np.int64) * (n_without + 1) + c
    unique_keys, inverse = np.unique(key, return_inverse=True)
    unique_p = np.empty(unique_keys.size)
    for i, k in enumerate(unique_keys):
        ai, ci = divmod(int(k), n_without + 1)
        unique_p[i] = fisher_exact([[ai, n_with - ai], [ci, n_without - ci]])[1]
    p_value = unique_p[inverse]
    q_value = benjamini_hochberg(p_value)

    order = np.lexsort((position, p_value, -youden))
    letters = [chr(b) for b in BASES]
    with open(out_path, "w") as out:
        out.write("\t".join([
            "rank", "aln_col", "ref_allele", "enriched_allele", "n_alleles",
            "plasmid_with_allele", "plasmid_total", "nonplasmid_with_allele",
            "nonplasmid_total", "sensitivity", "specificity", "youden_j",
            "fisher_p", "fisher_q", "allele_counts_plasmid", "allele_counts_nonplasmid",
        ]) + "\n")
        for rank, i in enumerate(order, start=1):
            ref = chr(reference[position[i] - 1]) if reference is not None else ""
            fmt = lambda counts: ",".join(f"{l}:{n}" for l, n in zip(letters, counts) if n)
            out.write("\t".join(map(str, [
                rank, position[i], ref, letters[best[i]], n_alleles[i],
                a[i], n_with, c[i], n_without,
                f"{sensitivity[i]:.4f}", f"{specificity[i]:.4f}", f"{youden[i]:.4f}",
                f"{p_value[i]:.3e}", f"{q_value[i]:.3e}",
                fmt(counts_with[i]), fmt(counts_without[i]),
            ])) + "\n")

    print(f"sites with youden_j == 1: {(youden >= 1 - 1e-12).sum()}")
    for cutoff in (0.99, 0.95, 0.9):
        print(f"sites with youden_j >= {cutoff}: {(youden >= cutoff).sum()}")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
