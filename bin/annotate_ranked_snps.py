#!/usr/bin/env python3
"""Add Bakta gene annotations to the ranked SNP table.

Each SNP's aln_col is a position on the Reference sequence (the core
alignment has no gap columns), so it is matched to the features in the
Bakta GFF3 of that sequence.

Columns added:
  feature_type, locus_tag, gene, product, strand   the feature containing the SNP
  other_allele, codon_pos, other_aa, enriched_aa, effect   for SNPs inside a CDS:
      other_allele is the most common other allele at the site, codon_pos the
      position in the codon (1-3), and other_aa / enriched_aa the amino acids
      encoded when the site carries the other or the plasmid-enriched allele.
      effect is the change from other to enriched allele: synonymous, missense,
      or nonsense (a stop gained or lost).
  nearest_upstream, nearest_downstream             locus tags of the closest
      features on either side, for SNPs outside any feature (intergenic)

Requires: Python 3 standard library only

Usage:
  python3 bin/annotate_ranked_snps.py RANKED.tsv REFERENCE.gff3 REFERENCE.fasta OUTPUT.tsv

Example:
  python3 bin/annotate_ranked_snps.py \
      results/ranked_snps/ranked_snps.tsv \
      results/reference_annotation/Reference.gff3 \
      data/reference/Reference.fasta \
      results/ranked_snps/ranked_snps.annotated.tsv
"""

import bisect
import csv
import sys
import urllib.parse

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: AMINO[16 * i + 4 * j + k]
               for i, a in enumerate(BASES)
               for j, b in enumerate(BASES)
               for k, c in enumerate(BASES)}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def read_fasta(path):
    return "".join(line.strip() for line in open(path) if not line.startswith(">")).upper()


def read_features(path):
    """Return features sorted by start: (start, end, type, locus_tag, gene, product, strand)."""
    features = []
    for line in open(path):
        if line.startswith("#") or line.startswith(">") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] in ("region", "gene"):
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        attrs = {k: urllib.parse.unquote(v) for k, v in attrs.items()}
        features.append((int(f[3]), int(f[4]), f[2], attrs.get("locus_tag", ""),
                         attrs.get("Name", attrs.get("gene", "")) if f[2] == "CDS" else attrs.get("Name", ""),
                         attrs.get("product", ""), f[6]))
    features.sort()
    return features


def codon_effect(sequence, feature, position, other, allele):
    """Codon position and amino acids for the other and enriched alleles."""
    start, end, strand = feature[0], feature[1], feature[6]
    if strand == "+":
        offset = position - start
        codon_start = start + offset - offset % 3
        codon = sequence[codon_start - 1:codon_start + 2]
        index = position - codon_start
        first = codon[:index] + other + codon[index + 1:]
        second = codon[:index] + allele + codon[index + 1:]
    else:
        offset = end - position
        codon_start = end - (offset - offset % 3)
        codon = sequence[codon_start - 3:codon_start][::-1].translate(COMPLEMENT)
        index = offset % 3
        first = codon[:index] + other.translate(COMPLEMENT) + codon[index + 1:]
        second = codon[:index] + allele.translate(COMPLEMENT) + codon[index + 1:]
    if len(codon) != 3 or "N" in codon:
        return "", "", "", ""
    other_aa, new_aa = CODON_TABLE[first], CODON_TABLE[second]
    if other_aa == new_aa:
        effect = "synonymous"
    elif new_aa == "*" or other_aa == "*":
        effect = "nonsense"
    else:
        effect = "missense"
    return index + 1, other_aa, new_aa, effect


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    ranked_path, gff_path, fasta_path, out_path = sys.argv[1:]

    sequence = read_fasta(fasta_path)
    features = read_features(gff_path)
    starts = [f[0] for f in features]
    print(f"{len(features)} features, sequence length {len(sequence)}")

    rows = list(csv.DictReader(open(ranked_path), delimiter="\t"))
    new_columns = ["feature_type", "locus_tag", "gene", "product", "strand",
                   "other_allele", "codon_pos", "other_aa", "enriched_aa", "effect",
                   "nearest_upstream", "nearest_downstream"]
    counts = {}
    with open(out_path, "w") as out:
        out.write("\t".join(list(rows[0].keys()) + new_columns) + "\n")
        for row in rows:
            position = int(row["aln_col"])
            allele = row["enriched_allele"]
            totals = {}
            for field in ("allele_counts_plasmid", "allele_counts_nonplasmid"):
                for item in row[field].split(","):
                    base, number = item.split(":")
                    totals[base] = totals.get(base, 0) + int(number)
            others = [b for b in sorted(totals, key=totals.get, reverse=True) if b != allele]
            other = others[0] if others else ""
            # Features that start at or before the SNP; the SNP is inside if it ends after.
            idx = bisect.bisect_right(starts, position)
            hits = [f for f in features[max(0, idx - 200):idx] if f[1] >= position]
            cds = [f for f in hits if f[2] == "CDS"]
            chosen = cds[-1] if cds else (hits[-1] if hits else None)
            extra = {k: "" for k in new_columns}
            extra["other_allele"] = other
            if chosen:
                extra.update(feature_type=chosen[2], locus_tag=chosen[3], gene=chosen[4],
                             product=chosen[5], strand=chosen[6])
                if chosen[2] == "CDS" and other:
                    pos, other_aa, new_aa, effect = codon_effect(sequence, chosen, position, other, allele)
                    extra.update(other_allele=other, codon_pos=pos, other_aa=other_aa,
                                 enriched_aa=new_aa, effect=effect)
                counts[chosen[2] + ("/" + extra["effect"] if extra["effect"] else "")] = \
                    counts.get(chosen[2] + ("/" + extra["effect"] if extra["effect"] else ""), 0) + 1
            else:
                extra["feature_type"] = "intergenic"
                before = [f for f in features[:idx] if f[1] < position]
                after = features[idx:idx + 1]
                extra["nearest_upstream"] = before[-1][3] if before else ""
                extra["nearest_downstream"] = after[0][3] if after else ""
                counts["intergenic"] = counts.get("intergenic", 0) + 1
            out.write("\t".join([row[k] for k in rows[0].keys()]
                                + [str(extra[k]) for k in new_columns]) + "\n")
    for key in sorted(counts):
        print(f"  {key}: {counts[key]}")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
