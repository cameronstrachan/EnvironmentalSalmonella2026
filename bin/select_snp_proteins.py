#!/usr/bin/env python3
"""Select proteins affected by, or next to, the filtered plasmid-lineage SNPs.

Takes the filtered ranked SNP table and the Bakta annotation of the Reference.
Selected genes (CDS):
  missense / nonsense   the SNP lies in the CDS and changes (or gains/loses a
                        stop in) the encoded protein
  adjacent              for each intergenic SNP, the nearest CDS on each side
Synonymous SNPs were already removed from the filtered table.

"Nonsense" rows: effect compares the plasmid-enriched allele with the most
common other allele. If the Reference codon is a stop codon at the end of the
CDS, the lineage protein ends there and the other allele reads through, so
the other-allele protein is longer. Both cases are recorded: the Reference
protein goes in the main FASTA, and when the other allele removes a stop the
read-through protein is written to a second FASTA.

Requires: Python 3 standard library only

Usage:
  python3 bin/select_snp_proteins.py FILTERED.tsv REFERENCE.gff3 REFERENCE.faa \
      REFERENCE.fasta OUTPUT_PREFIX

Outputs:
  <prefix>.faa                 Reference protein sequences of the selected genes
  <prefix>.metadata.tsv        one row per protein
  <prefix>.readthrough.faa     read-through proteins for stop-loss cases

Example:
  python3 bin/select_snp_proteins.py \
      results/ranked_snps/ranked_snps.filtered.tsv \
      results/reference_annotation/Reference.gff3 \
      results/reference_annotation/Reference.faa \
      data/reference/Reference.fasta \
      results/snp_proteins/snp_proteins
"""

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
    records, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name:
            records[name].append(line)
    return {k: "".join(v) for k, v in records.items()}


def read_cds(path):
    cds = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        attrs = {k: urllib.parse.unquote(v) for k, v in
                 (kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)}
        cds.append({"locus_tag": attrs["locus_tag"], "start": int(f[3]), "end": int(f[4]),
                    "strand": f[6], "product": attrs.get("product", ""),
                    "gene": attrs.get("Name", "")})
    cds.sort(key=lambda c: c["start"])
    return cds


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq):
    protein = ""
    for i in range(0, len(seq) - 2, 3):
        aa = CODON_TABLE.get(seq[i:i + 3], "X")
        if aa == "*":
            break
        protein += aa
    return protein


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    filtered_path, gff_path, faa_path, fasta_path, prefix = sys.argv[1:]

    genome = "".join(read_fasta(fasta_path).values()).upper()
    proteins = read_fasta(faa_path)
    cds = read_cds(gff_path)
    by_tag = {c["locus_tag"]: c for c in cds}
    snps = list(csv.DictReader(open(filtered_path), delimiter="\t"))

    selected = {}  # locus tag -> {reasons, snps}

    def add(tag, reason, snp, detail=""):
        entry = selected.setdefault(tag, {"reasons": set(), "snps": []})
        entry["reasons"].add(reason)
        entry["snps"].append((snp, reason, detail))

    for snp in snps:
        position = int(snp["aln_col"])
        if snp["feature_type"] == "CDS" and snp["effect"] in ("missense", "nonsense"):
            add(snp["locus_tag"], snp["effect"], snp)
        elif snp["feature_type"] == "intergenic":
            left = [c for c in cds if c["end"] < position]
            right = [c for c in cds if c["start"] > position]
            if left:
                add(left[-1]["locus_tag"], "adjacent_to_intergenic_snp", snp,
                    f"{position - left[-1]['end']} bp after gene end")
            if right:
                add(right[0]["locus_tag"], "adjacent_to_intergenic_snp", snp,
                    f"{right[0]['start'] - position} bp before gene start")

    rows, readthrough = [], []
    with open(f"{prefix}.faa", "w") as faa:
        for tag in sorted(selected, key=lambda t: by_tag[t]["start"]):
            gene, entry = by_tag[tag], selected[tag]
            protein = proteins[tag]
            faa.write(f">{tag} {gene['product']}\n")
            for i in range(0, len(protein), 80):
                faa.write(protein[i:i + 80] + "\n")

            changes, notes, terminal_stop = [], [], ""
            for snp, reason, detail in entry["snps"]:
                position = int(snp["aln_col"])
                if reason in ("missense", "nonsense"):
                    if gene["strand"] == "+":
                        codon_number = (position - gene["start"]) // 3 + 1
                    else:
                        codon_number = (gene["end"] - position) // 3 + 1
                    changes.append(f"{snp['other_aa']}{codon_number}{snp['enriched_aa']}"
                                   f" ({snp['other_allele']}>{snp['enriched_allele']} at {position})")
                    if reason == "nonsense":
                        is_last = codon_number == (gene["end"] - gene["start"] + 1) // 3
                        terminal_stop = "yes" if is_last else "no"
                        if is_last and snp["enriched_aa"] == "*":
                            # Reference ends here; the other allele reads through.
                            mutated = genome[:position - 1] + snp["other_allele"] + genome[position:]
                            if gene["strand"] == "+":
                                extended = translate(mutated[gene["start"] - 1:gene["end"] + 3000])
                            else:
                                extended = translate(reverse_complement(
                                    mutated[max(0, gene["start"] - 1 - 3000):gene["end"]]))
                            readthrough.append((tag, extended))
                            notes.append(f"reference stop codon is lost by the other allele; "
                                         f"read-through protein {len(extended)} aa vs {len(protein)} aa")
                else:
                    changes.append(f"intergenic SNP {position} ({detail})")

            rows.append({
                "protein_id": tag, "contig": "Reference", "start": gene["start"], "end": gene["end"],
                "strand": gene["strand"], "length_aa": len(protein),
                "bakta_light_product": gene["product"], "bakta_light_gene": gene["gene"],
                "selection_reason": ";".join(sorted(entry["reasons"])),
                "snp_ranks": ",".join(sorted({s["rank"] for s, _, _ in entry["snps"]}, key=int)),
                "snp_positions": ",".join(sorted({s["aln_col"] for s, _, _ in entry["snps"]}, key=int)),
                "changes": "; ".join(changes),
                "snp_at_terminal_stop_codon": terminal_stop,
                "notes": "; ".join(notes),
            })

    with open(f"{prefix}.readthrough.faa", "w") as out:
        for tag, protein in readthrough:
            out.write(f">{tag}_readthrough\n")
            for i in range(0, len(protein), 80):
                out.write(protein[i:i + 80] + "\n")

    with open(f"{prefix}.metadata.tsv", "w") as out:
        columns = list(rows[0].keys())
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in columns) + "\n")
    print(f"{len(rows)} proteins selected from {len(snps)} filtered SNPs")
    print(f"  with SNP in CDS: {sum(1 for r in rows if 'missense' in r['selection_reason'] or 'nonsense' in r['selection_reason'])}"
          f", adjacent to intergenic SNPs: {sum(1 for r in rows if 'adjacent' in r['selection_reason'])}")
    print(f"  read-through proteins written: {len(readthrough)}")


if __name__ == "__main__":
    main()
