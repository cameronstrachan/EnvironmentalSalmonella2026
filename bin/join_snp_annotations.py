#!/usr/bin/env python3
"""Join the filtered ranked SNPs with the compiled protein annotations.

One row per SNP. For a SNP inside a CDS, the annotation columns describe that
gene. For an intergenic SNP, the columns describe the nearest CDS on each side
(gene_1 = upstream on the genome, gene_2 = downstream), with the distance in bp.
Proteins are taken from the metadata table's "changes" field ("intergenic SNP
<pos> (<n> bp after gene end / before gene start)").

Requires: Python 3 standard library only

Usage:
  python3 bin/join_snp_annotations.py FILTERED_SNPS.tsv PROTEIN_ANNOTATIONS.tsv OUTPUT.tsv

Example:
  python3 bin/join_snp_annotations.py \
      results/ranked_snps/ranked_snps.filtered.tsv \
      results/snp_proteins/snp_protein_annotations.tsv \
      results/ranked_snps/ranked_snps.filtered.annotated_full.tsv
"""

import csv
import re
import sys

SNP_COLUMNS = ["rank", "aln_col", "ref_allele", "enriched_allele", "other_allele",
               "plasmid_with_allele", "plasmid_total", "nonplasmid_with_allele",
               "nonplasmid_total", "sensitivity", "specificity", "youden_j",
               "fisher_p", "fisher_q", "feature_type", "effect"]
GENE_COLUMNS = ["protein_id", "start", "end", "strand", "best_product",
                "best_product_source", "bakta_full_gene", "bakta_full_product",
                "bakta_light_product", "eggnog_Preferred_name", "eggnog_Description",
                "eggnog_COG_category", "eggnog_KEGG_ko", "eggnog_PFAMs", "bakta_full_ec",
                "bakta_full_go", "bakta_full_uniref", "changes"]


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    snp_path, ann_path, out_path = sys.argv[1:]
    snps = list(csv.DictReader(open(snp_path), delimiter="\t"))
    proteins = {r["protein_id"]: r for r in csv.DictReader(open(ann_path), delimiter="\t")}

    # Which proteins sit next to which intergenic SNP, and how far away.
    neighbours = {}
    for protein in proteins.values():
        for match in re.finditer(r"intergenic SNP (\d+) \((\d+) bp (after gene end|before gene start)\)",
                                 protein.get("changes", "")):
            position, distance, where = match.groups()
            side = 1 if where == "after gene end" else 2
            neighbours.setdefault(position, {})[side] = (protein, distance)

    columns = SNP_COLUMNS + ["snp_gene_or_neighbour"]
    columns += [f"gene_{c}" for c in GENE_COLUMNS]
    columns += [f"gene_2_{c}" for c in GENE_COLUMNS] + ["gene_1_distance_bp", "gene_2_distance_bp"]
    with open(out_path, "w") as out:
        out.write("\t".join(columns) + "\n")
        for snp in snps:
            row = {c: snp[c] for c in SNP_COLUMNS}
            if snp["feature_type"] == "CDS":
                gene = proteins[snp["locus_tag"]]
                row["snp_gene_or_neighbour"] = "in gene"
                row.update({f"gene_{c}": gene.get(c, "") for c in GENE_COLUMNS})
            else:
                found = neighbours.get(snp["aln_col"], {})
                row["snp_gene_or_neighbour"] = "intergenic, flanking genes"
                for side, prefix in ((1, "gene_"), (2, "gene_2_")):
                    if side in found:
                        protein, distance = found[side]
                        row.update({f"{prefix}{c}": protein.get(c, "") for c in GENE_COLUMNS})
                        row[f"gene_{side}_distance_bp"] = distance
            out.write("\t".join(row.get(c, "") for c in columns) + "\n")
    print(f"wrote {len(snps)} SNPs to {out_path}")


if __name__ == "__main__":
    main()
