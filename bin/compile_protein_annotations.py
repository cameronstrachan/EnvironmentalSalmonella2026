#!/usr/bin/env python3
"""Compile annotations of the SNP-associated proteins into one table.

One row per protein (plus the read-through variants, which only have eggNOG
results). Columns come from:
  metadata   selection reason, SNP positions and amino-acid changes
  bakta_light / bakta_full   Bakta JSON of the Reference genome with the light
             and full databases. The two runs can number locus tags
             differently, so genes are matched by start, stop and strand.
  eggnog     eggNOG-mapper annotations (the .emapper.annotations file)
  prodigal   Prodigal gene calls on the Reference genome, matched to the Bakta
             CDS by strand and stop position. Prodigal predicts genes from DNA,
             so it is used to check the gene model (start codon, partial,
             same protein) rather than to assign a function.
best_product is the first informative product (not "hypothetical protein") of
Bakta full, eggNOG description, Bakta light.

Requires: Python 3 standard library only

Usage:
  python3 bin/compile_protein_annotations.py METADATA.tsv BAKTA_LIGHT.json \
      BAKTA_FULL.json EGGNOG.emapper.annotations PRODIGAL.gff PRODIGAL.faa \
      OUTPUT.tsv

Example:
  python3 bin/compile_protein_annotations.py \
      results/snp_proteins/snp_proteins.metadata.tsv \
      results/reference_annotation/Reference.json \
      results/reference_annotation_full/Reference.json \
      results/snp_proteins/eggnog/snp_proteins.emapper.annotations \
      results/prodigal/Reference.gff results/prodigal/Reference.faa \
      results/snp_proteins/snp_protein_annotations.tsv
"""

import csv
import json
import sys

HYPOTHETICAL = {"", "hypothetical protein", "putative protein", "uncharacterized protein"}
XREF_PREFIXES = ("EC", "GO", "COG", "KEGG", "UniRef", "RefSeq", "NCBIProtein",
                 "BlastRules", "VFDB", "PFAM", "CDD", "Pfam", "ISfinder", "MOB")


def load_bakta(path):
    """Map (start, stop, strand) -> CDS feature dict from a Bakta JSON file."""
    features = json.load(open(path))["features"]
    return {(f["start"], f["stop"], f["strand"]): f for f in features if f["type"] == "cds"}


def xrefs_by_prefix(feature):
    xrefs = list(feature.get("db_xrefs", []))
    for key in ("psc", "pscc"):
        xrefs += feature.get(key, {}).get("db_xrefs", [])
    grouped = {}
    for x in xrefs:
        prefix, _, value = x.partition(":")
        if prefix in XREF_PREFIXES:
            grouped.setdefault(prefix, [])
            if value not in grouped[prefix]:
                grouped[prefix].append(value)
    return grouped


def bakta_columns(feature, name):
    if feature is None:
        return {f"{name}_{k}": "" for k in ("gene", "product", "ec", "go", "cog", "kegg",
                                            "uniref", "other_xrefs", "locus_tag")}
    grouped = xrefs_by_prefix(feature)
    uniref = grouped.pop("UniRef", [])
    other = [f"{k}:{v}" for k, vals in grouped.items() if k not in ("EC", "GO", "COG", "KEGG")
             for v in vals]
    return {
        f"{name}_locus_tag": feature.get("locus", ""),
        f"{name}_gene": feature.get("gene") or "",
        f"{name}_product": feature.get("product") or "",
        f"{name}_ec": ",".join(grouped.get("EC", [])),
        f"{name}_go": ",".join(grouped.get("GO", [])),
        f"{name}_cog": ",".join(grouped.get("COG", [])),
        f"{name}_kegg": ",".join(grouped.get("KEGG", [])),
        f"{name}_uniref": ",".join(uniref),
        f"{name}_other_xrefs": ",".join(other),
    }


def read_eggnog(path):
    rows, header = {}, None
    for line in open(path):
        if line.startswith("##") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if line.startswith("#query"):
            header = [h.lstrip("#") for h in fields]
            continue
        if header:
            rows[fields[0]] = dict(zip(header, fields))
    return rows


def read_prodigal(gff_path, faa_path):
    """Map (strand, stop) -> dict with start, end, attributes, protein."""
    proteins, name = {}, None
    for line in open(faa_path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            proteins[name] = ""
        else:
            proteins[name] += line
    genes = {}
    for line in open(gff_path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        start, end, strand = int(f[3]), int(f[4]), f[6]
        number = attrs["ID"].split("_")[1]
        protein = proteins.get(f"{f[0]}_{number}", "").rstrip("*")
        stop = end if strand == "+" else start
        genes[(strand, stop)] = {"start": start, "end": end, "attrs": attrs, "protein": protein}
    return genes


def main():
    if len(sys.argv) != 8:
        sys.exit(__doc__)
    meta_path, light_path, full_path, eggnog_path, pgff, pfaa, out_path = sys.argv[1:]

    metadata = list(csv.DictReader(open(meta_path), delimiter="\t"))
    light, full = load_bakta(light_path), load_bakta(full_path)
    eggnog = read_eggnog(eggnog_path)
    prodigal = read_prodigal(pgff, pfaa)

    egg_columns = ["seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
                   "COG_category", "Description", "Preferred_name", "GOs", "EC",
                   "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction",
                   "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"]

    def egg_values(query):
        row = eggnog.get(query, {})
        return {f"eggnog_{c}": ("" if row.get(c, "-") == "-" else row.get(c, "")) for c in egg_columns}

    rows = []
    for meta in metadata:
        key = (int(meta["start"]), int(meta["end"]), meta["strand"])
        light_feature, full_feature = light.get(key), full.get(key)
        row = dict(meta)
        row.update(bakta_columns(light_feature, "bakta_light"))
        row.update(bakta_columns(full_feature, "bakta_full"))
        row.update(egg_values(meta["protein_id"]))

        stop = int(meta["end"]) if meta["strand"] == "+" else int(meta["start"])
        gene = prodigal.get((meta["strand"], stop))
        if gene:
            same_start = (gene["start"], gene["end"]) == (int(meta["start"]), int(meta["end"]))
            bakta_protein = light_feature["aa"] if light_feature else ""
            row.update({
                "prodigal_start": gene["start"], "prodigal_end": gene["end"],
                "prodigal_same_coordinates": "yes" if same_start else "no",
                "prodigal_length_aa": len(gene["protein"]),
                "prodigal_same_protein": "yes" if gene["protein"] == bakta_protein else "no",
                "prodigal_start_type": gene["attrs"].get("start_type", ""),
                "prodigal_partial": gene["attrs"].get("partial", ""),
                "prodigal_rbs_motif": gene["attrs"].get("rbs_motif", ""),
                "prodigal_confidence": gene["attrs"].get("conf", ""),
            })
        else:
            row.update({"prodigal_same_coordinates": "no gene call with this stop"})

        candidates = [("bakta_full", row["bakta_full_product"]),
                      ("eggnog", row.get("eggnog_Description", "")),
                      ("bakta_light", row["bakta_light_product"])]
        best = next(((s, p) for s, p in candidates if p.strip().lower() not in HYPOTHETICAL),
                    ("none", "hypothetical protein"))
        row["best_product_source"], row["best_product"] = best
        rows.append(row)

        # Read-through variant of a stop-loss protein: eggNOG only.
        read_through = eggnog.get(f"{meta['protein_id']}_readthrough")
        if read_through:
            variant = {"protein_id": f"{meta['protein_id']}_readthrough",
                       "selection_reason": "read-through variant of " + meta["protein_id"],
                       "snp_positions": meta["snp_positions"],
                       "notes": meta["notes"]}
            variant.update(egg_values(f"{meta['protein_id']}_readthrough"))
            variant["best_product_source"] = "eggnog"
            variant["best_product"] = variant.get("eggnog_Description", "")
            rows.append(variant)

    columns = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with open(out_path, "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")
    print(f"wrote {len(rows)} rows x {len(columns)} columns to {out_path}")
    informative = sum(1 for r in rows if r["best_product_source"] != "none")
    print(f"{informative} of {len(rows)} rows have an informative product")


if __name__ == "__main__":
    main()
