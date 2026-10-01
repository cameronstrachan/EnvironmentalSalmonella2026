#!/usr/bin/env python3
"""Filter the annotated ranked SNP table.

Rows are dropped when
  - effect is "synonymous" (a synonymous SNP inside a CDS; intergenic and
    other non-CDS sites are kept), or
  - sensitivity is below --min-sensitivity, or
  - specificity is below --min-specificity.

Requires: Python 3 standard library only

Usage:
  python3 bin/filter_ranked_snps.py INPUT.tsv OUTPUT.tsv \
      [--min-sensitivity 1.0] [--min-specificity 0.9]

Example:
  python3 bin/filter_ranked_snps.py \
      results/ranked_snps/ranked_snps.annotated.tsv \
      results/ranked_snps/ranked_snps.filtered.tsv
"""

import csv
import sys


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        sys.exit(__doc__)
    in_path, out_path = args[:2]
    min_sens, min_spec = 1.0, 0.9
    rest = args[2:]
    while rest:
        flag = rest.pop(0)
        if flag == "--min-sensitivity":
            min_sens = float(rest.pop(0))
        elif flag == "--min-specificity":
            min_spec = float(rest.pop(0))
        else:
            sys.exit(f"unknown option {flag}\n{__doc__}")

    rows = list(csv.DictReader(open(in_path), delimiter="\t"))
    kept = [r for r in rows
            if r["effect"] != "synonymous"
            and float(r["sensitivity"]) >= min_sens
            and float(r["specificity"]) >= min_spec]
    with open(out_path, "w") as out:
        out.write("\t".join(rows[0].keys()) + "\n")
        for r in kept:
            out.write("\t".join(r.values()) + "\n")
    print(f"{len(rows)} rows in, {len(kept)} kept "
          f"(sensitivity >= {min_sens}, specificity >= {min_spec}, no synonymous)")


if __name__ == "__main__":
    main()
