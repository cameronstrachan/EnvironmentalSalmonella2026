#!/usr/bin/env python3
"""Extract one sequence from an alignment and write it as an ungapped FASTA.

The core alignment has no gap columns, so a position in the extracted
sequence equals the 1-based alignment column.

Usage:
  python3 bin/extract_sequence.py ALIGNMENT.aln SEQUENCE_ID OUTPUT.fasta

Example:
  python3 bin/extract_sequence.py \
      data/alignments/infantis_core_alignment_noambiguous.aln Reference \
      data/reference/Reference.fasta
"""

import sys


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    aln_path, wanted, out_path = sys.argv[1:]

    parts, found = [], False
    with open(aln_path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if found:
                    break
                found = line[1:] == wanted
            elif found:
                parts.append(line)
    if not parts:
        sys.exit(f"sequence {wanted!r} not found in {aln_path}")

    seq = "".join(parts).replace("-", "")
    with open(out_path, "w") as out:
        out.write(f">{wanted}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i + 80] + "\n")
    print(f"wrote {wanted}: {len(seq)} bases to {out_path}")


if __name__ == "__main__":
    main()
