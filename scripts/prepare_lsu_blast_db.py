#!/usr/bin/env python3
"""Write a BLAST-friendly LSU FASTA and a map back to original FASTA headers."""

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--fasta-out", required=True, type=Path)
    parser.add_argument("--id-map-out", required=True, type=Path)
    args = parser.parse_args()

    args.fasta_out.parent.mkdir(parents=True, exist_ok=True)
    args.id_map_out.parent.mkdir(parents=True, exist_ok=True)

    record_count = 0
    with args.input.open(encoding="utf-8") as source, \
            args.fasta_out.open("w", encoding="utf-8", newline="\n") as fasta_out, \
            args.id_map_out.open("w", encoding="utf-8", newline="") as map_out:
        mapping = csv.writer(map_out, delimiter="\t", lineterminator="\n")
        mapping.writerow(("blast_id", "original_header"))

        header = None
        sequence = []

        def write_record() -> None:
            nonlocal record_count
            if header is None:
                return
            seq = "".join(sequence)
            if not seq:
                raise ValueError(f"Empty sequence for FASTA header: {header}")
            record_count += 1
            blast_id = f"lcl|symLSU{record_count:04d}"
            fasta_out.write(f">{blast_id} {header}\n")
            for offset in range(0, len(seq), 80):
                fasta_out.write(seq[offset:offset + 80] + "\n")
            mapping.writerow((blast_id, header))

        for line_number, raw_line in enumerate(source, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                write_record()
                header = line[1:]
                if not header:
                    raise ValueError(f"Empty FASTA header on line {line_number}")
                sequence = []
            elif header is None:
                raise ValueError(f"Sequence data before first FASTA header on line {line_number}")
            else:
                sequence.append("".join(line.split()))

        write_record()

    if record_count == 0:
        raise ValueError(f"No FASTA records found in {args.input}")

    print(f"Wrote {record_count} records to {args.fasta_out}")
    print(f"Wrote original-header map to {args.id_map_out}")


if __name__ == "__main__":
    main()
