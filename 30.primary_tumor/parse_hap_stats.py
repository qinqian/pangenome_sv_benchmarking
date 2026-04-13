#!/usr/bin/env python3

import argparse
import re
import sys
from pathlib import Path


def parse_stats_file(path: Path):
    sz = None
    n50 = None

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            fields = line.strip().split()
            if not fields:
                continue

            if fields[0] == "SZ" and len(fields) >= 2:
                sz = fields[1]
            elif fields[0] == "NL" and len(fields) >= 3 and fields[1] == "50":
                n50 = fields[2]

    return sz, n50


def split_sample_and_hap(filename: str):
    """
    Examples supported:
      sample1_hap1.txt
      sample1.hap2.stats
      sample1-hap1-report.txt
      sample1hap2.txt

    Returns:
      (sample, hap) where hap is 'hap1' or 'hap2'

    Raises:
      ValueError if hap1/hap2 cannot be identified.
    """
    name = Path(filename).name

    m = re.search(r"^(.*?)[._-]?(hap[12])(?:[._-].*|$)", name, flags=re.IGNORECASE)
    if not m:
        raise ValueError(f"Cannot determine haplotype from filename: {filename}")

    sample = m.group(1).rstrip("._-")
    hap = m.group(2).lower()

    if not sample:
        raise ValueError(f"Cannot determine sample name from filename: {filename}")

    return sample, hap


def build_rows(files):
    rows = {}

    for file in files:
        path = Path(file)
        sample, hap = split_sample_and_hap(path.name)
        sz, n50 = parse_stats_file(path)

        if sample not in rows:
            rows[sample] = {
                "hap1_sz": "",
                "hap1_n50": "",
                "hap2_sz": "",
                "hap2_n50": "",
            }

        if hap == "hap1":
            rows[sample]["hap1_sz"] = sz or ""
            rows[sample]["hap1_n50"] = n50 or ""
        elif hap == "hap2":
            rows[sample]["hap2_sz"] = sz or ""
            rows[sample]["hap2_n50"] = n50 or ""

    return rows


def write_output(rows, out_handle):
    header = ["sample", "hap1_sz", "hap1_n50", "hap2_sz", "hap2_n50"]
    out_handle.write("\t".join(header) + "\n")

    for sample in sorted(rows):
        r = rows[sample]
        out_handle.write(
            "\t".join(
                [
                    sample,
                    str(r["hap1_sz"]),
                    str(r["hap1_n50"]),
                    str(r["hap2_sz"]),
                    str(r["hap2_n50"]),
                ]
            )
            + "\n"
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse assembly stats files and merge hap1/hap2 SZ and N50 into one table."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input stats text files (hap1/hap2 identified from filenames).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Output TSV file. Default: stdout",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    try:
        rows = build_rows(args.inputs)

        if args.output == "-":
            write_output(rows, sys.stdout)
        else:
            with open(args.output, "w", encoding="utf-8") as out:
                write_output(rows, out)

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
