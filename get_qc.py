#!/usr/bin/env python3

import argparse
import sys


REQUIRED_KEYS = (
    "total",
    "total_unmapped",
    "total_mapped",
    "total_dups",
    "total_nodups",
    "cis",
    "trans",
    "cis_1kb+",
    "cis_10kb+",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a compact text QC summary from pairtools stats output."
    )
    parser.add_argument(
        "-p",
        "--pairtools-stats",
        required=True,
        help="Pairtools stats output file, for example stats/SAMPLE.txt",
    )
    return parser.parse_args()


def read_stats(path):
    output_dict = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            attrs = line.split()
            if len(attrs) >= 2:
                output_dict[attrs[0]] = attrs[1]
    return output_dict


def require_int(stats, key):
    try:
        return int(stats[key])
    except KeyError:
        print(f"Missing required pairtools stats key: {key}", file=sys.stderr)
        sys.exit(1)
    except ValueError:
        print(f"Pairtools stats key is not an integer: {key}={stats[key]}", file=sys.stderr)
        sys.exit(1)


def percent(numerator, denominator):
    if denominator == 0:
        return "NA"
    return f"{round(numerator * 100.0 / denominator, 2)}%"


def print_table(table):
    widths = [max(len(row[index]) for row in table) for index in range(3)]
    for row in table:
        print("  ".join(value.ljust(widths[index]) for index, value in enumerate(row)))


def main():
    args = parse_args()
    stats = read_stats(args.pairtools_stats)

    missing_keys = [key for key in REQUIRED_KEYS if key not in stats]
    if missing_keys:
        print(
            "Missing required pairtools stats keys: " + ", ".join(missing_keys),
            file=sys.stderr,
        )
        sys.exit(1)

    total_reads = require_int(stats, "total")
    unmapped_reads = require_int(stats, "total_unmapped")
    mapped_reads = require_int(stats, "total_mapped")
    dup_reads = require_int(stats, "total_dups")
    nodup_reads = require_int(stats, "total_nodups")
    cis_reads = require_int(stats, "cis")
    trans_reads = require_int(stats, "trans")
    cis_gt1kb = require_int(stats, "cis_1kb+")
    cis_gt10kb = require_int(stats, "cis_10kb+")

    cis_lt1kb = cis_reads - cis_gt1kb
    valid_read_pairs = trans_reads + cis_gt1kb

    table = [
        ["Total Read Pairs", format(total_reads, ",d"), "100%" if total_reads else "NA"],
        ["Unmapped Read Pairs", format(unmapped_reads, ",d"), percent(unmapped_reads, total_reads)],
        ["Mapped Read Pairs", format(mapped_reads, ",d"), percent(mapped_reads, total_reads)],
        ["PCR Dup Read Pairs", format(dup_reads, ",d"), percent(dup_reads, total_reads)],
        ["No-Dup Read Pairs", format(nodup_reads, ",d"), percent(nodup_reads, total_reads)],
        ["No-Dup Cis Read Pairs", format(cis_reads, ",d"), percent(cis_reads, nodup_reads)],
        ["No-Dup Trans Read Pairs", format(trans_reads, ",d"), percent(trans_reads, nodup_reads)],
        [
            "No-Dup Valid Read Pairs (cis >= 1kb + trans)",
            format(valid_read_pairs, ",d"),
            percent(valid_read_pairs, nodup_reads),
        ],
        ["No-Dup Cis Read Pairs < 1kb", format(cis_lt1kb, ",d"), percent(cis_lt1kb, nodup_reads)],
        ["No-Dup Cis Read Pairs >= 1kb", format(cis_gt1kb, ",d"), percent(cis_gt1kb, nodup_reads)],
        ["No-Dup Cis Read Pairs >= 10kb", format(cis_gt10kb, ",d"), percent(cis_gt10kb, nodup_reads)],
    ]
    print_table(table)


if __name__ == "__main__":
    main()
