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
        description="Create a compact QC summary from pairtools stats output."
    )
    parser.add_argument(
        "-p",
        "--pairtools-stats",
        required=True,
        help="Pairtools stats output file, for example stats/SAMPLE.pairtools.stats.txt",
    )
    parser.add_argument(
        "--format",
        choices=("text", "tsv"),
        default="text",
        help="Output format. Default: text (legacy-compatible).",
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


def percent_value(numerator, denominator):
    if denominator == 0:
        return "NA"
    return str(round(numerator * 100.0 / denominator, 2))


def percent(numerator, denominator):
    value = percent_value(numerator, denominator)
    if value == "NA":
        return value
    return f"{value}%"


def print_table(table):
    widths = [max(len(row[index]) for row in table) for index in range(3)]
    for row in table:
        print("  ".join(value.ljust(widths[index]) for index, value in enumerate(row)))


def qc_rows(stats):
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

    return [
        ("total_read_pairs", total_reads, "100" if total_reads else "NA", "Total Read Pairs", total_reads),
        ("unmapped_read_pairs", unmapped_reads, percent_value(unmapped_reads, total_reads), "Unmapped Read Pairs", unmapped_reads),
        ("mapped_read_pairs", mapped_reads, percent_value(mapped_reads, total_reads), "Mapped Read Pairs", mapped_reads),
        ("pcr_dup_read_pairs", dup_reads, percent_value(dup_reads, total_reads), "PCR Dup Read Pairs", dup_reads),
        ("no_dup_read_pairs", nodup_reads, percent_value(nodup_reads, total_reads), "No-Dup Read Pairs", nodup_reads),
        ("no_dup_cis_read_pairs", cis_reads, percent_value(cis_reads, nodup_reads), "No-Dup Cis Read Pairs", cis_reads),
        ("no_dup_trans_read_pairs", trans_reads, percent_value(trans_reads, nodup_reads), "No-Dup Trans Read Pairs", trans_reads),
        (
            "no_dup_valid_read_pairs_cis_ge_1kb_plus_trans",
            valid_read_pairs,
            percent_value(valid_read_pairs, nodup_reads),
            "No-Dup Valid Read Pairs (cis >= 1kb + trans)",
            valid_read_pairs,
        ),
        ("no_dup_cis_read_pairs_lt_1kb", cis_lt1kb, percent_value(cis_lt1kb, nodup_reads), "No-Dup Cis Read Pairs < 1kb", cis_lt1kb),
        ("no_dup_cis_read_pairs_ge_1kb", cis_gt1kb, percent_value(cis_gt1kb, nodup_reads), "No-Dup Cis Read Pairs >= 1kb", cis_gt1kb),
        ("no_dup_cis_read_pairs_ge_10kb", cis_gt10kb, percent_value(cis_gt10kb, nodup_reads), "No-Dup Cis Read Pairs >= 10kb", cis_gt10kb),
    ]


def print_text(rows):
    table = []
    for _, _, pct, label, value in rows:
        text_pct = pct if pct == "NA" else f"{pct}%"
        table.append([label, format(value, ",d"), text_pct])
    print_table(table)


def print_tsv(rows):
    print("metric\tvalue\tpercent")
    for metric, value, pct, _, _ in rows:
        print(f"{metric}\t{value}\t{pct}")


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

    rows = qc_rows(stats)
    if args.format == "tsv":
        print_tsv(rows)
    else:
        print_text(rows)


if __name__ == "__main__":
    main()
