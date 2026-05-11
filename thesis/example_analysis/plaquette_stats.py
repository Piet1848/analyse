#!/usr/bin/env python3
"""Mean and blocked error for KLFT plaquette output."""

import argparse
import sys

from stats_common import blocked_error, mean


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Read a KLFT plaquette file and print mean and blocked error. "
            "Rows with step < THERM are discarded."
        )
    )
    parser.add_argument("file", help="plaquette output file")
    parser.add_argument("THERM", type=int, help="thermalization step cutoff")
    parser.add_argument("BLOCK", type=int, help="blocking size in measurements")
    return parser.parse_args()


def read_plaquette_values(filename, therm):
    values = []
    with open(filename, "r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.split()
            if len(fields) < 2:
                raise ValueError(
                    f"{filename}:{line_number}: expected at least 2 columns"
                )

            try:
                step = int(fields[0])
                plaquette = float(fields[1])
            except ValueError as exc:
                raise ValueError(
                    f"{filename}:{line_number}: could not parse step/plaquette"
                ) from exc

            if step >= therm:
                values.append(plaquette)

    return values


def main():
    args = parse_args()
    if args.THERM < 0:
        sys.exit("THERM must be non-negative")
    if args.BLOCK < 1:
        sys.exit("BLOCK must be >= 1")

    try:
        values = read_plaquette_values(args.file, args.THERM)
    except OSError as exc:
        sys.exit(f"failed to read {args.file}: {exc}")
    except ValueError as exc:
        sys.exit(str(exc))

    if not values:
        sys.exit("no plaquette measurements remain after thermalization cut")

    value_mean = mean(values)
    error, nblocks = blocked_error(values, args.BLOCK)
    tail = len(values) - nblocks * args.BLOCK

    print("# observable n_measurements block_size n_blocks tail mean error")
    print(
        "plaquette "
        f"{len(values)} {args.BLOCK} {nblocks} {tail} "
        f"{value_mean:.12g} {error:.12g}"
    )


if __name__ == "__main__":
    main()
