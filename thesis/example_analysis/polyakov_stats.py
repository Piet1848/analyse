#!/usr/bin/env python3
"""Mean and blocked error for KLFT Polyakov-loop outputs."""

import argparse
import sys

from stats_common import blocked_error, mean


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Read a KLFT Polyakov-loop or Polyakov-correlator file and print "
            "mean and blocked error. Rows with step < THERM are discarded. "
            "Polyakov-loop rows are '# step RePolyakov ImPolyakov'; "
            "correlator rows are '# step R real imaginary'."
        )
    )
    parser.add_argument("file", help="Polyakov-loop or correlator output file")
    parser.add_argument("THERM", type=int, help="thermalization step cutoff")
    parser.add_argument("BLOCK", type=int, help="blocking size in measurements")
    return parser.parse_args()


def parse_int_token(token, filename, line_number, name):
    try:
        value = float(token)
    except ValueError as exc:
        raise ValueError(
            f"{filename}:{line_number}: could not parse {name}"
        ) from exc

    integer_value = int(value)
    if value != integer_value:
        raise ValueError(f"{filename}:{line_number}: {name} is not an integer")
    return integer_value


def read_polyakov_values(filename, therm):
    loop_values = []
    correlator_values = {}
    saw_loop = False
    saw_correlator = False

    with open(filename, "r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.split()
            if len(fields) == 3:
                saw_loop = True
                try:
                    step = int(fields[0])
                    real_part = float(fields[1])
                    imaginary_part = float(fields[2])
                except ValueError as exc:
                    raise ValueError(
                        f"{filename}:{line_number}: could not parse "
                        "step/RePolyakov/ImPolyakov"
                    ) from exc

                if step >= therm:
                    loop_values.append((real_part, imaginary_part))
            elif len(fields) == 4:
                saw_correlator = True
                try:
                    step = int(fields[0])
                    real_part = float(fields[2])
                    imaginary_part = float(fields[3])
                except ValueError as exc:
                    raise ValueError(
                        f"{filename}:{line_number}: could not parse "
                        "step/real/imaginary"
                    ) from exc

                if step < therm:
                    continue

                r_value = parse_int_token(fields[1], filename, line_number, "R")
                correlator_values.setdefault(r_value, []).append(
                    (real_part, imaginary_part)
                )
            else:
                raise ValueError(
                    f"{filename}:{line_number}: expected either 3 columns "
                    "for Polyakov loop or 4 columns for correlator"
                )

    if saw_loop and saw_correlator:
        raise ValueError(
            f"{filename}: mixed Polyakov-loop and correlator row formats"
        )
    if saw_loop:
        return "loop", loop_values
    if saw_correlator:
        return "correlator", correlator_values
    return "empty", []


def summarize_complex_pairs(values, block_size):
    real_values = [value[0] for value in values]
    imaginary_values = [value[1] for value in values]
    real_error, nblocks = blocked_error(real_values, block_size)
    imaginary_error, imag_nblocks = blocked_error(imaginary_values, block_size)
    if imag_nblocks != nblocks:
        raise RuntimeError("internal error: inconsistent block counts")

    tail = len(values) - nblocks * block_size
    return (
        len(values),
        nblocks,
        tail,
        mean(real_values),
        real_error,
        mean(imaginary_values),
        imaginary_error,
    )


def main():
    args = parse_args()
    if args.THERM < 0:
        sys.exit("THERM must be non-negative")
    if args.BLOCK < 1:
        sys.exit("BLOCK must be >= 1")

    try:
        file_kind, values = read_polyakov_values(args.file, args.THERM)
    except OSError as exc:
        sys.exit(f"failed to read {args.file}: {exc}")
    except ValueError as exc:
        sys.exit(str(exc))

    if file_kind == "empty":
        sys.exit("no Polyakov measurements found")

    print(
        "# observable R n_measurements block_size n_blocks tail "
        "real_mean real_error imaginary_mean imaginary_error"
    )
    if file_kind == "loop":
        if not values:
            sys.exit("no Polyakov-loop measurements remain after thermalization cut")
        nvalues, nblocks, tail, rmean, rerr, imean, ierr = summarize_complex_pairs(
            values, args.BLOCK
        )
        print(
            "polyakov_loop - "
            f"{nvalues} {args.BLOCK} {nblocks} {tail} "
            f"{rmean:.12g} {rerr:.12g} {imean:.12g} {ierr:.12g}"
        )
    else:
        if not values:
            sys.exit(
                "no Polyakov-correlator measurements remain after thermalization cut"
            )
        for r_value in sorted(values):
            entries = values[r_value]
            nvalues, nblocks, tail, rmean, rerr, imean, ierr = (
                summarize_complex_pairs(entries, args.BLOCK)
            )
            print(
                f"polyakov_correlator {r_value} "
                f"{nvalues} {args.BLOCK} {nblocks} {tail} "
                f"{rmean:.12g} {rerr:.12g} {imean:.12g} {ierr:.12g}"
            )


if __name__ == "__main__":
    main()
