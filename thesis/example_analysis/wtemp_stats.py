#!/usr/bin/env python3
"""Mean and blocked error for KLFT temporal Wilson-loop output."""

import argparse
import sys

from stats_common import blocked_error, mean


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Read a KLFT temporal Wilson-loop file and print mean and blocked "
            "error for each (L, T). It also supports gradient-flow temporal "
            "Wilson-loop files with columns 'conf_id t_over_a2 L T W_temp'. "
            "Rows with step/conf_id < THERM are discarded."
        )
    )
    parser.add_argument("file", help="temporal Wilson-loop output file")
    parser.add_argument("THERM", type=int, help="thermalization step cutoff")
    parser.add_argument("BLOCK", type=int, help="blocking size in measurements")
    return parser.parse_args()


def parse_loop_size(token, filename, line_number, name):
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


def read_wtemp_values(filename, therm):
    values_by_loop = {}
    saw_standard = False
    saw_gradient_flow = False
    with open(filename, "r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.split()
            if len(fields) not in (4, 5):
                raise ValueError(
                    f"{filename}:{line_number}: expected 4 columns for normal "
                    "W_temp or 5 columns for gradient-flow W_temp"
                )

            try:
                conf_id = int(fields[0])
                if len(fields) == 4:
                    saw_standard = True
                    t_over_a2 = None
                    l_index = 1
                    w_index = 3
                else:
                    saw_gradient_flow = True
                    t_over_a2 = float(fields[1])
                    l_index = 2
                    w_index = 4
                wilson_loop = float(fields[w_index])
            except ValueError as exc:
                raise ValueError(
                    f"{filename}:{line_number}: could not parse W_temp row"
                ) from exc

            if conf_id < therm:
                continue

            loop_size = (
                t_over_a2,
                parse_loop_size(fields[l_index], filename, line_number, "L"),
                parse_loop_size(fields[l_index + 1], filename, line_number, "T"),
            )
            values_by_loop.setdefault(loop_size, []).append(wilson_loop)

    if saw_standard and saw_gradient_flow:
        raise ValueError(f"{filename}: mixed normal and gradient-flow W_temp formats")
    return values_by_loop


def main():
    args = parse_args()
    if args.THERM < 0:
        sys.exit("THERM must be non-negative")
    if args.BLOCK < 1:
        sys.exit("BLOCK must be >= 1")

    try:
        values_by_loop = read_wtemp_values(args.file, args.THERM)
    except OSError as exc:
        sys.exit(f"failed to read {args.file}: {exc}")
    except ValueError as exc:
        sys.exit(str(exc))

    if not values_by_loop:
        sys.exit("no temporal Wilson-loop measurements remain after thermalization cut")

    print("# t_over_a2 L T n_measurements block_size n_blocks tail mean error")
    for loop_size in sorted(values_by_loop):
        values = values_by_loop[loop_size]
        value_mean = mean(values)
        error, nblocks = blocked_error(values, args.BLOCK)
        tail = len(values) - nblocks * args.BLOCK
        t_over_a2 = "-" if loop_size[0] is None else f"{loop_size[0]:.12g}"
        print(
            f"{t_over_a2} {loop_size[1]} {loop_size[2]} "
            f"{len(values)} {args.BLOCK} {nblocks} {tail} "
            f"{value_mean:.12g} {error:.12g}"
        )


if __name__ == "__main__":
    main()
