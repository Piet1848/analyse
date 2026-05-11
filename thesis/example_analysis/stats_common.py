"""Shared statistics helpers for KLFT analysis scripts."""

import math


def mean(values):
    return sum(values) / len(values)


def blocked_error(values, block_size):
    nblocks = len(values) // block_size
    if nblocks < 2:
        return float("nan"), nblocks

    block_means = []
    for block_index in range(nblocks):
        start = block_index * block_size
        block = values[start : start + block_size]
        block_means.append(mean(block))

    block_mean = mean(block_means)
    variance = sum((value - block_mean) ** 2 for value in block_means)
    variance /= nblocks - 1
    return math.sqrt(variance / nblocks), nblocks
