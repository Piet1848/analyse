from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

FlowKey = Tuple[Optional[float], int, int]


class EmptyNumericDataError(ValueError):
    pass


class _NumericSeries:
    def __init__(self, values: np.ndarray):
        self.values = values


class _NumericTable:
    def __init__(self, columns: List[str], data: Dict[str, np.ndarray]):
        self.columns = columns
        self._data = data

    def __getitem__(self, name: str) -> _NumericSeries:
        return _NumericSeries(self._data[name])


class VariableData:
    def __init__(self, name: str):
        self.name = name
        self.value = None
        self.error = None
        self.bootstrap_samples: np.ndarray | None = None
        self.bootstrap_invalid_count = 0
        self.bootstrap_finite_fraction = None
        self.parameters: Dict[str, Any] = {}

    def set_value(self, value: Any, bootstrap_samples: Any = None, **params):
        self.value = value
        if bootstrap_samples is not None:
            boot_arr = np.asarray(bootstrap_samples, dtype=float)
            if boot_arr.ndim == 0:
                boot_arr = boot_arr.reshape(1)

            finite_boots = boot_arr[np.isfinite(boot_arr)]
            repaired_boots = boot_arr
            self.bootstrap_invalid_count = int((~np.isfinite(boot_arr)).sum())
            self.bootstrap_finite_fraction = (
                float(finite_boots.size / boot_arr.size) if boot_arr.size > 0 else None
            )

            # Treat isolated failed bootstrap replicas as missing values and
            # replace them with the central estimate so downstream error
            # propagation stays usable.
            invalid_mask = ~np.isfinite(repaired_boots)
            if repaired_boots.size > 0 and np.any(invalid_mask):
                replacement = None
                try:
                    value_float = float(value)
                except (TypeError, ValueError):
                    value_float = np.nan

                if np.isfinite(value_float):
                    replacement = value_float
                elif finite_boots.size > 0:
                    replacement = float(np.mean(finite_boots))

                if replacement is not None and np.isfinite(replacement):
                    repaired_boots = repaired_boots.copy()
                    repaired_boots[invalid_mask] = replacement

            self.bootstrap_samples = repaired_boots
            repaired_finite = repaired_boots[np.isfinite(repaired_boots)]
            if repaired_finite.size > 0:
                self.error = float(np.std(repaired_finite))
            else:
                self.error = None
        else:
            self.bootstrap_samples = None
            self.error = None
            self.bootstrap_invalid_count = 0
            self.bootstrap_finite_fraction = None
        self.parameters.update(params)

    def get(self) -> Any:
        return self.value

    def err(self) -> Any:
        return self.error

    def bootstrap(self) -> np.ndarray | None:
        return self.bootstrap_samples

    def __repr__(self):
        return f"VariableData(name={self.name}, value={self.value}, error={self.error})"


class ObservableData:
    def __init__(self, name: str, values=None):
        self.name = name
        if values is not None:
            self.values = np.asarray(values, dtype=np.float32)
        else:
            self.values = np.array([], dtype=np.float32)

    def append(self, value: float):
        self.values = np.append(self.values, np.float32(value))

    def extend(self, values):
        self.values = np.concatenate((self.values, np.asarray(values, dtype=np.float32)))

    def slice(self, indices: List[int]):
        self.values = self.values[indices]

    def slice_from(self, start: int):
        self.values = self.values[start:]

    def __repr__(self):
        sample = self.values[:5]
        return f"ObservableData(name={self.name}, values={sample}{'...' if len(self.values) > 5 else ''})"

    def get_bootstrap_samples(self, n_bootstrap: int, seed: int):
        rng = np.random.default_rng(seed)
        arr = np.array(self.values)
        if len(arr) == 0:
            for _ in range(n_bootstrap):
                yield ObservableData(self.name, [])
            return

        for _ in range(n_bootstrap):
            indices = rng.integers(0, len(arr), size=len(arr))
            yield ObservableData(self.name, arr[indices])


class FileData:
    def __init__(self, path: str):
        self.name = Path(path).stem
        self.path = Path(path)
        self.observables: List[ObservableData] = []

    def read_file(self):
        if not self.path.exists():
            raise FileNotFoundError(f"{self.path} not found")

        try:
            df = _read_numeric_table(self.path)
        except (EmptyNumericDataError, ValueError):
            return self

        self.observables = []
        for name in df.columns:
            self.observables.append(ObservableData(name, df[name].values))
        return self

    def align_lengths(self) -> int:
        if not self.observables:
            return 0

        min_len = min(len(obs.values) for obs in self.observables)
        for obs in self.observables:
            if len(obs.values) != min_len:
                obs.values = obs.values[:min_len]
        return min_len

    def remove_thermalization(self, min_step: int):
        if not self.observables:
            return
        self.align_lengths()

        step_obs = next((o for o in self.observables if o.name in ["# step", "step", "conf_id"]), None)
        if step_obs is None:
            return

        steps = np.asarray(step_obs.values)
        if len(steps) == 0:
            return

        monotonic = np.all(steps[1:] >= steps[:-1])

        if monotonic:
            start = int(np.searchsorted(steps, min_step, side='left'))
            if start <= 0:
                return
            if start >= len(steps):
                for obs in self.observables:
                    obs.values = np.array([], dtype=np.float32)
                return
            for obs in self.observables:
                obs.slice_from(start)
            return

        print(f"Warning: Non-monotonic steps in {self.path}, Stopping thermalization removal.")
        return

    def get(self, name: str) -> ObservableData:
        for obs in self.observables:
            if obs.name == name:
                return obs
        raise ValueError(f"Observable '{name}' not found in {self.path}")

    def get_bootstrap(self, n_bootstrap: int, seed: int):
        for i in range(n_bootstrap):
            bootstrap_file = FileData(str(self.path))
            bootstrap_file.observables = [
                next(iter(obs.get_bootstrap_samples(1, seed + i))) for obs in self.observables
            ]
            yield bootstrap_file

    def get_blocked_bootstrap(self, n_bootstrap: int, block_size: int, seed: int):
        step_obs = next((o for o in self.observables if o.name in ["# step", "step"]), None)
        if step_obs is None or block_size <= 1:
            yield from self.get_bootstrap(n_bootstrap, seed)
            return

        steps = np.array(step_obs.values)
        unique_steps = np.unique(steps)
        n_steps = len(unique_steps)

        step_to_indices = defaultdict(list)
        for idx, s in enumerate(steps):
            step_to_indices[s].append(idx)

        rng = np.random.default_rng(seed)
        block_size = int(block_size)
        if block_size < 1:
            raise ValueError("block_size must be at least 1")

        block_slices = [
            unique_steps[start_idx:start_idx + block_size]
            for start_idx in range(0, n_steps, block_size)
        ]
        n_blocks = len(block_slices)

        for _ in range(n_bootstrap):
            sampled_blocks = rng.integers(0, n_blocks, size=n_blocks)

            selected_row_indices = []
            for block_idx in sampled_blocks:
                block_steps = block_slices[int(block_idx)]
                for s in block_steps:
                    selected_row_indices.extend(step_to_indices[s])

            new_fd = FileData(str(self.path))
            new_fd.observables = []
            for obs in self.observables:
                arr = np.array(obs.values)
                new_vals = arr[selected_row_indices].tolist()
                new_fd.observables.append(ObservableData(obs.name, new_vals))

            yield new_fd


class CompactWilsonData(FileData):
    def __init__(
        self,
        path: str,
        pair_order: List[Tuple[int, int]] | None = None,
        wilson_by_pair: Dict[Tuple[int, int], np.ndarray] | None = None,
        flow_pair_order: List[FlowKey] | None = None,
        wilson_by_flow_pair: Dict[FlowKey, np.ndarray] | None = None,
    ):
        super().__init__(path)
        if wilson_by_flow_pair is None:
            wilson_by_flow_pair = {
                (None, int(r), int(t)): values
                for (r, t), values in (wilson_by_pair or {}).items()
            }
        if flow_pair_order is None:
            flow_pair_order = [
                (None, int(r), int(t))
                for r, t in (pair_order or list(wilson_by_pair or {}))
            ]

        self.flow_pair_order = [
            (_normalize_flow_time(flow_time), int(r), int(t))
            for flow_time, r, t in flow_pair_order
        ]
        self.wilson_by_flow_pair = {
            (_normalize_flow_time(flow_time), int(r), int(t)): np.asarray(values, dtype=np.float32)
            for (flow_time, r, t), values in wilson_by_flow_pair.items()
        }
        self.pair_order = [(r, t) for flow_time, r, t in self.flow_pair_order if flow_time is None]
        self.wilson_by_pair = {
            (r, t): values
            for (flow_time, r, t), values in self.wilson_by_flow_pair.items()
            if flow_time is None
        }
        first = next(iter(self.wilson_by_flow_pair.values()), np.array([], dtype=np.float32))
        self.n_configurations = int(len(first))

    def available_flow_times(self) -> list[Optional[float]]:
        return sorted(
            {flow_time for flow_time, _, _ in self.flow_pair_order},
            key=lambda value: (-1.0 if value is None else float(value)),
        )


def _normalize_flow_time(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        flow_time = float(value)
    except (TypeError, ValueError):
        return None
    return round(flow_time, 12)


def _read_numeric_table(path: Path) -> _NumericTable:
    delimiter: str | None = None
    header_tokens: list[str] | None = None
    rows: list[list[float]] = []
    expected_width: int | None = None

    with path.open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                tokens = stripped[1:].strip().split()
                if tokens and "=" not in stripped and all(_looks_like_column_name(token) for token in tokens):
                    header_tokens = tokens
                continue

            if delimiter is None:
                delimiter = "," if "," in stripped else None
            tokens = _split_table_line(stripped, delimiter)
            if not tokens:
                continue

            if header_tokens is None and all(_looks_like_column_name(token) for token in tokens):
                header_tokens = tokens
                continue

            try:
                values = [float(token) for token in tokens]
            except ValueError:
                continue

            if expected_width is None:
                expected_width = len(values)
                if header_tokens is None:
                    header_tokens = [f"col{idx}" for idx in range(expected_width)]
                elif len(header_tokens) != expected_width:
                    raise ValueError(
                        f"{path}: header has {len(header_tokens)} columns but first data row has {expected_width}"
                    )
            elif len(values) != expected_width:
                continue

            rows.append(values)

    if not rows or header_tokens is None:
        raise EmptyNumericDataError("no numeric data")

    matrix = np.asarray(rows, dtype=np.float32)
    data = {
        name: matrix[:, idx]
        for idx, name in enumerate(header_tokens)
    }
    return _NumericTable(header_tokens, data)


def _split_table_line(line: str, delimiter: str | None) -> list[str]:
    if delimiter == ",":
        return [token.strip() for token in next(csv.reader([line], skipinitialspace=True))]
    return line.split()


def _looks_like_column_name(token: str) -> bool:
    return any(ch.isalpha() for ch in token) and all(
        ch.isalnum() or ch in {"_", "/"} for ch in token
    )


class ExperimentData:
    def __init__(self, path: str):
        self.path = Path(path)
        self.data: Dict[str, List[FileData]] = self._get_data()

    def _get_data(self) -> Dict[str, List[FileData]]:
        data_dict: Dict[str, List[FileData]] = {}
        if not self.path.exists():
            return data_dict

        for file_path in self.path.iterdir():
            if file_path.is_file() and file_path.suffix == ".out" and not file_path.name.startswith("_"):
                fd = FileData(str(file_path))
                fd.read_file()
                data_dict.setdefault(file_path.stem, []).append(fd)
        return data_dict


def _infer_rows_per_configuration(
    steps: Optional[np.ndarray],
    obs_L: np.ndarray,
    obs_T: np.ndarray,
    flow_times: Optional[np.ndarray] = None,
) -> int:
    if steps is not None and len(steps) > 0:
        change_idx = np.flatnonzero(steps != steps[0])
        if change_idx.size > 0:
            return int(change_idx[0])
        return int(len(steps))

    seen: set[FlowKey] = set()
    flow_iter = flow_times if flow_times is not None else [None] * len(obs_L)
    for i, (flow_time, l_val, t_val) in enumerate(zip(flow_iter, obs_L, obs_T)):
        pair = (_normalize_flow_time(flow_time), int(l_val), int(t_val))
        if pair in seen:
            return i
        seen.add(pair)
    return int(len(obs_L))


def load_compact_wilson_file(path: str, min_step: int = 0) -> CompactWilsonData | None: #marker4
    fd = FileData(path)
    fd.read_file()
    fd.align_lengths()
    if min_step > 0:
        fd.remove_thermalization(min_step)
    fd.align_lengths() # should be redundant but ensures consistency after thermalization removal

    if not fd.observables:
        return None
    if min((len(obs.values) for obs in fd.observables), default=0) <= 0:
        return None

    try:
        obs_w = np.asarray(fd.get("W_temp").values, dtype=np.float32)
        obs_L = np.asarray(fd.get("L").values)
        obs_T = np.asarray(fd.get("T").values)
    except ValueError:
        print(f"Warning: Required observables missing in {path}, skipping W_temp loading.")
        return None

    flow_obs = next((o for o in fd.observables if o.name == "t_over_a2"), None)
    flow_times = np.asarray(flow_obs.values) if flow_obs is not None else None
    step_obs = next((o for o in fd.observables if o.name in ["# step", "step", "conf_id"]), None)
    steps = np.asarray(step_obs.values) if step_obs is not None else None
    rows_per_cfg = _infer_rows_per_configuration(steps, obs_L, obs_T, flow_times=flow_times)

    if rows_per_cfg <= 0 or len(obs_w) % rows_per_cfg != 0:
        print(f"Warning: Unable to determine rows per configuration for {path}, skipping W_temp loading.")
        return None

    if flow_times is None:
        flow_order = [(None, int(l), int(t)) for l, t in zip(obs_L[:rows_per_cfg], obs_T[:rows_per_cfg])]
    else:
        flow_order = [
            (_normalize_flow_time(flow_time), int(l), int(t))
            for flow_time, l, t in zip(flow_times[:rows_per_cfg], obs_L[:rows_per_cfg], obs_T[:rows_per_cfg])
        ]
    if len(set(flow_order)) != len(flow_order):
        print(f"Warning: Duplicate pairs found in {path}, skipping W_temp loading.")
        return None

    w_matrix = obs_w.reshape(-1, rows_per_cfg)
    wilson_by_flow_pair = {pair: w_matrix[:, idx] for idx, pair in enumerate(flow_order)}
    return CompactWilsonData(path, flow_pair_order=flow_order, wilson_by_flow_pair=wilson_by_flow_pair)


def combine_compact_wilson_data(
    files: Iterable[CompactWilsonData],
    source_name: str = "combined",
) -> CompactWilsonData | None:
    flow_pair_order: Optional[List[FlowKey]] = None
    common_pairs: Optional[set[FlowKey]] = None
    chunks_by_pair: Dict[FlowKey, List[np.ndarray]] = defaultdict(list)

    for fd in files:
        if not fd.wilson_by_flow_pair:
            continue
        current_pairs = list(fd.flow_pair_order)
        current_set = set(current_pairs)
        if flow_pair_order is None:
            flow_pair_order = current_pairs
            common_pairs = set(current_pairs)
        else:
            previous_common = set(common_pairs)
            common_pairs &= current_set
            for pair in previous_common - common_pairs:
                chunks_by_pair.pop(pair, None)

        active_pairs = common_pairs if common_pairs is not None else current_set
        for pair in current_pairs:
            if pair in active_pairs:
                chunks_by_pair[pair].append(fd.wilson_by_flow_pair[pair])

    if not flow_pair_order or not common_pairs:
        return None

    final_order = [pair for pair in flow_pair_order if pair in common_pairs]
    if not final_order:
        return None

    combined_pairs = {
        pair: np.concatenate(chunks_by_pair[pair]).astype(np.float32, copy=False)
        for pair in final_order
    }
    return CompactWilsonData(
        f"{source_name}.out",
        flow_pair_order=final_order,
        wilson_by_flow_pair=combined_pairs,
    )


def combine_file_data(
    files: Iterable[FileData],
    min_step: int = 0,
    source_name: str = "combined",
) -> FileData | None:
    ordered_names: Optional[List[str]] = None
    common_names: Optional[set[str]] = None
    chunks_by_name: Dict[str, List[np.ndarray]] = defaultdict(list)

    for fd in files:
        fd.align_lengths()
        if min_step > 0:
            fd.remove_thermalization(min_step)
        fd.align_lengths()
        if not fd.observables:
            continue

        chunk_len = min((len(obs.values) for obs in fd.observables), default=0)
        if chunk_len <= 0:
            continue

        names = [obs.name for obs in fd.observables]
        if ordered_names is None:
            ordered_names = names
            common_names = set(names)
        else:
            common_names &= set(names)

        for obs in fd.observables:
            chunks_by_name[obs.name].append(obs.values[:chunk_len])

    if not ordered_names or not common_names:
        return None

    final_names = [name for name in ordered_names if name in common_names]
    if not final_names:
        return None

    combined = FileData(f"{source_name}.out")
    combined.observables = [
        ObservableData(name, np.concatenate(chunks_by_name[name]).astype(np.float32, copy=False))
        for name in final_names
    ]

    if min((len(obs.values) for obs in combined.observables), default=0) <= 0:
        return None

    combined.align_lengths()
    return combined
