from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List
from collections import defaultdict
import csv
import numpy as np


class VariableData:
    def __init__(self, name: str):
        self.name = name
        self.value = None
        self.error = None
        self.bootstrap_samples = None
        self.parameters: Dict[str, Any] = {}

    def set_value(self, value: Any, bootstrap_samples: List[Any] = None, **params):
        self.value = value
        if bootstrap_samples is not None:
            self.bootstrap_samples = bootstrap_samples
            self.error = np.std(bootstrap_samples)
        else:
            self.error = None
        self.parameters.update(params)

    def get(self) -> Any:
        return self.value

    def err(self) -> Any:
        return self.error
    
    def bootstrap(self) -> List[Any]:
        return self.bootstrap_samples
    
    def __repr__(self):
        return f"VariableData(name={self.name}, value={self.value}, error={self.error})"


class ObservableData:
    def __init__(self, name: str, values=None):
        self.name = name
        self.values = values if values is not None else []

    def append(self, value: float):
        if isinstance(self.values, np.ndarray):
            self.values = np.append(self.values, value)
        else:
            self.values.append(value)

    def extend(self, values):
        if isinstance(self.values, np.ndarray):
            self.values = np.concatenate((self.values, values))
        else:
            self.values.extend(values)
        
    def slice(self, indices: List[int]):
        """Keep only values at specific indices."""
        if isinstance(self.values, np.ndarray):
            self.values = self.values[indices]
        else:
            self.values = [self.values[i] for i in indices]

    def __repr__(self):
        sample = self.values[:5]
        return f"ObservableData(name={self.name}, values={sample}{'...' if len(self.values) > 5 else ''})"
    
    def get_bootstrap_samples(self, n_bootstrap: int, seed: int) -> List[ObservableData]:
        """Generate independent (IID) bootstrap samples."""
        rng = np.random.default_rng(seed)
        arr = np.array(self.values)
        if len(arr) == 0:
            return [ObservableData(self.name, []) for _ in range(n_bootstrap)]
            
        indices = rng.integers(0, len(arr), size=(n_bootstrap, len(arr)))
        return [ObservableData(self.name, arr[idx]) for idx in indices]


class FileData:
    def __init__(self, path: str):
        self.name = Path(path).stem
        self.path = Path(path)
        self.observables: List[ObservableData] = []

    def read_file(self):
        """Read CSV and populate observables."""
        if not self.path.exists():
            raise FileNotFoundError(f"{self.path} not found")

        data = np.genfromtxt(self.path, delimiter=',', names=True, encoding='utf-8-sig', ndmin=1)
        if data.dtype.names is None:
            return self

        self.observables = []
        for name in data.dtype.names:
            self.observables.append(ObservableData(name, data[name]))

        return self

    def align_lengths(self) -> int:
        """
        Truncate all observable arrays to the same minimum length.
        This guards against malformed rows or partial writes.
        """
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

        step_obs = next((o for o in self.observables if o.name in ["# step", "step"]), None)
        if step_obs is None: return

        keep_indices = [i for i, s in enumerate(step_obs.values) if s >= min_step]
        if len(keep_indices) == len(step_obs.values): return

        for obs in self.observables:
            obs.slice(keep_indices)

    def get(self, name: str) -> ObservableData:
        for obs in self.observables:
            if obs.name == name:
                return obs
        raise ValueError(f"Observable '{name}' not found in {self.path}")
    
    def get_bootstrap(self, n_bootstrap: int, seed: int) -> List[FileData]:
        """Legacy IID bootstrap (row-based)."""
        bootstrap_files: List[FileData] = []
        for i in range(n_bootstrap):
            bootstrap_file = FileData(self.path)
            bootstrap_file.observables = [
                obs.get_bootstrap_samples(1, seed + i)[0] for obs in self.observables
            ]
            bootstrap_files.append(bootstrap_file)
        return bootstrap_files

    def get_blocked_bootstrap(self, n_bootstrap: int, block_size: int, seed: int) -> List[FileData]:
        """
        Generate bootstrap samples using BLOCKED resampling of STEPS.
        1. Identifies unique steps.
        2. Resamples steps in blocks of size `block_size` to preserve autocorrelation.
        3. Selects all data rows corresponding to the sampled steps.
        """
        step_obs = next((o for o in self.observables if o.name in ["# step", "step"]), None)
        if step_obs is None or block_size <= 1:
            # Fallback to standard if no step info or block_size is 1
            return self.get_bootstrap(n_bootstrap, seed)
            
        steps = np.array(step_obs.values)
        unique_steps = np.unique(steps) # Sorted unique steps
        n_steps = len(unique_steps)
        
        # Map step -> list of row indices
        # (Optimization: We assume data might be interleaved, so we build a full map)
        step_to_indices = defaultdict(list)
        for idx, s in enumerate(steps):
            step_to_indices[s].append(idx)
            
        rng = np.random.default_rng(seed)
        bootstrap_files = []
        
        for _ in range(n_bootstrap):
            # Block Bootstrap the unique steps
            n_blocks = int(np.ceil(n_steps / block_size))
            max_start = max(1, n_steps - block_size + 1)
            
            # Pick random start positions for blocks
            starts = rng.integers(0, max_start, size=n_blocks)
            
            selected_row_indices = []
            
            for start_idx in starts:
                end_idx = min(start_idx + block_size, n_steps)
                block_steps = unique_steps[start_idx : end_idx]
                
                for s in block_steps:
                    selected_row_indices.extend(step_to_indices[s])
            
            # Construct new FileData with these rows
            # We convert to numpy for faster indexing
            new_fd = FileData(str(self.path))
            new_fd.observables = []
            
            for obs in self.observables:
                arr = np.array(obs.values)
                # Take slice
                new_vals = arr[selected_row_indices].tolist()
                new_fd.observables.append(ObservableData(obs.name, new_vals))
                
            bootstrap_files.append(new_fd)
            
        return bootstrap_files


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


def combine_file_data(files: List[FileData], min_step: int = 0, source_name: str = "combined") -> FileData | None:
    """
    Combine multiple FileData objects by:
    1. Optional thermalization cut per file.
    2. Taking only columns common to all files.
    3. Appending rows file-by-file.
    """
    if not files:
        return None

    prepared: List[FileData] = []
    for fd in files:
        fd.align_lengths()
        if min_step > 0:
            fd.remove_thermalization(min_step)
        fd.align_lengths()
        if not fd.observables:
            continue
        if min(len(obs.values) for obs in fd.observables) <= 0:
            continue
        prepared.append(fd)

    if not prepared:
        return None

    common_names = set(obs.name for obs in prepared[0].observables)
    for fd in prepared[1:]:
        common_names &= {obs.name for obs in fd.observables}

    if not common_names:
        return None

    ordered_names = [obs.name for obs in prepared[0].observables if obs.name in common_names]

    combined = FileData(f"{source_name}.out")
    combined.observables = [ObservableData(name, []) for name in ordered_names]
    combined_map = {obs.name: obs for obs in combined.observables}

    for fd in prepared:
        local_map = {obs.name: obs for obs in fd.observables}
        chunk_len = min(len(local_map[name].values) for name in ordered_names)
        if chunk_len <= 0:
            continue
        for name in ordered_names:
            combined_map[name].extend(local_map[name].values[:chunk_len])

    if not combined.observables:
        return None
    if min((len(obs.values) for obs in combined.observables), default=0) <= 0:
        return None

    combined.align_lengths()
    return combined
