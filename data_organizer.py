from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List
import csv
import numpy as np


class ObservableData:
    def __init__(self, name: str, values: List[float] = None):
        self.name = name
        self.values = values or []

    def append(self, value: float):
        self.values.append(value)

    def extend(self, values: List[float]):
        self.values.extend(values)
        
    def slice(self, indices: List[int]):
        """Keep only values at specific indices."""
        # Use numpy for faster indexing if possible, but list comp is fine here
        self.values = [self.values[i] for i in indices]

    def __repr__(self):
        sample = self.values[:5]
        return f"ObservableData(name={self.name}, values={sample}{'...' if len(self.values) > 5 else ''})"
    
    def get_bootstrap_samples(self, n_bootstrap: int, seed: int) -> List[ObservableData]:
        """Generate bootstrap samples for this observable."""
        rng = np.random.default_rng(seed)
        arr = np.array(self.values)
        if len(arr) == 0:
            return [ObservableData(self.name, []) for _ in range(n_bootstrap)]
            
        indices = rng.integers(0, len(arr), size=(n_bootstrap, len(arr)))
        return [ObservableData(self.name, arr[idx].tolist()) for idx in indices]


class FileData:
    def __init__(self, path: str):
        self.name = Path(path).stem
        self.path = Path(path)
        self.observables: List[ObservableData] = []

    def read_file(self):
        """Read CSV and populate observables. Skip empty files gracefully."""
        if not self.path.exists():
            raise FileNotFoundError(f"{self.path} not found")

        with self.path.open() as f:
            reader = csv.reader(f)
            headers = next(reader, None)
            if headers is None:
                return self

            headers = [h.strip().lstrip("\ufeff") for h in headers]
            self.observables = [ObservableData(name) for name in headers]

            for row_idx, row in enumerate(reader, start=1):
                if not row or all(not cell.strip() for cell in row):
                    continue
                if len(row) < len(self.observables):
                    continue # Skip malformed rows
                
                for obs, val in zip(self.observables, row):
                    try:
                        obs.append(float(val))
                    except ValueError:
                        pass # Skip non-numeric

        return self

    def remove_thermalization(self, min_step: int):
        """
        Discard all data points where 'step' < min_step.
        Looks for column named '# step' or 'step'.
        """
        step_obs = next((o for o in self.observables if o.name in ["# step", "step"]), None)
        
        if step_obs is None:
            # If no step column found, we cannot reliably thermalize by step count.
            # You might optionally warn here.
            return

        # Find indices where step >= min_step
        keep_indices = [i for i, s in enumerate(step_obs.values) if s >= min_step]
        
        # If no data needs to be removed
        if len(keep_indices) == len(step_obs.values):
            return

        # Update ALL observables in this file to keep only these rows
        for obs in self.observables:
            obs.slice(keep_indices)

    def get(self, name: str) -> ObservableData:
        for obs in self.observables:
            if obs.name == name:
                return obs
        raise ValueError(f"Observable '{name}' not found in {self.path}")
    
    def get_bootstrap(self, n_bootstrap: int, seed: int) -> List[FileData]:
        """Generate bootstrap samples for all observables in the file."""
        bootstrap_files: List[FileData] = []
        for i in range(n_bootstrap):
            bootstrap_file = FileData(self.path)
            bootstrap_file.observables = [
                obs.get_bootstrap_samples(1, seed + i)[0] for obs in self.observables
            ]
            bootstrap_files.append(bootstrap_file)
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

    def get_bootstrap_data(self, n_bootstrap: int, seed: int) -> List[ExperimentData]:
        bootstrap_experiments: List[ExperimentData] = []
        for i in range(n_bootstrap):
            bootstrap_exp = ExperimentData(str(self.path))
            bootstrap_exp.data = {}
            for obs_name, file_list in self.data.items():
                bootstrap_exp.data[obs_name] = []
                for file_data in file_list:
                    bootstrap_file = FileData(str(file_data.path))
                    bootstrap_file.observables = [
                        obs.get_bootstrap_samples(1, seed + i)[0] for obs in file_data.observables
                    ]
                    bootstrap_exp.data[obs_name].append(bootstrap_file)
            bootstrap_experiments.append(bootstrap_exp)
        return bootstrap_experiments

if __name__ == "__main__":
    exp_data = ExperimentData("../data/20251124/01/")
    # Test thermalization
    # for files in exp_data.data.values():
    #     for f in files: f.remove_thermalization(100)
    pass