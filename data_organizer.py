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

    def __repr__(self):
        sample = self.values[:5]
        return f"ObservableData(name={self.name}, values={sample}{'...' if len(self.values) > 5 else ''})"
    
    def get_bootstrap_samples(self, n_bootstrap: int, seed: int) -> List[ObservableData]:
        """Generate bootstrap samples for this observable."""
        rng = np.random.default_rng(seed)
        arr = np.array(self.values)
        indices = rng.integers(0, len(arr), size=(n_bootstrap, len(arr)))
        return [ObservableData(self.name, arr[idx].tolist()) for idx in indices]


class FileData:
    def __init__(self, path: str):
        self.name = Path(path).stem
        self.path = Path(path)
        self.observables: List[ObservableData] = []

    def read_file(self):
        """Read CSV and populate observables."""
        if not self.path.exists():
            raise FileNotFoundError(f"{self.path} not found")

        with self.path.open() as f:
            reader = csv.reader(f)
            headers = next(reader)
            self.observables = [ObservableData(name) for name in headers]

            for row in reader:
                for obs, val in zip(self.observables, row):
                    obs.append(float(val))

        return self

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
        """Scan the folder for .out files and return a dict of FileData lists."""
        data_dict: Dict[str, List[FileData]] = {}
        for file_path in self.path.iterdir():
            if file_path.is_file() and file_path.suffix == ".out":
                fd = FileData(str(file_path))
                fd.read_file()
                data_dict.setdefault(file_path.stem, []).append(fd)
        return data_dict


    def get_bootstrap_data(self, n_bootstrap: int, seed: int) -> List[ExperimentData]:
        """Generate bootstrap samples for all files and observables in the experiment."""
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
    bootstrap_samples = exp_data.get_bootstrap_data(n_bootstrap=1000, seed=42)
    for obs_name, samples in bootstrap_samples.items():
        print(f"Observable: {obs_name}, Number of bootstrap samples: {len(samples)}")