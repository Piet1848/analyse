#!/usr/bin/env python3
from typing import Any, FrozenSet, Callable
import data_organizer
import numpy as np

Key = tuple[str, FrozenSet[tuple[str, Any]]]

def make_key(name: str, params: dict[str, Any]) -> Key:
    return (name, frozenset(params.items()))

def get_key_name(key: Key) -> str:
    return key[0]

def get_key_params(key: Key) -> dict[str, Any]:
    return dict(key[1])


class Calculator:
    _registry: dict[str, Callable] = {}

    def __init__(self, file_data: data_organizer.FileData):
        self.file_data = file_data
        self.variables: dict[Key, data_organizer.VariableData] = {}

    @classmethod
    def register(cls, name: str):
        """Decorator to register a method as a calculator for a specific variable name."""
        def decorator(method):
            cls._registry[name] = method
            return method
        return decorator

    def get_variable(self, name: str, **params) -> data_organizer.VariableData:
        key = make_key(name, params)
        if key in self.variables:
            return self.variables[key]

        if name in self._registry:
            try:
                # Call the registered method (passing self and params)
                var_data = self._registry[name](self, **params)
            except TypeError as e:
                raise TypeError(f"Error calculating '{name}'. Check if params {params} match method signature.") from e
            
            self.variables[key] = var_data
            return var_data

        raise KeyError(f"No calculator registered for variable '{name}'")
    
    def get_observable(self, obs_name: str) -> data_organizer.ObservableData:
        for obs in self.file_data.observables:
            if obs.name == obs_name:
                return obs
        raise KeyError(f"Observable '{obs_name}' not found in file data.")
    
    ### Variable implementations ###

    @register("W_R_T")
    def _calc_W_R_T(self, R: int, T: int) -> data_organizer.VariableData:
        try:
            obs_val = self.file_data.get("W_temp").values 
            obs_L   = self.file_data.get("L").values
            obs_T   = self.file_data.get("T").values
        except ValueError as e:
            raise KeyError(f"Missing required columns (L, T, or W_temp) in file: {e}")

        arr_val = np.array(obs_val)
        arr_L   = np.array(obs_L)
        arr_T   = np.array(obs_T)

        # 3. Create a Mask (Where L==R AND T==T)
        mask = (arr_L == R) & (arr_T == T)
        selected_values = arr_val[mask]

        if len(selected_values) == 0:
            raise ValueError(f"No data found for R={R}, T={T}")

        var_data = data_organizer.VariableData("W_R_T")
        
        mean_val = np.mean(selected_values)
        var_data.set_value(mean_val, bootstrap_samples=selected_values, R=R, T=T)

        return var_data