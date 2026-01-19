# Gemini Code Assistant Context

This document provides context for the Gemini code assistant to help it understand the project and provide better assistance.

## Project Overview

This project is a Lattice QCD analysis suite focusing on extracting physical observables like the static potential $V(R)$, string tension $\sigma$, and the Sommer scale $r_0$ from Wilson loop data. It includes tools for searching across large datasets of simulation runs, caching analysis results, and performing robust statistical error analysis (bootstrap).

## Key Technologies

*   **Language:** Python 3.10+
*   **Core Libraries:**
    *   `numpy`: Heavy numerical lifting (fitting, array manipulation).
    *   `scipy`: Curve fitting (`curve_fit`) and interpolation.
    *   `matplotlib`: Plotting (in `plot_*.py` scripts).
    *   `PyYAML`: Parsing simulation parameters (`input.yaml`).
    *   `concurrent.futures`: Parallel processing in `search_data.py`.

## Project Structure

### Core Analysis Logic
*   **`calculator.py`**: The physics engine.
    *   Implements the `Calculator` class which lazily computes observables.
    *   Uses a decorator-based registry (`@register`) to define calculation methods (e.g., `_calc_r0`, `_calc_chi`).
    *   **Features**: Bootstrap error propagation, global weighted fitting for `r0_chi`, robust handling of negative/noisy data.
*   **`run_evaluation.py`**: The evaluation pipeline.
    *   Serves as the bridge between raw CSV data and the `Calculator`.
    *   Manages **caching** of results in `../data/calcResult/` (hashed by run path).
    *   Main entry point: `get_or_calculate(run_path)`.
*   **`data_organizer.py`**: Data structures.
    *   Classes: `FileData`, `ObservableData`, `VariableData`.
    *   Handles reading raw `.out` files from simulation runs.

### CLI Tools & Workflows
*   **`search_data.py`**: The primary tool for aggregating results.
    *   **Usage**: `python search_data.py <root_dir> [filters] [outputs]`
    *   Scans directories, matches parameters (e.g., `beta=2.1`), runs analysis in parallel, and prints a table.
*   **`analysis_notebook.ipynb`**: Interactive debugging.
    *   Use this to step through the analysis of a *single* problematic run.
    *   Great for visualizing why a specific fit (like `r0_chi`) might be failing.

### Plotting Scripts
*   `plot_potential.py`, `plot_history.py`, `plot_search.py`: Generate visualizations from the analyzed data.

## How to Run

### 1. Bulk Analysis & Search
To find all runs with $\beta=2.1$ and display their $r_0$ and $\chi$ values:

```bash
# Syntax: python search_data.py <data_root> <filters> <output_columns>
python search_data.py ../data "beta=2.1" "L0=43" r0 r0_err r0_chi
```

*   **Filters**: `KEY=VALUE` (e.g., `beta=2.1`).
*   **Outputs**: Column names (e.g., `r0`, `chi`, `a`, `tau_int`).

### 2. Debugging a Specific Run
1.  Open `analysis_notebook.ipynb` in VS Code.
2.  Set the `run_path` variable to the directory of the run you want to investigate.
3.  Run the cells to see step-by-step extraction of Wilson loops, potentials, and forces.

### 3. Plotting
(Example inferred)
```bash
python plot_potential.py --run-path ../data/20251221/25
```

## Important Logic Details
   

*   **r0_chi Calculation**:
    *   The code attempts to solve $r^2 \chi(r, t_{large}) = 1.65$.
    *   If data is noisy or negative (non-confining signal), it falls back to a **global weighted linear fit** of $r^2 F(r)$ vs $r^2$ to estimate the crossing point.
    *   This robustness is critical for analyzing runs with low statistics or large finite-volume effects.
*   **Caching**:
    *   Results are cached to JSON. If you change the analysis logic (e.g., in `calculator.py`), you must bump `CALC_VERSION` in `run_evaluation.py` to invalidate old caches.

## Coding Conventions

*   **Type Hinting**: Use Python's `typing` module (e.g., `List`, `Dict`, `Optional`) for function signatures.
*   **Error Handling**: Analysis functions should return `None` or `NaN` rather than crashing, allowing `search_data.py` to continue processing other runs.
*   **Modular Design**: Keep physics logic in `calculator.py` and data management in `data_organizer.py`.
