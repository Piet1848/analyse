# Gemini Code Assistant Context

This document provides context for the Gemini code assistant to help it understand the project and provide better assistance.

## Project Overview

This project consists of a collection of Python scripts for data analysis and visualization. The scripts are designed to process, analyze, and plot data as part of a research project, likely for a Master's thesis. The main tasks include searching data, organizing it, performing calculations, and generating various plots.

## Key Technologies

*   **Language:** Python 3
*   **Libraries (likely):**
    *   `numpy` and `pandas` for data manipulation and numerical operations.
    *   `matplotlib`, `seaborn`, or `plotly` for data visualization.
    *   `PyYAML` for reading YAML configuration files (`load_input_yaml.py`).

## Project Structure

The project is organized into several Python scripts, each with a specific purpose:

*   `analyze_wilson.py`: Script for a specific analysis, possibly related to a "Wilson" dataset or method.
*   `calculator.py`: Performs calculations on the data.
*   `data_organizer.py`: Handles organization and preprocessing of data.
*   `load_input_yaml.py`: Utility for loading configuration or input data from YAML files.
*   `plot_history.py`: Generates history plots.
*   `plot_potential.py`: Generates potential plots.
*   `plot_search.py`: Generates plots related to search results.
*   `run_evaluation.py`: Main script to run an evaluation pipeline.
*   `search_data.py`: Script to search through the dataset.
*   `__pycache__/`: Directory for Python's cached bytecode.
*   `.gitignore`: Specifies files to be ignored by Git.

## How to Run

The scripts are intended to be run from the command line.

```bash
# Example of running a script
python run_evaluation.py --input-file data.yaml
```

Please inspect individual scripts for specific command-line arguments.

## Coding Conventions

*   Please adhere to the **PEP 8 Style Guide for Python Code**.
*   Maintain the existing coding style and patterns found in the project's files.
*   When adding new features, ensure they are modular and well-documented.
*   Add comments to explain complex logic, not to describe what the code is doing.
