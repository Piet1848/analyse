import os
import sys
import pprint

# Add the current directory to the path to find our custom modules
if '.' not in sys.path:
    sys.path.insert(0, '.')

from load_input_yaml import load_params
import run_evaluation

print("Setup complete. Modules imported.")


run_path = "../data/20260205/14"  # Example path with N/A and nan values
calculated_data = run_evaluation.get_or_calculate(run_path, force_recalc=True)

if 'error' in calculated_data:
    print("Error:", calculated_data['error'])
else:
    print("--- Result of get_or_calculate ---")
    pprint.pprint(calculated_data)