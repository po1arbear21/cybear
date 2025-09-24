"""
Data loader module - handles FBS to MAT conversion and data loading
"""

import subprocess
import numpy as np
from pathlib import Path
from scipy.io import loadmat

class DataLoader:
    """Handles all data loading operations"""

    def __init__(self, flott_cli='flott-cli'):
        self.flott_cli = flott_cli

    def convert_fbs_to_mat(self, job_name, fbs_file):
        """Convert .fbs to .mat using flott-cli"""
        # Find the project root (where jobs/ directory is)
        current = Path.cwd()
        if (current / "jobs").exists():
            root = current
        elif (current.parent / "jobs").exists():
            root = current.parent
        else:
            root = current  # Fallback

        fbs_path = root / f"jobs/{job_name}/run/{fbs_file}"
        mat_path = fbs_path.with_suffix('.mat')

        # Skip if already converted and up to date
        if mat_path.exists() and mat_path.stat().st_mtime > fbs_path.stat().st_mtime:
            return mat_path

        # Run flott-cli conversion (suppress warnings)
        cmd = [self.flott_cli, '-f', str(fbs_path), 'export', '-a', '-r', '-o', str(mat_path)]
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check if the output file was created (more reliable than checking return code)
        if not mat_path.exists():
            # Only show actual errors, not warnings
            if result.returncode != 0 and result.stderr:
                # Filter out common warnings
                error_lines = []
                for line in result.stderr.split('\n'):
                    if line and not any(w in line for w in ['WARNING:', 'ERROR: Labels:', 'QSocketNotifier:']):
                        error_lines.append(line)
                if error_lines:
                    raise RuntimeError(f"flott-cli conversion failed: {' '.join(error_lines)}")
            raise RuntimeError(f"flott-cli failed to create output file")

        print(f"✓ Converted: {fbs_path.name} → {mat_path.name}")
        return mat_path

    def load_data(self, job_name, fbs_file):
        """Load data from FBS file (converts to MAT first)"""
        mat_path = self.convert_fbs_to_mat(job_name, fbs_file)

        # Load MAT file
        data = loadmat(str(mat_path))

        # Clean up MATLAB metadata
        cleaned = {}
        for key, value in data.items():
            if not key.startswith('__'):
                cleaned[key] = value

        return cleaned

    def detect_fields(self, data, mode='smart'):
        """Detect plottable fields from data

        mode: 'smart' - only important physics fields
              'all' - all plottable fields
        """

        if mode == 'all':
            # Return all plottable fields
            fields = []
            for key, value in data.items():
                # Skip coordinate/grid arrays
                skip_patterns = ['_vertices', '_length', '_surface', '_volume',
                               'grid_', 'dummy_', 'edge', 'cell', 'face']
                if any(p in key for p in skip_patterns):
                    continue

                # Check if plottable (not scalar)
                if hasattr(value, 'shape') and value.size > 1:
                    fields.append(key)
            return sorted(fields)

        # Smart mode - priority physics fields
        priority_fields = [
            'pot',      # Potential
            'ndens',    # Electron density
            'pdens',    # Hole density
            'Ex', 'Ey', 'Ez',  # E-field components
            'efield',   # E-field magnitude
            'current', 'jn', 'jp',  # Current densities
            'mobility', # Carrier mobility
            'recomb',   # Recombination
            'temperature', 'T'  # Temperature
        ]

        fields = []
        for field in priority_fields:
            if field in data:
                val = data[field]
                # Check if actually plottable
                if hasattr(val, 'shape') and val.size > 1:
                    fields.append(field)

        return fields

    def get_coordinates(self, data):
        """Extract coordinate arrays from data"""

        # Try standard names first
        x = data.get('x_vertices', data.get('vertices', None))
        y = data.get('y_vertices', None)

        # Fallback to searching for coordinate arrays
        if x is None:
            for key in data.keys():
                if 'x' in key.lower() and 'vertices' in key.lower():
                    x = data[key]
                    break

        if y is None:
            for key in data.keys():
                if 'y' in key.lower() and 'vertices' in key.lower():
                    y = data[key]
                    break

        # Convert to nm and flatten
        if x is not None:
            x = x.flatten() * 1e7  # m to nm
        if y is not None:
            y = y.flatten() * 1e7  # m to nm

        return x, y

    def get_dimensionality(self, data, field):
        """Determine if field is 1D, 2D, or 3D data"""
        field_data = data[field]
        shape = field_data.shape

        # Check for 1D
        if len(shape) == 1:
            return '1D'
        elif len(shape) == 2 and min(shape) == 1:
            return '1D'
        # Check for 2D
        elif len(shape) == 2 and min(shape) > 1:
            return '2D'
        elif len(shape) == 3:
            # Third dimension might be time/bias steps
            if shape[2] == 1:
                return '2D'
            else:
                return '2D-multi'  # Multiple frames

        return 'unknown'

    def detect_iv_curves(self, data):
        """Detect IV curve data (current-voltage pairs)

        Returns list of tuples: [(voltage_field, current_field, label), ...]
        """
        iv_curves = []

        # Find all voltage fields
        voltage_fields = []
        current_fields = []

        for key, val in data.items():
            if not hasattr(val, 'shape'):
                continue

            # Look for voltage data
            if any(v in key.upper() for v in ['V_', 'VOLT', 'BIAS', 'VG', 'VD', 'VS']):
                if val.size > 1:  # Must have multiple points
                    voltage_fields.append(key)

            # Look for current data
            elif any(i in key.upper() for i in ['I_', 'CURR', 'ID', 'IG', 'IS']):
                if val.size > 1:
                    current_fields.append(key)

        # Try to match voltage-current pairs
        matched = set()

        # First try exact name matching (V_GATE with I_GATE, etc.)
        for v_field in voltage_fields:
            base_name = v_field.replace('V_', '').replace('v_', '')
            for c_field in current_fields:
                if c_field not in matched:
                    c_base = c_field.replace('I_', '').replace('i_', '')
                    if base_name.lower() == c_base.lower():
                        # Found a match!
                        label = f"I-V ({base_name})"
                        iv_curves.append((v_field, c_field, label))
                        matched.add(c_field)
                        break

        # For unmatched currents, pair with any voltage if only one exists
        unmatched_currents = [c for c in current_fields if c not in matched]
        if len(unmatched_currents) > 0 and len(voltage_fields) == 1:
            for c_field in unmatched_currents:
                label = f"I-V ({c_field})"
                iv_curves.append((voltage_fields[0], c_field, label))

        return iv_curves