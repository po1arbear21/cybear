#!/usr/bin/env python3
"""
Quick plotting script for common Cybear simulation results
Usage: python3 quick_plot.py <job_name>
"""

import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')

# Add plottool to path
sys.path.insert(0, str(Path(__file__).parent))

from cybear_plotter import CybearPlotter

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 quick_plot.py <job_name>")
        print("Example: python3 quick_plot.py vfet_test")
        sys.exit(1)

    job_name = sys.argv[1]

    # Initialize plotter with Nature style (default)
    plotter = CybearPlotter(style=['science', 'ieee'])

    print(f"Quick plotting for {job_name}/SS.fbs")
    print("=" * 50)

    # Standard plots for semiconductor device analysis
    plots_to_generate = [
        ('pot', '2d', 'Electric potential field'),
        ('ndens', '2d', 'Electron density field'),
        ('pot', 'x_linecut', 'Potential along channel'),
        ('ndens', 'x_linecut', 'Density along channel'),
    ]

    generated_files = []

    for field, plot_type, description in plots_to_generate:
        try:
            if plot_type == '2d':
                output_path = f"{job_name}_{field}_2d"
                fig = plotter.plot_2d_field(job_name, 'SS.fbs', field, save_path=output_path)
            elif plot_type.endswith('_linecut'):
                direction = plot_type.split('_')[0]
                output_path = f"{job_name}_{field}_{direction}_linecut"
                fig = plotter.plot_line_cut(job_name, 'SS.fbs', field,
                                          direction=direction, save_path=output_path)

            generated_files.append(f"{output_path}.pdf")
            print(f"✓ {description}: {output_path}.pdf")

        except Exception as e:
            print(f"✗ Failed {description}: {e}")

    # Try convergence plot
    try:
        output_path = f"{job_name}_convergence"
        fig = plotter.plot_convergence(job_name, save_path=output_path)
        generated_files.append(f"{output_path}.pdf")
        print(f"✓ Convergence analysis: {output_path}.pdf")
    except Exception as e:
        print(f"✗ Convergence plot not available: {e}")

    print("=" * 50)
    print(f"Generated {len(generated_files)} publication-ready figures:")
    for file in generated_files:
        print(f"  {file}")

if __name__ == "__main__":
    main()
