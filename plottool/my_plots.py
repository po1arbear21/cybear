#!/usr/bin/env python3
"""
Personal plotting script - MATLAB style
Just modify the parameters below and run this file!
"""

import matplotlib
matplotlib.use('Agg')
import sys
from pathlib import Path

# Add plottool to path
sys.path.insert(0, str(Path(__file__).parent))

from cybear_plotter import CybearPlotter

# ============================================================================
# MODIFY THESE PARAMETERS - THEN JUST RUN THIS FILE!
# ============================================================================

# Basic settings
JOB_NAME = 'nw_2D'
FBS_FILE = 'SS.fbs'
STYLE = 'ieee'  # Options: 'nature', 'ieee', 'retro', 'high-vis'

# Fields to plot
FIELDS_TO_PLOT = ['pot', 'ndens']

# Plot types to generate
PLOT_TYPES = ['2d', 'linecut']  # Options: '2d', 'linecut', 'convergence'

# # Axis limits (set to None for auto-scaling)
XLIM = None    # X limits in nm, or None for auto
YLIM = None    # Y limits in nm, or None for auto

# Output settings
OUTPUT_PREFIX = 'my'   # Prefix for output files
SAVE_TO_FOLDER = '.'   # Output folder

# Field-specific settings (optional - customize per field)
FIELD_SETTINGS = {
    'ndens': {
        'xlim': XLIM,       # Spatial X limits (nm)
        'ylim': YLIM,       # Spatial Y limits (nm)
        'vmin': 1e16,       # Minimum density value (cm^-3)
        'vmax': 1e20,       # Maximum density value (cm^-3)
        'colormap': 'plasma'
    },
    'pot': {
        'xlim': XLIM,       # Spatial X limits (nm)
        'ylim': YLIM,       # Spatial Y limits (nm)
        'vmin': None,       # Auto colorbar limits
        'vmax': None,       # Auto colorbar limits
        'colormap': 'RdBu_r'
    }
}

# ============================================================================
# AUTO-EXECUTION - DON'T MODIFY BELOW THIS LINE
# ============================================================================

def main():
    print(f"Generating plots for {JOB_NAME}/{FBS_FILE}")
    print(f"Style: {STYLE}")
    print(f"Fields: {FIELDS_TO_PLOT}")
    print(f"Plot types: {PLOT_TYPES}")
    print("=" * 60)

    # Initialize plotter with chosen style
    if STYLE == 'nature':
        plotter = CybearPlotter(style=['science', 'nature'])
    elif STYLE == 'ieee':
        plotter = CybearPlotter(style=['science', 'ieee'])
    elif STYLE == 'retro':
        plotter = CybearPlotter(style=['science', 'retro'])
    elif STYLE == 'high-vis':
        plotter = CybearPlotter(style=['science', 'high-vis'])
    else:
        plotter = CybearPlotter(style=['science', 'nature'])

    output_dir = Path(SAVE_TO_FOLDER)
    output_dir.mkdir(exist_ok=True)

    generated_files = []

    # Generate 2D plots
    if '2d' in PLOT_TYPES:
        print("\nGenerating 2D field plots...")
        for field in FIELDS_TO_PLOT:
            try:
                # Get field-specific settings or use defaults
                settings = FIELD_SETTINGS.get(field, {})
                xlim = settings.get('xlim', XLIM)
                ylim = settings.get('ylim', YLIM)
                vmin = settings.get('vmin', None)  # Colorbar minimum
                vmax = settings.get('vmax', None)  # Colorbar maximum
                colormap = settings.get('colormap', 'RdBu_r')

                output_path = output_dir / f"{OUTPUT_PREFIX}_{JOB_NAME}_{field}_2d"

                fig = plotter.plot_2d_field(JOB_NAME, FBS_FILE, field,
                                          xlim=xlim,
                                          ylim=ylim,
                                          vmin=vmin,
                                          vmax=vmax,
                                          colormap=colormap,
                                          save_path=str(output_path))

                generated_files.append(f"{output_path}.pdf")
                print(f"  ✓ {field}: {output_path}.pdf")

            except Exception as e:
                print(f"  ✗ Failed {field}: {e}")

    # Generate line cuts
    if 'linecut' in PLOT_TYPES:
        print("\nGenerating line cut plots...")
        for field in FIELDS_TO_PLOT:
            for direction in ['x', 'y']:
                try:
                    output_path = output_dir / f"{OUTPUT_PREFIX}_{JOB_NAME}_{field}_{direction}_cut"

                    fig = plotter.plot_line_cut(JOB_NAME, FBS_FILE, field,
                                              direction=direction,
                                              save_path=str(output_path))

                    generated_files.append(f"{output_path}.pdf")
                    print(f"  ✓ {field} {direction.upper()}-cut: {output_path}.pdf")

                except Exception as e:
                    print(f"  ✗ Failed {field} {direction}-cut: {e}")

    # Generate convergence plot
    if 'convergence' in PLOT_TYPES:
        print("\nGenerating convergence plot...")
        try:
            output_path = output_dir / f"{OUTPUT_PREFIX}_{JOB_NAME}_convergence"
            fig = plotter.plot_convergence(JOB_NAME, save_path=str(output_path))
            generated_files.append(f"{output_path}.pdf")
            print(f"  ✓ Convergence: {output_path}.pdf")
        except Exception as e:
            print(f"  ✗ Failed convergence: {e}")

    print("=" * 60)
    print(f"✅ Generated {len(generated_files)} files:")
    for file in generated_files:
        print(f"   {file}")
    print("\n🎉 All done!")

if __name__ == "__main__":
    main()
