#!/usr/bin/env python3
"""
GALENE vs Cybear Comparison - Simple Configuration Style
Just modify the parameters below and run this file!
"""

import matplotlib
matplotlib.use('Agg')
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat

# Add plottool to path
sys.path.insert(0, str(Path(__file__).parent))

from cybear_plotter import CybearPlotter
from data_loader import CybearDataLoader

# ============================================================================
# MODIFY THESE PARAMETERS - THEN JUST RUN THIS FILE!
# ============================================================================

# Basic settings
CYBEAR_JOB = 'nw_2D'           # Cybear job name
CYBEAR_FBS = 'SS.fbs'         # Cybear FBS file
GALENE_DB = 'galene_database.mat'  # GALENE database file (relative to project root)
STYLE = 'ieee'                # Options: 'nature', 'ieee', 'retro', 'high-vis'

# Variables to compare (will generate separate plots for each)
VARIABLES_TO_COMPARE = ['potential', 'n_density']

# Output settings
OUTPUT_PREFIX = 'comparison'  # Prefix for output files
SAVE_TO_FOLDER = 'comparison_plots'  # Output folder

# Variable mapping between GALENE and Cybear naming
VARIABLE_MAPPING = {
    'potential': {'galene': 'potential', 'cybear': 'pot'},
    'n_density': {'galene': 'n_density', 'cybear': 'ndens'},
    'p_density': {'galene': 'p_density', 'cybear': 'pdens'},
}

# Field-specific plotting settings
FIELD_SETTINGS = {
    'potential': {
        'colormap': 'RdBu_r',
        'log_scale': False
    },
    'n_density': {
        'colormap': 'plasma', 
        'log_scale': True
    },
    'p_density': {
        'colormap': 'plasma',
        'log_scale': True
    }
}

# ============================================================================
# AUTO-MAGIC HELPER FUNCTIONS - DON'T MODIFY BELOW THIS LINE
# ============================================================================

def find_cybear_root():
    """Auto-find cybear project root"""
    current = Path.cwd()
    for path in [current] + list(current.parents):
        if (path / 'jobs').exists() and (path / 'plottool').exists():
            return path
    return current

def load_galene_data(galene_db_path):
    """Load and reshape GALENE data from .mat file"""
    print(f"Loading GALENE data from: {galene_db_path}")
    
    data = loadmat(str(galene_db_path))
    galene_data = {}
    
    # Extract all data from structured arrays
    for struct_name in ['sim', 'geo']:
        if struct_name in data:
            struct = data[struct_name][0, 0]
            for field_name in struct.dtype.names:
                if not field_name.endswith('_info'):
                    field_data = struct[field_name]
                    if hasattr(field_data, 'flatten'):
                        galene_data[field_name] = field_data.flatten()
                    else:
                        galene_data[field_name] = field_data
    
    # Reshape 2D fields using coordinate information
    if 'x_coord' in galene_data and 'y_coord' in galene_data:
        Nx, Ny = len(galene_data['x_coord']), len(galene_data['y_coord'])
        expected_size = Nx * Ny
        
        print(f"Grid dimensions: {Nx} × {Ny} = {expected_size}")
        
        # Fields that should be 2D
        field_variables = ['potential', 'n_density', 'p_density', 'electric_field']
        
        for field_name in field_variables:
            if field_name in galene_data:
                field_data = galene_data[field_name]
                if hasattr(field_data, 'size') and field_data.size == expected_size:
                    # Reshape with Fortran order (MATLAB-style)
                    galene_data[field_name] = field_data.reshape((Nx, Ny), order='F')
                    print(f"✓ Reshaped {field_name} to {galene_data[field_name].shape}")
    
    return galene_data

def load_cybear_data(job_name, fbs_file):
    """Load Cybear data using existing data loader"""
    print(f"Loading Cybear data from: {job_name}/{fbs_file}")
    
    data_loader = CybearDataLoader()
    cybear_data = data_loader.load_simulation_data(job_name, fbs_file)
    
    return cybear_data

def create_cybear_plot(job_name, fbs_file, var_name, plotter, output_dir):
    """Create Cybear plot using existing CybearPlotter methods - guaranteed to work!"""
    cybear_var = VARIABLE_MAPPING[var_name]['cybear']
    
    # Get plot settings
    settings = FIELD_SETTINGS.get(var_name.replace('_density', '').replace('_', ''), {})
    colormap = settings.get('colormap', 'viridis')
    
    # Use CybearPlotter's native 2D plotting method
    output_file = output_dir / f"{OUTPUT_PREFIX}_cybear_{var_name}"
    
    fig = plotter.plot_2d_field(job_name, fbs_file, cybear_var,
                              colormap=colormap,
                              save_path=str(output_file))
    
    return f"{output_file}.pdf"

def create_galene_plot(galene_data, var_name, plotter, output_dir):
    """Create GALENE plot manually but avoid LaTeX issues"""
    galene_var = VARIABLE_MAPPING[var_name]['galene']
    field_data = galene_data[galene_var]
    
    # Get plot settings
    settings = FIELD_SETTINGS.get(var_name.replace('_density', '').replace('_', ''), {})
    colormap = settings.get('colormap', 'viridis')
    log_scale = settings.get('log_scale', False)
    
    # Create plot with minimal matplotlib - avoid CybearPlotter's LaTeX styling
    plt.style.use('default')  # Reset to avoid LaTeX issues
    fig, ax = plt.subplots(figsize=(6, 5))
    
    # Get coordinates
    x_coord = galene_data.get('x_coord', np.arange(field_data.shape[0]))
    y_coord = galene_data.get('y_coord', np.arange(field_data.shape[1]))
    
    # Plot with appropriate scaling
    if log_scale and np.all(field_data > 0):
        im = ax.contourf(x_coord, y_coord, field_data.T, levels=50, cmap=colormap,
                        norm=plt.matplotlib.colors.LogNorm())
    else:
        im = ax.contourf(x_coord, y_coord, field_data.T, levels=50, cmap=colormap)
    
    # Add colorbar
    plt.colorbar(im, ax=ax, shrink=0.8)
    
    # Simple labels - no LaTeX
    ax.set_title(f'GALENE: {var_name}')
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    
    # Save manually
    output_file = output_dir / f"{OUTPUT_PREFIX}_galene_{var_name}.pdf"
    fig.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close(fig)
    
    return str(output_file)

# ============================================================================
# AUTO-EXECUTION - DON'T MODIFY BELOW THIS LINE
# ============================================================================

def main():
    print("🎯 GALENE vs Cybear Comparison")
    print(f"Cybear job: {CYBEAR_JOB}/{CYBEAR_FBS}")
    print(f"GALENE database: {GALENE_DB}")
    print(f"Style: {STYLE}")
    print(f"Variables: {VARIABLES_TO_COMPARE}")
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

    # Auto-locate files
    cybear_root = find_cybear_root()
    galene_db_path = cybear_root / GALENE_DB
    output_dir = Path(SAVE_TO_FOLDER)
    output_dir.mkdir(exist_ok=True)

    try:
        # Load data
        print("\nLoading simulation data...")
        galene_data = load_galene_data(galene_db_path)
        cybear_data = load_cybear_data(CYBEAR_JOB, CYBEAR_FBS)
        
        print(f"✓ GALENE: {len(galene_data)} variables")
        print(f"✓ Cybear: {len(cybear_data)} variables")

        # Generate individual plots for each variable and simulator
        print(f"\nGenerating comparison plots...")
        generated_files = []
        
        for var_name in VARIABLES_TO_COMPARE:
            if var_name not in VARIABLE_MAPPING:
                print(f"⚠ Variable '{var_name}' not in mapping, skipping...")
                continue
            
            galene_var = VARIABLE_MAPPING[var_name]['galene']
            cybear_var = VARIABLE_MAPPING[var_name]['cybear']
            
            # Check if variables exist
            if galene_var not in galene_data:
                print(f"⚠ GALENE variable '{galene_var}' not found, skipping...")
                continue
            
            if cybear_var not in cybear_data:
                print(f"⚠ Cybear variable '{cybear_var}' not found, skipping...")
                continue
            
            try:
                # Create GALENE plot (manual plotting to avoid LaTeX issues)
                galene_file = create_galene_plot(galene_data, var_name, plotter, output_dir)
                generated_files.append(galene_file)
                print(f"✓ GALENE {var_name}: {Path(galene_file).name}")
                
                # Create Cybear plot (use CybearPlotter - works perfectly)
                cybear_file = create_cybear_plot(CYBEAR_JOB, CYBEAR_FBS, var_name, plotter, output_dir)
                generated_files.append(cybear_file)
                print(f"✓ Cybear {var_name}: {Path(cybear_file).name}")
                
            except Exception as e:
                print(f"✗ Failed {var_name}: {e}")
                import traceback
                traceback.print_exc()

        print("=" * 60)
        print(f"✅ Generated {len(generated_files)} comparison plots:")
        for file in generated_files:
            print(f"   {file}")
        print(f"\n🎉 All plots saved to: {output_dir}")
        
        # Generate simple summary
        summary_file = output_dir / "comparison_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("GALENE vs Cybear Comparison Summary\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Cybear Job: {CYBEAR_JOB}/{CYBEAR_FBS}\n")
            f.write(f"GALENE Database: {GALENE_DB}\n")
            f.write(f"Variables Compared: {', '.join(VARIABLES_TO_COMPARE)}\n")
            f.write(f"Generated Files: {len(generated_files)}\n\n")
            f.write("Generated Files:\n")
            for file in generated_files:
                f.write(f"  - {Path(file).name}\n")
        
        print(f"📝 Summary saved to: {summary_file}")
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()