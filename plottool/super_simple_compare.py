#!/usr/bin/env python3
"""
Super Simple GALENE vs Cybear Comparison - CLEAN VERSION

Uses CybearPlotter styling but keeps the comparison logic dead simple.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

from cybear_plotter import CybearPlotter
from data_loader import CybearDataLoader

# ============================================================================
# SIMPLE CONFIGURATION
# ============================================================================

CYBEAR_JOB = 'nw_2D'           
VARIABLES = ['potential', 'n_density']  
STYLE = 'ieee'                 

# ============================================================================
# AUTO-MAGIC
# ============================================================================

def find_cybear_root():
    """Auto-find cybear project root"""
    current = Path.cwd()
    for path in [current] + list(current.parents):
        if (path / 'jobs').exists() and (path / 'plottool').exists():
            return path
    return current

def load_galene_data(galene_db_path):
    """Load and reshape GALENE data"""
    print(f"Loading GALENE data from: {galene_db_path}")
    
    data = loadmat(str(galene_db_path))
    galene_data = {}
    
    # Extract all data
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
    
    # Reshape 2D fields with Fortran order (MATLAB-style)
    if 'x_coord' in galene_data and 'y_coord' in galene_data:
        Nx, Ny = len(galene_data['x_coord']), len(galene_data['y_coord'])
        expected_size = Nx * Ny
        
        for field_name in ['potential', 'n_density']:
            if field_name in galene_data:
                field_data = galene_data[field_name]
                if hasattr(field_data, 'size') and field_data.size == expected_size:
                    galene_data[field_name] = field_data.reshape((Nx, Ny), order='F')
                    print(f"✓ Reshaped {field_name} to {galene_data[field_name].shape}")
    
    return galene_data

def create_comparison_plot(galene_data, cybear_data, var_name, output_dir):
    """Create simple side-by-side comparison with CybearPlotter styling"""
    
    # Variable mapping
    var_map = {
        'potential': {'galene': 'potential', 'cybear': 'pot'},
        'n_density': {'galene': 'n_density', 'cybear': 'ndens'}
    }
    
    if var_name not in var_map:
        return None
        
    galene_var = var_map[var_name]['galene']
    cybear_var = var_map[var_name]['cybear']
    
    if galene_var not in galene_data or cybear_var not in cybear_data:
        return None
    
    # Initialize CybearPlotter for proper styling
    plotter = CybearPlotter(style=STYLE)
    
    # Create figure with CybearPlotter styling
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    plotter._apply_style()
    
    # GALENE plot using CybearPlotter's field plotting logic
    galene_field = galene_data[galene_var]
    x_g = galene_data.get('x_coord', np.arange(galene_field.shape[0]))
    y_g = galene_data.get('y_coord', np.arange(galene_field.shape[1]))
    
    # Use proper field visualization from CybearPlotter
    if 'density' in galene_var.lower():
        im1 = ax1.contourf(x_g, y_g, galene_field.T, levels=50, cmap='plasma', 
                          norm=plt.matplotlib.colors.LogNorm())
    else:
        im1 = ax1.contourf(x_g, y_g, galene_field.T, levels=50, cmap='viridis')
    
    cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.8)
    ax1.set_title(f'GALENE: {plotter.labels.get_label(galene_var)[0]}')
    ax1.set_xlabel(plotter.labels.get_label('x')[0])
    ax1.set_ylabel(plotter.labels.get_label('y')[0])
    
    # Cybear plot using CybearPlotter's field plotting logic
    cybear_field = cybear_data[cybear_var]
    x_c = cybear_data.get('x_vertices', np.arange(cybear_field.shape[1]))
    y_c = cybear_data.get('y_vertices', np.arange(cybear_field.shape[0]))
    
    # Use proper field visualization from CybearPlotter
    if 'dens' in cybear_var.lower():
        im2 = ax2.contourf(x_c, y_c, cybear_field.T, levels=50, cmap='plasma',
                          norm=plt.matplotlib.colors.LogNorm())
    else:
        im2 = ax2.contourf(x_c, y_c, cybear_field.T, levels=50, cmap='viridis')
    
    cbar2 = plt.colorbar(im2, ax=ax2, shrink=0.8)
    ax2.set_title(f'Cybear: {plotter.labels.get_label(cybear_var)[0]}')
    ax2.set_xlabel(plotter.labels.get_label('x')[0])
    ax2.set_ylabel(plotter.labels.get_label('y')[0])
    
    plt.tight_layout()
    
    # Save using CybearPlotter's save method
    output_file = output_dir / f"proper_{var_name}.pdf"
    plotter._save_figure(fig, str(output_file))
    plt.close(fig)
    
    return str(output_file)

def main():
    print("🎯 Clean GALENE vs Cybear Comparison")
    print("=" * 50)
    
    # Auto-locate files
    cybear_root = find_cybear_root()
    galene_db_path = cybear_root / 'galene_database.mat'
    output_dir = cybear_root / 'clean_comparison'
    output_dir.mkdir(exist_ok=True)
    
    try:
        # Load data
        galene_data = load_galene_data(galene_db_path)
        
        data_loader = CybearDataLoader()
        cybear_data = data_loader.load_simulation_data(CYBEAR_JOB, 'SS.fbs')
        
        print(f"✓ GALENE: {len([k for k in galene_data.keys() if not k.startswith('_')])} variables")
        print(f"✓ Cybear: {len(cybear_data)} variables")
        
        # Generate plots
        generated_files = []
        for var in VARIABLES:
            output_file = create_comparison_plot(galene_data, cybear_data, var, output_dir)
            if output_file:
                generated_files.append(output_file)
                print(f"✓ {var}: {Path(output_file).name}")
        
        print(f"\n🎉 Generated {len(generated_files)} clean comparison plots in:")
        print(f"   {output_dir}")
            
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    main()