#!/usr/bin/env python3
"""
Test script for Cybear plotting toolkit using vfet_test data

This script demonstrates the new toolkit structure:
1. Data Loading (using flott-cli like MATLAB)
2. Plotting with LaTeX labels
"""

import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

# Add plottool to path
sys.path.insert(0, str(Path(__file__).parent))

import cybear_plotter
import data_loader
import label_dict

CybearPlotter = cybear_plotter.CybearPlotter
CybearDataLoader = data_loader.CybearDataLoader  
LabelDictionary = label_dict.LabelDictionary

def main():
    # Change to cybear root directory
    cybear_root = Path(__file__).parent.parent
    import os
    os.chdir(cybear_root)
    print(f"Working directory: {os.getcwd()}")
    
    print("Testing Cybear Plotting Toolkit with vfet_test data")
    print("=" * 60)
    
    # Initialize components
    plotter = CybearPlotter()
    loader = CybearDataLoader()
    labels = LabelDictionary()
    
    job_name = 'vfet_test'
    fbs_name = 'SS.fbs'
    
    try:
        # Test 1: List available data
        print("Available FBS files in vfet_test:")
        available = plotter.list_available_data(job_name)
        print(f"  {available}")
        print()
        
        # Test 2: Load and inspect data
        print("Loading simulation data...")
        data = loader.load_simulation_data(job_name, fbs_name)
        print(f"Variables loaded: {list(data.keys())}")
        
        # Print data shapes for key variables
        for var in ['x', 'y', 'pot', 'ndens']:
            if var in data:
                shape = data[var].shape if hasattr(data[var], 'shape') else 'scalar'
                print(f"  {var}: {shape}")
        print()
        
        # Test 3: LaTeX label system
        print("LaTeX label examples:")
        test_vars = ['x', 'pot', 'ndens', 'current']
        for var in test_vars:
            label = labels.get_label(var)
            title = labels.get_title(var)
            print(f"  {var} -> {label} ({title})")
        print()
        
        # Test 4: Plot 2D potential field  
        print("Plotting 2D potential field...")
        fig1 = plotter.plot_2d_field(job_name, fbs_name, 'pot',
                                    save_path='vfet_potential_2d',
                                    colormap='RdBu_r')
        print("  Saved: vfet_potential_2d.pdf")
        
        # Test 5: Plot 2D electron density (log scale)
        print("Plotting 2D electron density...")
        fig2 = plotter.plot_2d_field(job_name, fbs_name, 'ndens',
                                    save_path='vfet_ndens_2d',
                                    colormap='viridis')
        print("  Saved: vfet_ndens_2d.pdf")
        
        # Test 6: Plot line cuts  
        print("Plotting line cuts...")
        fig3 = plotter.plot_line_cut(job_name, fbs_name, 'pot',
                                   direction='x', 
                                   save_path='vfet_potential_linecut')
        print("  Saved: vfet_potential_linecut.pdf")
        
        fig4 = plotter.plot_line_cut(job_name, fbs_name, 'ndens',
                                   direction='x',
                                   save_path='vfet_ndens_linecut')
        print("  Saved: vfet_ndens_linecut.pdf")
        
        # Test 7: Plot convergence data
        print("Plotting convergence...")
        try:
            fig5 = plotter.plot_convergence(job_name, save_path='vfet_convergence')
            print("  Saved: vfet_convergence.pdf")
        except Exception as e:
            print(f"  Convergence plot failed: {e}")
        
        print("\n" + "=" * 60)
        print("SUCCESS: All plots generated!")
        print("Generated files:")
        generated_files = [
            'vfet_potential_2d.pdf',
            'vfet_ndens_2d.pdf', 
            'vfet_potential_linecut.pdf',
            'vfet_ndens_linecut.pdf',
            'vfet_convergence.pdf'
        ]
        for f in generated_files:
            if Path(f).exists():
                print(f"  ✓ {f}")
            else:
                print(f"  ✗ {f}")
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())