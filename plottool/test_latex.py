#!/usr/bin/env python3
"""
Quick test of LaTeX rendering in plottool
"""

import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import os

# Change to cybear root directory
cybear_root = Path(__file__).parent.parent
os.chdir(cybear_root)

# Add plottool to path
sys.path.insert(0, str(Path(__file__).parent))

import cybear_plotter
import data_loader
import label_dict

def test_latex():
    print("Testing LaTeX rendering...")
    
    # Initialize plotter with LaTeX enabled
    plotter = cybear_plotter.CybearPlotter()
    
    # Test with actual data
    job_name = 'vfet_test'
    fbs_name = 'SS.fbs'
    
    try:
        # Generate one plot to test LaTeX
        print("Plotting 2D potential with LaTeX labels...")
        fig = plotter.plot_2d_field(job_name, fbs_name, 'pot',
                                   save_path='latex_test_potential',
                                   colormap='RdBu_r')
        print("✓ LaTeX test successful!")
        print("Generated: latex_test_potential.pdf")
        
        return True
        
    except Exception as e:
        print(f"✗ LaTeX test failed: {e}")
        return False

if __name__ == "__main__":
    success = test_latex()
    sys.exit(0 if success else 1)