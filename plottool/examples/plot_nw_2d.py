#!/usr/bin/env python3
"""
Example: Plotting 2D Nanowire FET Results

This script demonstrates how to use the Cybear plotting toolkit
to create publication-quality figures from simulation results.

Usage:
    python plot_nw_2d.py [path_to_simulation_output]
"""

import sys
from pathlib import Path

# Add the plotting directory to path
plotting_dir = Path(__file__).parent.parent
sys.path.insert(0, str(plotting_dir))

from cybear_plot import CybearPlot
from data_reader import SimulationData
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

def main():
    # Default to the nw_2D simulation results
    if len(sys.argv) > 1:
        data_path = sys.argv[1]
    else:
        # Go up to cybear root directory
        cybear_root = Path(__file__).parent.parent.parent.parent
        data_path = cybear_root / "jobs" / "nw_2D" / "run"
    
    print(f"Loading simulation data from: {data_path}")
    
    try:
        # Load simulation data
        data = SimulationData(data_path)
        data.info()
        
        # Initialize plotter with Nature-style formatting
        plotter = CybearPlot(style="cybear_nature")
        
        # Example 1: Plot convergence data if available
        convergence = data.get_convergence_data()
        if convergence:
            plot_convergence(plotter, convergence, save_path="convergence_plot")
        
        # Example 2: Demonstrate with synthetic I-V data
        # (In real usage, this would come from data.get_iv_data())
        plot_synthetic_iv(plotter)
        
        # Example 3: Show configuration summary
        config = data.get_config()
        if config:
            print("\nSimulation Configuration Summary:")
            for section, params in config.items():
                print(f"[{section}]")
                for key, value in params.items():
                    if isinstance(value, dict) and 'value' in value:
                        print(f"  {key} = {value['value']} {value['unit']}")
                    else:
                        print(f"  {key} = {value}")
                print()
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

def plot_convergence(plotter, convergence_data, save_path=None):
    """Plot Newton-Raphson convergence data"""
    print("Plotting convergence data...")
    
    # Extract NLPE convergence data
    nlpe_data = [entry for entry in convergence_data if entry['type'] == 'NLPE']
    
    if nlpe_data:
        iterations = [entry['iteration'] for entry in nlpe_data]
        residuals = [entry['residual'] for entry in nlpe_data]
        
        fig, ax = plt.subplots(figsize=plotter.figsize)
        ax.semilogy(iterations, residuals, 'o-', linewidth=1.2, markersize=4)
        ax.set_xlabel('Newton Iteration')
        ax.set_ylabel('Residual Norm')
        ax.set_title('NLPE Convergence')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plotter._save_figure(fig, save_path)
        
        plt.close(fig)  # Close instead of showing

def plot_synthetic_iv(plotter):
    """Create synthetic I-V plots to demonstrate the toolkit"""
    print("Creating synthetic I-V characteristics...")
    
    # Generate synthetic data similar to FET characteristics
    vds = np.linspace(0, 2.0, 50)  # Drain voltage [V]
    vgs_values = np.array([-0.5, 0.0, 0.5, 1.0, 1.5])  # Gate voltages [V]
    
    # Simple FET model for demonstration
    ids_data = np.zeros((len(vgs_values), len(vds)))
    
    for i, vgs in enumerate(vgs_values):
        # Simplified FET I-V model
        for j, vd in enumerate(vds):
            if vgs < 0.5:  # Below threshold
                ids_data[i, j] = 1e-12 * np.exp(vgs/0.1) * vd
            else:  # Above threshold
                vgs_eff = vgs - 0.5
                if vd < vgs_eff:  # Linear region
                    ids_data[i, j] = 1e-4 * vgs_eff * vd * (1 - vd/(2*vgs_eff))
                else:  # Saturation region
                    ids_data[i, j] = 1e-4 * vgs_eff**2 * 0.5
    
    # Plot I-V characteristics
    fig1 = plotter.plot_iv_characteristics(
        vds=vds,
        ids=ids_data, 
        vgs_values=vgs_values,
        save_path="iv_characteristics"
    )
    
    # Plot transfer characteristics
    vgs = np.linspace(-1.0, 2.0, 100)  # Gate voltage [V]
    vds_values = np.array([0.1, 0.5, 1.0, 2.0])  # Drain voltages [V]
    
    ids_transfer = np.zeros((len(vds_values), len(vgs)))
    
    for i, vd in enumerate(vds_values):
        for j, vg in enumerate(vgs):
            if vg < 0.5:  # Below threshold
                ids_transfer[i, j] = 1e-12 * np.exp(vg/0.1) * vd
            else:  # Above threshold
                vg_eff = vg - 0.5
                if vd < vg_eff:  # Linear region
                    ids_transfer[i, j] = 1e-4 * vg_eff * vd * (1 - vd/(2*vg_eff))
                else:  # Saturation region
                    ids_transfer[i, j] = 1e-4 * vg_eff**2 * 0.5
    
    fig2 = plotter.plot_transfer_characteristics(
        vgs=vgs,
        ids=ids_transfer,
        vds_values=vds_values,
        log_scale=True,
        save_path="transfer_characteristics"
    )
    
    # Example 2D field plot
    x = np.linspace(0, 100, 50)  # nm
    y = np.linspace(0, 20, 20)   # nm
    X, Y = np.meshgrid(x, y)
    
    # Synthetic potential field
    potential = np.sin(X/20) * np.exp(-Y/10) + 0.5 * (X/100)**2
    
    fig3 = plotter.plot_2d_field(
        x=x, y=y, field=potential,
        field_name="Electric Potential",
        unit="V",
        save_path="potential_field"
    )
    
    # Close all figures
    plt.close('all')
    
    print("Example plots created successfully!")
    print("Generated files:")
    print("  - iv_characteristics.pdf")
    print("  - transfer_characteristics.pdf") 
    print("  - potential_field.pdf")

if __name__ == "__main__":
    sys.exit(main())