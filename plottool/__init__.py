"""
Cybear Scientific Plotting Toolkit

A publication-quality plotting toolkit for Cybear semiconductor device simulation results.
Redesigned following your MATLAB workflow with separate data loading and plotting modules.

Main Components:
- CybearDataLoader: FlatBuffers data conversion using flott-cli (like MATLAB init_flott/load_flott)
- CybearPlotter: Publication-quality plotting with LaTeX labels  
- LabelDictionary: LaTeX label and unit management

Example Usage:
    from plottool import CybearPlotter
    
    plotter = CybearPlotter()
    
    # Plot 2D potential field
    plotter.plot_2d_field('vfet_test', 'SS.fbs', 'pot')
    
    # Plot electron density line cut
    plotter.plot_line_cut('vfet_test', 'SS.fbs', 'ndens', direction='x')
"""

from .data_loader import CybearDataLoader
from .cybear_plotter import CybearPlotter  
from .label_dict import LabelDictionary

__version__ = "1.0.0"
__author__ = "Cybear Project - Chenyang Yu"

# Main interfaces
def plot_2d_field(job_name: str, fbs_name: str, field_name: str, **kwargs):
    """
    Quick interface to plot 2D field distribution
    
    Parameters:
    -----------
    job_name : str
        Job directory name (e.g., 'vfet_test', 'nw_2D')
    fbs_name : str  
        FBS file name (e.g., 'SS.fbs')
    field_name : str
        Field variable name ('pot', 'ndens', etc.)
    **kwargs
        Additional plotting options
    """
    plotter = CybearPlotter()
    return plotter.plot_2d_field(job_name, fbs_name, field_name, **kwargs)

def plot_line_cut(job_name: str, fbs_name: str, field_name: str, **kwargs):
    """
    Quick interface to plot 1D line cut
    
    Parameters:
    -----------
    job_name : str
        Job directory name
    fbs_name : str
        FBS file name
    field_name : str
        Field variable name
    **kwargs
        Additional plotting options (direction, cut_position, etc.)
    """
    plotter = CybearPlotter()
    return plotter.plot_line_cut(job_name, fbs_name, field_name, **kwargs)