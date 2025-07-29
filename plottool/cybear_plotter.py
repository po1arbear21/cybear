"""
Publication-Quality Plotting Module for Cybear Simulation Results

Uses SciencePlots for professional journal-ready styling.
Integrates with CybearDataLoader for seamless data handling.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import warnings

# Import SciencePlots for publication-quality styling
try:
    import scienceplots
    SCIENCEPLOTS_AVAILABLE = True
except ImportError:
    SCIENCEPLOTS_AVAILABLE = False
    warnings.warn("SciencePlots not available. Install with: pip install scienceplots")

from data_loader import CybearDataLoader
from label_dict import LabelDictionary

class CybearPlotter:
    """
    High-level plotting interface for Cybear simulation results.
    
    Uses SciencePlots for professional journal styling:
    - Nature, IEEE, and other journal standards
    - Automatic LaTeX integration
    - Color-blind friendly palettes
    - Publication-ready vector output
    
    Usage:
        plotter = CybearPlotter()
        plotter.plot_2d_field('vfet_test', 'SS.fbs', 'pot')
        plotter.plot_line_cut('vfet_test', 'SS.fbs', 'ndens', direction='x')
    """
    
    def __init__(self, style: str = "nature", figsize: Tuple[float, float] = (3.5, 2.5)):
        """
        Initialize plotter with SciencePlots styling
        
        Parameters:
        -----------
        style : str
            Journal style ('nature', 'ieee', 'science', or custom list)
        figsize : tuple
            Figure size in inches (3.5, 2.5) ≈ 88mm width for Nature single column
        """
        self.style = style
        self.figsize = figsize
        self.data_loader = CybearDataLoader()
        self.labels = LabelDictionary()
        
        # Apply SciencePlots styling
        self._apply_style()
        
    def _apply_style(self):
        """Apply SciencePlots styling - much cleaner than manual config!"""
        
        if not SCIENCEPLOTS_AVAILABLE:
            warnings.warn("SciencePlots not available, using matplotlib defaults")
            # Minimal fallback configuration
            plt.rcParams.update({
                'figure.figsize': self.figsize,
                'font.size': 8,
                'axes.grid': True,
                'grid.alpha': 0.3,
            })
            return
        
        # Apply SciencePlots style - handles all the LaTeX, fonts, colors automatically!
        if self.style == "nature":
            plt.style.use(['science', 'nature'])
        elif self.style == "ieee":
            plt.style.use(['science', 'ieee'])
        elif self.style == "science":
            plt.style.use('science')
        elif isinstance(self.style, list):
            plt.style.use(self.style)  # Custom style combination
        else:
            plt.style.use(['science', 'nature'])  # Default fallback
        
        # Override figure size (SciencePlots doesn't set this)
        plt.rcParams['figure.figsize'] = self.figsize
        
        # Enable LaTeX rendering for proper SciencePlots styling
        plt.rcParams['text.usetex'] = True
    
    def plot_2d_field(self, job_name: str, fbs_name: str, field_name: str,
                     frame: int = 0, colormap: str = 'RdBu_r',
                     xlim: Optional[Tuple[float, float]] = None,
                     ylim: Optional[Tuple[float, float]] = None,
                     vmin: Optional[float] = None,
                     vmax: Optional[float] = None,
                     save_path: Optional[str] = None, **kwargs) -> plt.Figure:
        """
        Plot 2D field distribution (potential, electron density, etc.)
        
        Parameters:
        -----------
        job_name : str
            Job directory name (e.g., 'vfet_test', 'nw_2D')
        fbs_name : str
            FBS file name (e.g., 'SS.fbs')
        field_name : str
            Field variable name ('pot', 'ndens', etc.)
        frame : int
            Frame index for time-dependent data (default: 0)
        colormap : str
            Matplotlib colormap name (default: 'RdBu_r')
        xlim : tuple of float, optional
            X-axis limits as (xmin, xmax) in nm (default: None, auto)
        ylim : tuple of float, optional
            Y-axis limits as (ymin, ymax) in nm (default: None, auto)
        vmin : float, optional
            Minimum value for colorbar (default: None, auto)
        vmax : float, optional
            Maximum value for colorbar (default: None, auto)
        save_path : str, optional
            Path to save figure
        **kwargs
            Additional arguments for imshow/contourf
            
        Returns:
        --------
        plt.Figure
            The created figure
        """
        # Load simulation data
        data = self.data_loader.load_simulation_data(job_name, fbs_name)
        field_data = self.data_loader.get_2d_field_data(data, field_name, frame)
        
        # Extract arrays
        x = field_data['x']
        y = field_data['y']
        field = field_data['field']
        
        # Create figure
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # Choose visualization method based on field type
        if 'dens' in field_name.lower():
            # Use logarithmic normalization for densities to show actual values on colorbar
            from matplotlib.colors import LogNorm
            # Avoid zero/negative values for log scale
            field_safe = np.maximum(field, 1e-20)
            
            # Create LogNorm with optional vmin/vmax
            if vmin is not None or vmax is not None:
                norm = LogNorm(vmin=vmin, vmax=vmax)
            else:
                norm = LogNorm()
                
            im = ax.imshow(field_safe.T, extent=[x.min(), x.max(), y.min(), y.max()],
                          origin='lower', aspect='auto', cmap=colormap, 
                          norm=norm, **kwargs)
            # Colorbar with logarithmic scale showing actual values
            cbar = plt.colorbar(im, ax=ax, shrink=0.8)
            # Label shows the actual quantity (not log), but colorbar scale is logarithmic
            symbol = self.labels._quantities.get(field_name.lower(), {}).get('symbol', field_name)
            unit = self.labels._quantities.get(field_name.lower(), {}).get('unit', '')
            if unit:
                cbar.set_label(r'$' + symbol + r'$ ($\mathrm{' + unit + r'}$)', rotation=270, labelpad=15)
            else:
                cbar.set_label(r'$' + symbol + r'$', rotation=270, labelpad=15)
        else:
            # Linear scale for other fields
            im = ax.imshow(field.T, extent=[x.min(), x.max(), y.min(), y.max()],
                          origin='lower', aspect='auto', cmap=colormap, 
                          vmin=vmin, vmax=vmax, **kwargs)
            # Colorbar
            cbar = plt.colorbar(im, ax=ax, shrink=0.8)
            cbar.set_label(self.labels.get_label(field_name), rotation=270, labelpad=15)
        
        # Labels using dictionary
        xlabel, ylabel = self.labels.get_labels('x', 'y')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(self.labels.get_title(field_name))
        
        # Set axis limits if provided
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        
        # Formatting
        plt.tight_layout()
        
        if save_path:
            self._save_figure(fig, save_path)
        
        return fig
    
    def plot_line_cut(self, job_name: str, fbs_name: str, field_name: str,
                     direction: str = 'x', cut_position: Optional[float] = None,
                     frame: int = 0, log_scale: Optional[bool] = None,
                     save_path: Optional[str] = None, **kwargs) -> plt.Figure:
        """
        Plot 1D line cut through 2D field data
        
        Parameters:
        -----------
        job_name : str
            Job directory name
        fbs_name : str
            FBS file name
        field_name : str
            Field variable name
        direction : str
            Cut direction ('x' or 'y')
        cut_position : float, optional
            Position of cut (default: center)
        frame : int
            Frame index for time-dependent data
        log_scale : bool, optional
            Use logarithmic scale (auto-detect for densities)
        save_path : str, optional
            Path to save figure
        **kwargs
            Additional arguments for plot()
            
        Returns:
        --------
        plt.Figure
            The created figure
        """
        # Load simulation data
        data = self.data_loader.load_simulation_data(job_name, fbs_name)
        line_data = self.data_loader.get_line_cut(data, field_name, direction, cut_position, frame)
        
        # Extract arrays
        coord = line_data['coord']
        field = line_data['field']
        coord_name = line_data['coord_name']
        cut_pos = line_data['cut_position']
        
        # Auto-detect log scale for densities
        if log_scale is None:
            log_scale = 'dens' in field_name.lower()
        
        # Create figure
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # Plot with appropriate scale
        if log_scale:
            ax.semilogy(coord, np.maximum(field, 1e-20), **kwargs)
        else:
            ax.plot(coord, field, **kwargs)
        
        # Labels using dictionary
        xlabel, ylabel = self.labels.get_labels(coord_name, field_name)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        
        # Title with cut information
        other_coord = 'y' if coord_name == 'x' else 'x'
        title = f"{self.labels.get_title(field_name)} at {other_coord} = {cut_pos:.1f} nm"
        ax.set_title(title)
        
        # Formatting
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if save_path:
            self._save_figure(fig, save_path)
        
        return fig
    
    def plot_iv_characteristics(self, job_name: str, fbs_name: str,
                              voltage_var: str = 'v_drn', current_var: str = 'i_drn',
                              save_path: Optional[str] = None, **kwargs) -> plt.Figure:
        """
        Plot I-V characteristics from simulation data
        
        Parameters:
        -----------
        job_name : str
            Job directory name  
        fbs_name : str
            FBS file name
        voltage_var : str
            Voltage variable name (default: 'v_drn')
        current_var : str
            Current variable name (default: 'i_drn')
        save_path : str, optional
            Path to save figure
        **kwargs
            Additional arguments for plot()
            
        Returns:
        --------
        plt.Figure
            The created figure
        """
        # Load simulation data
        data = self.data_loader.load_simulation_data(job_name, fbs_name)
        
        if voltage_var not in data or current_var not in data:
            available = list(data.keys())
            raise ValueError(f"Variables '{voltage_var}', '{current_var}' not found. Available: {available}")
        
        voltage = data[voltage_var]
        current = data[current_var]
        
        # Create figure
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # Plot I-V curve
        ax.plot(voltage, current, 'o-', linewidth=1.2, markersize=4, **kwargs)
        
        # Labels using dictionary
        xlabel, ylabel = self.labels.get_labels(voltage_var, current_var)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title('I-V Characteristics')
        
        # Formatting
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if save_path:
            self._save_figure(fig, save_path)
        
        return fig
    
    def plot_convergence(self, job_name: str, save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot Newton-Raphson convergence from log file
        
        Parameters:
        -----------
        job_name : str
            Job directory name
        save_path : str, optional
            Path to save figure
            
        Returns:
        --------
        plt.Figure
            The created figure
        """
        # Read log file directly (not from FBS)
        log_path = Path.cwd() / "jobs" / job_name / "run" / "stdout.txt"
        
        if not log_path.exists():
            raise FileNotFoundError(f"Log file not found: {log_path}")
        
        # Parse convergence data
        convergence_data = []
        with open(log_path, 'r') as f:
            for line in f:
                line = line.strip()
                if 'NLPE:' in line:
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            iteration = int(parts[1])
                            residual = float(parts[2])
                            convergence_data.append((iteration, residual))
                        except (ValueError, IndexError):
                            continue
        
        if not convergence_data:
            raise ValueError(f"No convergence data found in {log_path}")
        
        # Extract arrays
        iterations, residuals = zip(*convergence_data)
        
        # Create figure
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # Plot convergence
        ax.semilogy(iterations, residuals, 'o-', linewidth=1.2, markersize=4)
        
        # Labels using dictionary
        xlabel, ylabel = self.labels.get_labels('iteration', 'residual')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title('Newton-Raphson Convergence')
        
        # Formatting
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if save_path:
            self._save_figure(fig, save_path)
        
        return fig
    
    def plot_comparison(self, data_list: List[Dict], x_var: str, y_var: str,
                       labels: Optional[List[str]] = None,
                       save_path: Optional[str] = None, **kwargs) -> plt.Figure:
        """
        Plot comparison of multiple datasets
        
        Parameters:
        -----------
        data_list : list of dict
            List of simulation data dictionaries
        x_var, y_var : str
            Variable names for x and y axes
        labels : list of str, optional
            Legend labels for each dataset
        save_path : str, optional
            Path to save figure
        **kwargs
            Additional arguments for plot()
            
        Returns:
        --------
        plt.Figure
            The created figure
        """
        fig, ax = plt.subplots(figsize=self.figsize)
        
        for i, data in enumerate(data_list):
            if x_var not in data or y_var not in data:
                warnings.warn(f"Variables '{x_var}', '{y_var}' not found in dataset {i}")
                continue
            
            x = data[x_var]
            y = data[y_var]
            
            label = labels[i] if labels and i < len(labels) else f"Dataset {i+1}"
            ax.plot(x, y, 'o-', linewidth=1.2, markersize=4, label=label, **kwargs)
        
        # Labels using dictionary
        xlabel, ylabel = self.labels.get_labels(x_var, y_var)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f'{self.labels.get_title(y_var)} vs {self.labels.get_title(x_var)}')
        
        # Legend and formatting
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if save_path:
            self._save_figure(fig, save_path)
        
        return fig
    
    def _save_figure(self, fig: plt.Figure, save_path: str):
        """
        Save figure with publication-quality settings
        
        Parameters:
        -----------
        fig : plt.Figure
            Figure to save
        save_path : str
            Output path (will add .pdf extension if not present)
        """
        save_path = Path(save_path)
        if not save_path.suffix:
            save_path = save_path.with_suffix('.pdf')
        
        fig.savefig(save_path,
                   bbox_inches='tight',
                   pad_inches=0.05,
                   dpi=300,
                   facecolor='white',
                   edgecolor='none')
        
        print(f"Figure saved: {save_path}")
    
    def set_custom_label(self, var_name: str, symbol: str, unit: str, description: str):
        """
        Add custom label to the dictionary
        
        Parameters:
        -----------
        var_name : str
            Variable name
        symbol : str
            LaTeX symbol (without $)
        unit : str
            Physical unit
        description : str
            Descriptive name
        """
        self.labels.add_custom_label(var_name, symbol, unit, description)
    
    def list_available_data(self, job_name: str) -> List[str]:
        """
        List available .fbs files in a job directory
        
        Parameters:
        -----------
        job_name : str
            Job directory name
            
        Returns:
        --------
        list
            List of available .fbs file names
        """
        return self.data_loader.list_available_data(job_name)