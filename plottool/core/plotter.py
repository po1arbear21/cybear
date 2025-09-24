"""
Main plotting module - handles all visualization
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from .loader import DataLoader
from .labels import format_axis_label, get_label

# Try to use SciencePlots for journal styling
try:
    import scienceplots
    plt.style.use(['science', 'ieee'])
    LATEX_AVAILABLE = True
except:
    LATEX_AVAILABLE = False

# Try to import tikzplotlib for TikZ export
try:
    import tikzplotlib
    TIKZ_AVAILABLE = True
except ImportError as e:
    TIKZ_AVAILABLE = False
    TIKZ_IMPORT_ERROR = str(e)

class Plotter:
    """Main plotter class - handles all visualization"""

    def __init__(self, style='ieee', export_tikz=False):
        self.loader = DataLoader()
        self.export_tikz = export_tikz
        self.apply_style(style)

        if export_tikz and not TIKZ_AVAILABLE:
            raise ImportError(f"TikZ export failed: {TIKZ_IMPORT_ERROR}")

    def apply_style(self, style):
        """Apply plot styling"""
        if LATEX_AVAILABLE and style:
            if style == 'nature':
                plt.style.use(['science', 'nature'])
            elif style == 'ieee':
                plt.style.use(['science', 'ieee'])
            else:
                plt.style.use('science')

        # Set figure size
        plt.rcParams['figure.figsize'] = (4, 3)

    def plot(self, job_name, fbs_file, fields='auto', output_dir='.'):
        """Main plotting function - the entry point"""

        print(f"\n{'='*50}")
        print(f"Plotting: {job_name}/{fbs_file}")
        print(f"{'='*50}\n")

        # Load data
        try:
            data = self.loader.load_data(job_name, fbs_file)
        except Exception as e:
            print(f"✗ Error loading data: {e}")
            return

        # FIRST: Check for IV curves (these don't need spatial coordinates)
        iv_curves = self.loader.detect_iv_curves(data)
        if iv_curves:
            print(f"📈 Found {len(iv_curves)} IV curve(s)!")
            # Generate output prefix
            fbs_base = Path(fbs_file).stem
            output_prefix = Path(output_dir) / f"{job_name}_{fbs_base}"
            Path(output_dir).mkdir(exist_ok=True)

            for v_field, i_field, label in iv_curves:
                self._plot_iv_curve(data, v_field, i_field, label, output_prefix)

        # Get coordinates for spatial plots
        x, y = self.loader.get_coordinates(data)

        # If no coordinates, check if we already plotted IV curves
        if x is None:
            if iv_curves:
                print("\n✅ IV curves plotted successfully!")
                return
            else:
                print("✗ Could not find coordinate arrays or IV data")
                return

        # Determine what fields to plot
        if fields == 'auto':
            fields_to_plot = self.loader.detect_fields(data, mode='smart')
            if fields_to_plot:
                print(f"Auto-detected fields: {', '.join(fields_to_plot)}")
        elif fields == 'all':
            fields_to_plot = self.loader.detect_fields(data, mode='all')
            print(f"Plotting all fields: {', '.join(fields_to_plot)}")
        elif isinstance(fields, list):
            fields_to_plot = fields
            print(f"Plotting specified fields: {', '.join(fields_to_plot)}")
        else:
            print("✗ Invalid fields specification")
            return

        if not fields_to_plot:
            print("✗ No plottable fields found")
            return

        # Create output directory if needed
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)

        # Generate output prefix
        fbs_base = Path(fbs_file).stem
        output_prefix = output_path / f"{job_name}_{fbs_base}"

        # Plot each field
        for field in fields_to_plot:
            if field not in data:
                print(f"  ✗ Field '{field}' not found")
                continue

            print(f"\n📊 Plotting {field}...")

            # Determine dimensionality
            dim = self.loader.get_dimensionality(data, field)

            if dim == '1D':
                self._plot_1d(data, field, x, output_prefix)
            elif dim == '2D' or dim == '2D-multi':
                self._plot_2d(data, field, x, y, output_prefix)
                self._plot_linecuts(data, field, x, y, output_prefix)
            else:
                print(f"  ✗ Cannot plot {field} (unknown dimensionality)")

        print(f"\n{'='*50}")
        print("✓ Done! Check your PDFs.")
        print(f"{'='*50}\n")

    def _plot_1d(self, data, field, x, output_prefix):
        """Plot 1D field data"""
        field_data = data[field].flatten()

        fig, ax = plt.subplots()

        # Use log scale for density fields
        if 'dens' in field and field_data.min() > 0:
            ax.semilogy(x, field_data)
        else:
            ax.plot(x, field_data)

        # Labels
        ax.set_xlabel(format_axis_label('x'))
        ax.set_ylabel(format_axis_label(field))
        ax.grid(True, alpha=0.3)

        # Save
        output_file_base = f"{output_prefix}_{field}_1d"
        self._save_figure(fig, output_file_base)
        plt.close()

    def _plot_2d(self, data, field, x, y, output_prefix):
        """Plot 2D field distribution"""
        if y is None:
            print(f"  ✗ Cannot plot 2D without y coordinates")
            return

        field_data = data[field]

        # Handle 3D data (take first frame)
        if len(field_data.shape) == 3:
            field_data = field_data[:, :, 0]

        # Ensure correct orientation
        if field_data.shape[0] != len(x) or field_data.shape[1] != len(y):
            if field_data.shape[1] == len(x) and field_data.shape[0] == len(y):
                field_data = field_data.T

        fig, ax = plt.subplots()

        # Log scale for density fields
        use_log = 'dens' in field
        if use_log and field_data.min() > 0:
            field_data = np.log10(field_data)
            colorbar_label = f"log₁₀ {get_label(field)[0]} ({get_label(field)[1]})"
        else:
            colorbar_label = format_axis_label(field)

        # Select colormap
        cmap = self._get_colormap(field)

        # Create mesh plot
        X, Y = np.meshgrid(x, y, indexing='ij')
        im = ax.pcolormesh(X, Y, field_data, cmap=cmap, shading='auto')

        # Labels and colorbar
        ax.set_xlabel(format_axis_label('x'))
        ax.set_ylabel(format_axis_label('y'))
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(colorbar_label)

        # Save
        output_file_base = f"{output_prefix}_{field}_2d"
        self._save_figure(fig, output_file_base)
        plt.close()

    def _plot_linecuts(self, data, field, x, y, output_prefix):
        """Plot X and Y line cuts through center of 2D data"""
        if y is None:
            return

        field_data = data[field]

        # Handle 3D data
        if len(field_data.shape) == 3:
            field_data = field_data[:, :, 0]

        # Ensure correct orientation
        if field_data.shape[0] != len(x) or field_data.shape[1] != len(y):
            if field_data.shape[1] == len(x) and field_data.shape[0] == len(y):
                field_data = field_data.T

        use_log = 'dens' in field

        # X-cut (through center of y)
        fig, ax = plt.subplots()
        y_center = len(y) // 2
        line_data = field_data[:, y_center]

        if use_log and line_data.min() > 0:
            ax.semilogy(x, line_data)
        else:
            ax.plot(x, line_data)

        ax.set_xlabel(format_axis_label('x'))
        ax.set_ylabel(format_axis_label(field))
        ax.grid(True, alpha=0.3)

        output_file_base = f"{output_prefix}_{field}_xcut"
        self._save_figure(fig, output_file_base)
        plt.close()

        # Y-cut (through center of x)
        fig, ax = plt.subplots()
        x_center = len(x) // 2
        line_data = field_data[x_center, :]

        if use_log and line_data.min() > 0:
            ax.semilogy(y, line_data)
        else:
            ax.plot(y, line_data)

        ax.set_xlabel(format_axis_label('y'))
        ax.set_ylabel(format_axis_label(field))
        ax.grid(True, alpha=0.3)

        output_file_base = f"{output_prefix}_{field}_ycut"
        self._save_figure(fig, output_file_base)
        plt.close()

    def _save_figure(self, fig, output_file_base, dpi=300):
        """Save figure as PDF and optionally as TikZ

        Args:
            fig: matplotlib figure
            output_file_base: base filename without extension
            dpi: DPI for PDF output
        """
        # Always save PDF
        pdf_file = f"{output_file_base}.pdf"
        plt.tight_layout()
        plt.savefig(pdf_file, dpi=dpi, bbox_inches='tight')
        print(f"  ✓ Saved: {Path(pdf_file).name}")

        # Optionally save TikZ
        if self.export_tikz:
            tex_file = f"{output_file_base}.tex"
            import tikzplotlib
            tikzplotlib.save(
                tex_file,
                figure=fig,
                axis_width='\\figwidth',  # Use LaTeX variable for width
                axis_height='\\figheight',  # Use LaTeX variable for height
                textsize=10.0,
                float_format=".3g"  # Limit decimal places
            )
            print(f"  ✓ Saved TikZ: {Path(tex_file).name}")

            # Also create a standalone LaTeX file for testing
            standalone_file = f"{output_file_base}_standalone.tex"
            with open(standalone_file, 'w') as f:
                f.write(r"""\documentclass{standalone}
\usepackage{pgfplots}
\pgfplotsset{compat=1.17}
\usepackage{amsmath}
\newlength{\figwidth}
\newlength{\figheight}
\setlength{\figwidth}{8cm}
\setlength{\figheight}{6cm}
\begin{document}
\input{""" + Path(tex_file).name + r"""}
\end{document}
""")
            print(f"  ✓ Saved standalone: {Path(standalone_file).name}")

    def _get_colormap(self, field):
        """Get appropriate colormap for field type"""
        if 'pot' in field or 'phi' in field:
            return 'RdBu_r'
        elif 'dens' in field:
            return 'plasma'
        elif 'field' in field or field in ['Ex', 'Ey', 'Ez']:
            return 'viridis'
        elif 'current' in field or field in ['jn', 'jp']:
            return 'coolwarm'
        else:
            return 'viridis'

    def _plot_iv_curve(self, data, v_field, i_field, label, output_prefix):
        """Plot an IV curve"""
        voltage = data[v_field].flatten()
        current = data[i_field].flatten()

        # Convert current to appropriate units
        current_unit = 'A'

        # Auto-scale current
        max_current = np.abs(current).max()
        if max_current < 1e-9:
            current = current * 1e12  # Convert to pA
            current_unit = 'pA'
        elif max_current < 1e-6:
            current = current * 1e9   # Convert to nA
            current_unit = 'nA'
        elif max_current < 1e-3:
            current = current * 1e6   # Convert to µA
            current_unit = 'µA'
        elif max_current < 1:
            current = current * 1e3   # Convert to mA
            current_unit = 'mA'

        # Create figure
        fig, ax = plt.subplots()

        # Plot with markers and lines
        ax.plot(voltage, current, 'o-', markersize=4, linewidth=1.5)

        # Labels
        ax.set_xlabel('Voltage (V)')
        ax.set_ylabel(f'Current ({current_unit})')
        ax.grid(True, alpha=0.3)

        # Add zero lines for reference
        ax.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
        ax.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)

        # Title
        ax.set_title(label, fontsize=10)

        # Save
        clean_label = label.replace(' ', '_').replace('(', '').replace(')', '')
        output_file_base = f"{output_prefix}_{clean_label}"
        self._save_figure(fig, output_file_base)
        plt.close()