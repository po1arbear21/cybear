# Cybear Scientific Plotting Toolkit

A publication-quality plotting toolkit for Cybear semiconductor device simulation results. Redesigned following your MATLAB workflow with separate data loading and plotting modules.

## 🎯 Overview

This toolkit provides two main components:
1. **Data Loading** (`CybearDataLoader`) - FlatBuffers conversion using `flott-cli` (equivalent to your MATLAB `init_flott`/`load_flott`)
2. **Plotting** (`CybearPlotter`) - Publication-quality plots with LaTeX labels and journal standards

## 🚀 Quick Start

### Basic Usage

```python
from plottool import CybearPlotter

# Initialize plotter
plotter = CybearPlotter()

# Plot 2D potential field
plotter.plot_2d_field('vfet_test', 'SS.fbs', 'pot', save_path='potential_plot')

# Plot electron density line cut  
plotter.plot_line_cut('vfet_test', 'SS.fbs', 'ndens', direction='x', save_path='density_cut')
```

### Example: vfet_test Analysis

```python
# Test the toolkit
python3 test_vfet.py
```

This generates:
- `vfet_potential_2d.pdf` - 2D potential distribution
- `vfet_ndens_2d.pdf` - 2D electron density (log scale)
- `vfet_potential_linecut.pdf` - 1D potential line cut
- `vfet_ndens_linecut.pdf` - 1D density line cut  
- `vfet_convergence.pdf` - Newton-Raphson convergence

## 📁 Architecture

### Data Loading Module (`data_loader.py`)
Handles FlatBuffers conversion exactly like your MATLAB workflow:

```python
from data_loader import CybearDataLoader

loader = CybearDataLoader()

# Convert .fbs to .mat (like MATLAB init_flott)
mat_file = loader.init_flott('jobs/vfet_test/run/SS.fbs')

# Load specific variables (like MATLAB load_flott)
data = loader.load_flott('jobs/vfet_test/run/SS.fbs', ['pot', 'ndens'])
```

### Label Dictionary (`label_dict.py`)
LaTeX symbols and units for scientific variables:

```python
from label_dict import LabelDictionary

labels = LabelDictionary()
xlabel, ylabel = labels.get_labels('x', 'ndens')
# Returns: ('$x$ (nm)', '$n$ (cm$^{-3}$)')
```

## 🔬 Successfully Tested Features

✅ **FlatBuffers Data Loading**: Uses your `flott-cli` workflow  
✅ **2D Field Visualization**: Potential and density plots with proper scaling  
✅ **Line Cut Analysis**: 1D cuts through 2D data  
✅ **Convergence Plotting**: Newton-Raphson iteration analysis  
✅ **LaTeX Label Dictionary**: Automatic scientific notation  
✅ **Publication-Quality Output**: Vector PDF with proper formatting  
✅ **Auto-Detection**: Coordinate arrays (`x_vertices`, `y_vertices`)  
✅ **Log Scale**: Automatic for density fields  

## 🧪 Test Results

**Successfully generated from `vfet_test/SS.fbs`**:
- All plots generated with proper LaTeX labels
- Vector PDF output ready for publication
- Handles your actual simulation data structure
- Compatible with your existing `flott-cli` workflow

---

**Ready for your perovskite VFET research! 🚀**