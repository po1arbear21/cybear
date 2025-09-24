# Cybear Plotter - Click & Run!

## 🚀 Quick Start (No Command Line!)

### Method 1: Edit & Run (Recommended)
1. Open `quickplot.py` in any text editor
2. Edit these 3 lines:
   ```python
   JOB = 'vfet_test'    # Your job name
   FILE = 'SS.fbs'      # Your file
   FIELDS = 'auto'      # What to plot
   ```
3. Save and double-click the file (or run `python3 quickplot.py`)
4. Done! PDFs are generated automatically.

### Method 2: Interactive Menu
Just run:
```bash
python3 interactive.py
```
Follow the prompts - no typing commands!

## 📊 Field Options

- `'auto'` - Smart detection of important physics fields (pot, ndens, Ex, Ey)
- `'all'` - Plot everything available
- `['pot', 'ndens']` - Specific fields only

## 📁 Structure

```
plottool/
├── quickplot.py      # ← MAIN: Edit & run this!
├── interactive.py    # Alternative: Menu-driven
├── core/
│   ├── loader.py    # Data handling
│   ├── plotter.py   # Plotting logic
│   └── labels.py    # LaTeX labels
└── README.md        # This file
```

## ✨ Features

- **Zero command-line** - Just click and run!
- **Smart defaults** - Knows what physicists want to plot
- **Auto-conversion** - Handles .fbs → .mat automatically
- **1D/2D support** - Works with all simulation types
- **Publication-ready** - IEEE journal styling with LaTeX

## 📝 Examples

**Plot VFET simulation:**
```python
# In quickplot.py, just change:
JOB = 'vfet_test'
FILE = 'SS.fbs'
# Then run it!
```

**Plot Schottky diode:**
```python
JOB = 'schottky_test'
FILE = 'SS_Test.fbs'
```

**Plot everything:**
```python
FIELDS = 'all'  # Instead of 'auto'
```

That's it! No more command lines, just edit and click! 🎉