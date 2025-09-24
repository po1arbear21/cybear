#!/usr/bin/env python3
"""
CLICK AND RUN PLOTTER FOR CYBEAR
Just edit the variables below and run this file!
No command line arguments needed.
"""

# ==============================================================================
#                           EDIT THESE SETTINGS
# ==============================================================================

JOB = 'schottky_test'           # Job name (folder in jobs/)
FILE = 'IV_Sweep.fbs'              # FBS file to plot

# What to plot:
#   'auto'  = Smart detection (pot, ndens, Ex, Ey, etc.)
#   'all'   = Everything available
#   ['pot', 'ndens']  = Specific fields

# FIELDS = 'auto'
FIELDS =  ['pot', 'ndens']

OUTPUT_DIR = '.'             # Where to save PDFs ('.' = current directory)

# Export options:
EXPORT_TIKZ = True           # Also export TikZ/LaTeX files for direct inclusion in papers!
                            # Creates: .tex files + standalone LaTeX for testing
                            # Requires: pip install tikzplotlib


# ==============================================================================
#                          DON'T EDIT BELOW THIS LINE
# ==============================================================================

if __name__ == "__main__":
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent))
    from core import Plotter

    # Create plotter and run
    plotter = Plotter(style='ieee', export_tikz=EXPORT_TIKZ)
    plotter.plot(JOB, FILE, FIELDS, OUTPUT_DIR)

    # That's it! Check your PDFs.
