#!/usr/bin/env python3
"""
Demo script showing different ways to use the plotter
"""

print("""
╔══════════════════════════════════════════════════════════════╗
║                    CYBEAR PLOTTER DEMO                      ║
╚══════════════════════════════════════════════════════════════╝

This plotter has THREE ways to use it:

1. QUICKPLOT (Edit & Run):
   - Edit quickplot.py
   - Change JOB, FILE, FIELDS variables
   - Run it (no arguments needed!)

2. INTERACTIVE (Menu-driven):
   - Run: python3 interactive.py
   - Follow the menus
   - Option 3 shows you all available fields!

3. COMMAND LINE (For scripts):
   - Import: from core import Plotter
   - Use: plotter.plot(job, file, fields, output)

═══════════════════════════════════════════════════════════════

Example field options:
  • 'auto'              → Smart detection (pot, ndens, Ex, Ey)
  • 'all'               → Everything available
  • ['pot', 'ndens']    → Specific fields only
  • ['Ex', 'Ey']        → Just electric field components

The plotter handles:
  ✓ 1D simulations (Schottky diodes)
  ✓ 2D simulations (MOSFETs, VFETs)
  ✓ Auto-detection of dimensionality
  ✓ Log scale for density fields
  ✓ Appropriate colormaps for each field type
  ✓ Publication-ready IEEE journal styling

Try it now:
  python3 interactive.py

""")