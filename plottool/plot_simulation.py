#!/usr/bin/env python3
"""
Practical plotting script for Cybear simulation results
Usage: python3 plot_simulation.py <job_name> [options]
"""

import sys
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

# Add plottool to path
sys.path.insert(0, str(Path(__file__).parent))

from cybear_plotter import CybearPlotter
from data_loader import CybearDataLoader

def main():
    parser = argparse.ArgumentParser(description='Plot Cybear simulation results')
    parser.add_argument('job_name', help='Job name (e.g., vfet_test)')
    parser.add_argument('--fbs', default='SS.fbs', help='FBS file name (default: SS.fbs)')
    parser.add_argument('--style', default='ieee', choices=['nature', 'ieee', 'retro', 'high-vis'],
                       help='Plot style (default: ieee)')
    parser.add_argument('--output-dir', default='.', help='Output directory (default: current)')
    parser.add_argument('--prefix', default='', help='Output file prefix')
    parser.add_argument('--fields', nargs='+', default=['pot', 'ndens'],
                       help='Fields to plot (default: pot ndens)')
    parser.add_argument('--plots', nargs='+', default=['2d', 'linecut'],
                       choices=['2d', 'linecut', 'convergence'],
                       help='Plot types (default: 2d linecut)')

    args = parser.parse_args()

    # Initialize plotter with chosen style
    if args.style == 'nature':
        plotter = CybearPlotter(style=['science', 'nature'])
    elif args.style == 'ieee':
        plotter = CybearPlotter(style=['science', 'ieee'])
    elif args.style == 'retro':
        plotter = CybearPlotter(style=['science', 'retro'])
    elif args.style == 'high-vis':
        plotter = CybearPlotter(style=['science', 'high-vis'])

    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    prefix = f"{args.prefix}_" if args.prefix else ""

    print(f"Plotting {args.job_name}/{args.fbs} with {args.style} style")
    print(f"Fields: {args.fields}")
    print(f"Plot types: {args.plots}")
    print(f"Output directory: {output_dir}")
    print("=" * 60)

    generated_files = []

    # Plot 2D fields
    if '2d' in args.plots:
        for field in args.fields:
            try:
                output_path = output_dir / f"{prefix}{args.job_name}_{field}_2d"
                fig = plotter.plot_2d_field(args.job_name, args.fbs, field,
                                          save_path=str(output_path))
                generated_files.append(f"{output_path}.pdf")
                print(f"✓ 2D {field}: {output_path}.pdf")
            except Exception as e:
                print(f"✗ Failed 2D {field}: {e}")

    # Plot line cuts
    if 'linecut' in args.plots:
        for field in args.fields:
            for direction in ['x', 'y']:
                try:
                    output_path = output_dir / f"{prefix}{args.job_name}_{field}_{direction}_linecut"
                    fig = plotter.plot_line_cut(args.job_name, args.fbs, field,
                                              direction=direction, save_path=str(output_path))
                    generated_files.append(f"{output_path}.pdf")
                    print(f"✓ {direction.upper()}-cut {field}: {output_path}.pdf")
                except Exception as e:
                    print(f"✗ Failed {direction.upper()}-cut {field}: {e}")

    # Plot convergence (if available)
    if 'convergence' in args.plots:
        try:
            output_path = output_dir / f"{prefix}{args.job_name}_convergence"
            fig = plotter.plot_convergence(save_path=str(output_path))
            generated_files.append(f"{output_path}.pdf")
            print(f"✓ Convergence: {output_path}.pdf")
        except Exception as e:
            print(f"✗ Failed convergence: {e}")

    print("=" * 60)
    print(f"Generated {len(generated_files)} files:")
    for file in generated_files:
        print(f"  {file}")

if __name__ == "__main__":
    main()
