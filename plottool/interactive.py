#!/usr/bin/env python3
"""
Interactive plotter for Cybear - Menu-driven interface
Just run this file and follow the prompts!
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from core import Plotter, DataLoader

def list_jobs():
    """List available simulation jobs"""
    # Find the project root
    current = Path.cwd()
    if (current / "jobs").exists():
        jobs_dir = current / "jobs"
    elif (current.parent / "jobs").exists():
        jobs_dir = current.parent / "jobs"
    else:
        return []

    jobs = []
    for job_path in jobs_dir.iterdir():
        if job_path.is_dir() and (job_path / "run").exists():
            jobs.append(job_path.name)
    return sorted(jobs)

def list_fbs_files(job_name):
    """List FBS files in a job"""
    # Find the project root
    current = Path.cwd()
    if (current / "jobs").exists():
        root = current
    elif (current.parent / "jobs").exists():
        root = current.parent
    else:
        return []

    run_dir = root / f"jobs/{job_name}/run"
    if not run_dir.exists():
        return []

    fbs_files = [f.name for f in run_dir.glob("*.fbs")]
    return sorted(fbs_files)

def main():
    print("\n" + "="*60)
    print(" " * 15 + "🔬 CYBEAR PLOTTER 🔬")
    print(" " * 18 + "Interactive Mode")
    print("="*60)
    print("   Publication-quality plots for semiconductor simulations")
    print("="*60 + "\n")

    # Step 1: Select job
    jobs = list_jobs()
    if not jobs:
        print("✗ No simulation jobs found in jobs/ directory")
        return

    print("Available simulation jobs:")
    print("-" * 30)
    for i, job in enumerate(jobs, 1):
        print(f"  {i}. {job}")

    while True:
        try:
            choice = input(f"\nSelect job (1-{len(jobs)}): ")
            job_idx = int(choice) - 1
            if 0 <= job_idx < len(jobs):
                selected_job = jobs[job_idx]
                break
            print("Invalid choice, please try again.")
        except (ValueError, KeyboardInterrupt):
            print("\n\nCancelled.")
            return

    # Step 2: Select FBS file
    fbs_files = list_fbs_files(selected_job)
    if not fbs_files:
        print(f"✗ No FBS files found in jobs/{selected_job}/run/")
        return

    print(f"\nFBS files in {selected_job}:")
    print("-" * 30)
    for i, fbs in enumerate(fbs_files, 1):
        print(f"  {i}. {fbs}")

    while True:
        try:
            choice = input(f"\nSelect file (1-{len(fbs_files)}): ")
            file_idx = int(choice) - 1
            if 0 <= file_idx < len(fbs_files):
                selected_file = fbs_files[file_idx]
                break
            print("Invalid choice, please try again.")
        except (ValueError, KeyboardInterrupt):
            print("\n\nCancelled.")
            return

    # Step 3: Select what to plot
    print("\nWhat would you like to plot?")
    print("-" * 30)
    print("  1. Auto (smart detection - recommended)")
    print("  2. All available fields")
    print("  3. Let me check what's available first")

    while True:
        try:
            choice = input("\nSelect option (1-3): ")
            if choice == '1':
                fields = 'auto'
                break
            elif choice == '2':
                fields = 'all'
                break
            elif choice == '3':
                # Show available fields
                print("\nChecking available fields...")
                loader = DataLoader()
                try:
                    # Suppress flott-cli warnings by redirecting stderr temporarily
                    import os, sys
                    old_stderr = sys.stderr
                    sys.stderr = open(os.devnull, 'w')

                    data = loader.load_data(selected_job, selected_file)

                    # Restore stderr
                    sys.stderr.close()
                    sys.stderr = old_stderr

                    smart_fields = loader.detect_fields(data, mode='smart')
                    all_fields = loader.detect_fields(data, mode='all')

                    print("\n" + "="*50)
                    print("📊 FIELD ANALYSIS")
                    print("="*50)

                    print("\n🎯 Smart detection (recommended):")
                    if smart_fields:
                        for field in smart_fields:
                            print(f"    • {field}")
                    else:
                        print("    (no priority physics fields found)")

                    print("\n📋 All available fields:")
                    if all_fields:
                        for field in all_fields:
                            marker = "★" if field in smart_fields else " "
                            print(f"    {marker} {field}")
                    else:
                        print("    (no plottable fields found)")

                    print("\n" + "-"*50)
                    print("Options:")
                    print("  a - Plot auto (smart detection)")
                    print("  l - Plot all fields")
                    print("  s - Select specific fields")
                    print("  c - Cancel")

                    while True:
                        choice2 = input("\nYour choice [a/l/s/c]: ").strip().lower()
                        if choice2 == 'a':
                            fields = 'auto'
                            break
                        elif choice2 == 'l':
                            fields = 'all'
                            break
                        elif choice2 == 's':
                            # Let user select specific fields
                            print("\nEnter field names separated by spaces:")
                            print(f"Available: {', '.join(all_fields)}")
                            field_input = input("\nFields: ").strip()
                            if field_input:
                                fields = field_input.split()
                                # Validate fields
                                invalid = [f for f in fields if f not in all_fields]
                                if invalid:
                                    print(f"✗ Invalid fields: {', '.join(invalid)}")
                                    continue
                                break
                            else:
                                print("No fields entered.")
                        elif choice2 == 'c':
                            print("Cancelled.")
                            return
                        else:
                            print("Invalid choice. Please enter a, l, s, or c.")
                    break
                except Exception as e:
                    # Try to give a cleaner error message
                    error_msg = str(e)
                    if 'flott-cli conversion failed' in error_msg:
                        print(f"\n✗ Error: Could not load simulation data")
                        print(f"   This might be a temporary issue with the data file.")
                        print(f"   Try option 1 (Auto) or 2 (All) instead.")
                    else:
                        print(f"✗ Error: {error_msg}")

                    # Offer to continue anyway
                    cont = input("\nContinue with auto mode? [Y/n]: ").strip().lower()
                    if cont != 'n':
                        fields = 'auto'
                        break
                    else:
                        return
            else:
                print("Invalid choice, please try again.")
        except KeyboardInterrupt:
            print("\n\nCancelled.")
            return

    # Step 4: Output directory
    print("\n📁 Output directory options:")
    print("  . (current directory)")
    print("  plots/ (create plots folder)")
    print("  Or enter custom path")
    output = input("\nOutput directory [.]: ").strip() or '.'

    # Create directory if needed
    from pathlib import Path
    output_path = Path(output)
    if not output_path.exists():
        create = input(f"\nDirectory '{output}' doesn't exist. Create it? [Y/n]: ").strip().lower()
        if create != 'n':
            output_path.mkdir(parents=True, exist_ok=True)
            print(f"✓ Created directory: {output}")

    # Step 5: Plot!
    print("\n" + "="*60)
    print("Ready to plot!")
    print(f"  Job: {selected_job}")
    print(f"  File: {selected_file}")
    print(f"  Fields: {fields}")
    print(f"  Output: {output}")
    print("="*60)

    confirm = input("\nProceed? [Y/n]: ").strip().lower()
    if confirm == 'n':
        print("Cancelled.")
        return

    # Create plotter and run
    print("\n🚀 Starting plot generation...\n")
    plotter = Plotter(style='ieee')
    plotter.plot(selected_job, selected_file, fields, output)

    # Show summary
    print("\n" + "="*60)
    print("📊 PLOTTING COMPLETE!")
    print("="*60)
    print(f"✓ Job: {selected_job}")
    print(f"✓ File: {selected_file}")
    print(f"✓ Mode: {fields if isinstance(fields, str) else 'custom'}")
    if isinstance(fields, list):
        print(f"✓ Fields: {', '.join(fields)}")
    print(f"✓ Output: {output}")
    print("\n🎉 Your PDFs are ready!")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n✗ Interrupted by user")
    except Exception as e:
        print(f"\n✗ Error: {e}")