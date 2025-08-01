"""
Data Loading Module for Cybear Simulation Results

This module handles FlatBuffers conversion using flott-cli, similar to the MATLAB workflow.
Provides clean separation between data loading and plotting functionality.
"""

import os
import subprocess
import numpy as np
from pathlib import Path
from typing import Dict, Any, Optional, Union, List
import warnings
from scipy.io import loadmat
import tempfile

class CybearDataLoader:
    """
    Data loader for Cybear simulation results using flott-cli conversion.
    
    Mimics the MATLAB workflow:
    1. init_flott() -> Converts .fbs to .mat using flott-cli
    2. load_flott() -> Loads specific variables from converted .mat file
    """
    
    def __init__(self, flott_cli_path: str = "flott-cli"):
        """
        Initialize data loader
        
        Parameters:
        -----------
        flott_cli_path : str
            Path to flott-cli executable (default: "flott-cli" in PATH)
        """
        self.flott_cli_path = flott_cli_path
        self._verify_flott_cli()
        
    def _verify_flott_cli(self):
        """Verify that flott-cli is available"""
        try:
            result = subprocess.run([self.flott_cli_path, "--version"], 
                                  capture_output=True, text=True, check=True)
            print(f"Found flott-cli: {result.stdout.strip()}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError(f"flott-cli not found at: {self.flott_cli_path}")
    
    def init_flott(self, fbs_file: Union[str, Path], mat_file: Optional[Union[str, Path]] = None) -> Path:
        """
        Convert .fbs file to .mat using flott-cli (equivalent to MATLAB init_flott.m)
        
        Parameters:
        -----------
        fbs_file : str or Path
            Path to input .fbs file
        mat_file : str or Path, optional
            Output .mat file path (default: same name as .fbs with .mat extension)
            
        Returns:
        --------
        Path
            Path to created .mat file
        """
        fbs_path = Path(fbs_file)
        
        if not fbs_path.exists():
            raise FileNotFoundError(f"FBS file not found: {fbs_path}")
            
        if fbs_path.suffix.lower() != '.fbs':
            raise ValueError(f"Input file must have .fbs extension, got: {fbs_path.suffix}")
        
        # Generate output filename if not provided
        if mat_file is None:
            mat_path = fbs_path.with_suffix('.mat')
        else:
            mat_path = Path(mat_file)
            
        # Skip conversion if .mat file already exists and is newer
        if mat_path.exists() and mat_path.stat().st_mtime > fbs_path.stat().st_mtime:
            print(f"Using existing .mat file: {mat_path}")
            return mat_path
        
        # Construct flott-cli command (matching MATLAB version)
        command = [
            self.flott_cli_path,
            "-f", str(fbs_path),
            "export", "-a", "-r",
            "-o", str(mat_path)
        ]
        
        print(f"Converting: {fbs_path.name} -> {mat_path.name}")
        
        try:
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            
            # Verify output file was created
            if not mat_path.exists():
                raise RuntimeError(f"flott-cli completed but output file not created: {mat_path}")
                
            print(f"Successfully created: {mat_path}")
            return mat_path
            
        except subprocess.CalledProcessError as e:
            # Clean up partial file if it exists
            if mat_path.exists():
                mat_path.unlink()
                print(f"Cleaned up partial file: {mat_path}")
            
            error_msg = f"flott-cli error (code {e.returncode}):\n{e.stderr}"
            raise RuntimeError(error_msg)
    
    def load_flott(self, fbs_file: Union[str, Path], 
                   var_names: Optional[List[str]] = None,
                   mat_file: Optional[Union[str, Path]] = None) -> Dict[str, np.ndarray]:
        """
        Load variables from .fbs file (equivalent to MATLAB load_flott.m)
        
        Parameters:
        -----------
        fbs_file : str or Path
            Path to .fbs file
        var_names : list of str, optional
            List of variable names to load (default: load all)
        mat_file : str or Path, optional
            Pre-converted .mat file path
            
        Returns:
        --------
        dict
            Dictionary with variable names as keys and numpy arrays as values
        """
        fbs_path = Path(fbs_file)
        
        # Convert to .mat if needed
        if mat_file is None:
            mat_path = self.init_flott(fbs_path)
        else:
            mat_path = Path(mat_file)
            if not mat_path.exists():
                mat_path = self.init_flott(fbs_path, mat_file)
        
        # Load .mat file
        try:
            data = loadmat(str(mat_path), squeeze_me=True, struct_as_record=False)
        except Exception as e:
            raise RuntimeError(f"Failed to load .mat file {mat_path}: {e}")
        
        # Remove MATLAB metadata
        data = {k: v for k, v in data.items() if not k.startswith('__')}
        
        # Filter requested variables
        if var_names is not None:
            filtered_data = {}
            for var_name in var_names:
                if var_name in data:
                    filtered_data[var_name] = data[var_name]
                else:
                    warnings.warn(f"Variable '{var_name}' not found in {mat_path}")
                    filtered_data[var_name] = None
            return filtered_data
        
        return data
    
    def load_simulation_data(self, job_name: str, fbs_name: str) -> Dict[str, np.ndarray]:
        """
        Load simulation data by job name and FBS file name
        
        Parameters:
        -----------
        job_name : str
            Name of the job directory (e.g., 'vfet_test', 'nw_2D')
        fbs_name : str
            Name of the FBS file (e.g., 'SS.fbs')
            
        Returns:
        --------
        dict
            Dictionary containing simulation data
        """
        # Auto-locate cybear root and construct path
        current = Path.cwd()
        cybear_root = current
        
        # Find cybear root by looking for jobs and plottool directories
        for path in [current] + list(current.parents):
            if (path / 'jobs').exists() and (path / 'plottool').exists():
                cybear_root = path
                break
        
        fbs_path = cybear_root / "jobs" / job_name / "run" / fbs_name
        
        if not fbs_path.exists():
            raise FileNotFoundError(f"Simulation data not found: {fbs_path}")
        
        # Load all data
        data = self.load_flott(fbs_path)
        
        print(f"Loaded simulation data from {job_name}/{fbs_name}")
        print(f"Available variables: {list(data.keys())}")
        
        return data
    
    def get_2d_field_data(self, data: Dict[str, np.ndarray], 
                         field_name: str, frame: int = 0) -> Dict[str, np.ndarray]:
        """
        Extract 2D field data in format suitable for plotting
        
        Parameters:
        -----------
        data : dict
            Simulation data dictionary
        field_name : str
            Name of field variable (e.g., 'pot', 'ndens')
        frame : int
            Frame index for time-dependent data (default: 0)
            
        Returns:
        --------
        dict
            Dictionary with 'x', 'y', 'field' arrays ready for plotting
        """
        if field_name not in data:
            raise ValueError(f"Field '{field_name}' not found in data. Available: {list(data.keys())}")
        
        # Get coordinate arrays (try different possible names)
        x = data.get('x', data.get('x_vertices', data.get('x_coords', None)))
        y = data.get('y', data.get('y_vertices', data.get('y_coords', None)))
        
        if x is None or y is None:
            available_coords = [k for k in data.keys() if any(coord in k.lower() for coord in ['x', 'y', 'coord', 'vertex'])]
            raise ValueError(f"Coordinate arrays not found. Available coordinate-like variables: {available_coords}")
        
        # Get field data
        field = data[field_name]
        
        # Handle 3D data (take specified frame)
        if field.ndim == 3:
            field = field[:, :, frame]
        elif field.ndim != 2:
            raise ValueError(f"Expected 2D or 3D field data, got {field.ndim}D")
        
        return {
            'x': x,
            'y': y, 
            'field': field
        }
    
    def get_line_cut(self, data: Dict[str, np.ndarray], 
                    field_name: str, cut_direction: str = 'x',
                    cut_position: Optional[float] = None, frame: int = 0) -> Dict[str, np.ndarray]:
        """
        Extract 1D line cut from 2D field data
        
        Parameters:
        -----------
        data : dict
            Simulation data dictionary
        field_name : str
            Name of field variable
        cut_direction : str
            Direction of cut ('x' or 'y')
        cut_position : float, optional
            Position of cut (default: center)
        frame : int
            Frame index for time-dependent data
            
        Returns:
        --------
        dict
            Dictionary with coordinate and field arrays for 1D plot
        """
        field_data = self.get_2d_field_data(data, field_name, frame)
        
        x = field_data['x']
        y = field_data['y']
        field = field_data['field']
        
        if cut_direction.lower() == 'x':
            # Cut along x-direction at specified y position
            if cut_position is None:
                cut_idx = len(y) // 2  # Center
            else:
                cut_idx = np.argmin(np.abs(y - cut_position))
            
            return {
                'coord': x,
                'field': field[:, cut_idx],
                'coord_name': 'x',
                'cut_position': y[cut_idx]
            }
        
        elif cut_direction.lower() == 'y':
            # Cut along y-direction at specified x position
            if cut_position is None:
                cut_idx = len(x) // 2  # Center
            else:
                cut_idx = np.argmin(np.abs(x - cut_position))
            
            return {
                'coord': y,
                'field': field[cut_idx, :],
                'coord_name': 'y',
                'cut_position': x[cut_idx]
            }
        
        else:
            raise ValueError(f"cut_direction must be 'x' or 'y', got: {cut_direction}")
    
    def list_available_data(self, job_name: str) -> List[str]:
        """
        List available .fbs files in a job directory
        
        Parameters:
        -----------
        job_name : str
            Name of the job directory
            
        Returns:
        --------
        list
            List of available .fbs file names
        """
        # Auto-locate cybear root and construct path
        current = Path.cwd()
        cybear_root = current
        
        # Find cybear root by looking for jobs and plottool directories
        for path in [current] + list(current.parents):
            if (path / 'jobs').exists() and (path / 'plottool').exists():
                cybear_root = path
                break
        
        job_path = cybear_root / "jobs" / job_name / "run"
        
        if not job_path.exists():
            raise FileNotFoundError(f"Job directory not found: {job_path}")
        
        fbs_files = list(job_path.glob("*.fbs"))
        return [f.name for f in fbs_files]