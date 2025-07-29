"""
LaTeX Label Dictionary for Scientific Plotting

Provides consistent, publication-quality LaTeX labels with proper units
for semiconductor device simulation quantities.
"""

from typing import Dict, Tuple, Optional

class LabelDictionary:
    """
    Dictionary for LaTeX labels and units in semiconductor device simulation.

    Usage:
        labels = LabelDictionary()
        xlabel, ylabel = labels.get_labels('x', 'ndens')
        # Returns: ('$x$ (nm)', '$n$ (cm$^{-3}$)')
    """

    def __init__(self):
        """Initialize the label dictionary with common variables"""

        # Coordinate labels
        self._coordinates = {
            'x': {'symbol': r'x', 'unit': 'nm', 'description': 'Position'},
            'y': {'symbol': r'y', 'unit': 'nm', 'description': 'Position'},
            'z': {'symbol': r'z', 'unit': 'nm', 'description': 'Position'},
            'r': {'symbol': r'r', 'unit': 'nm', 'description': 'Radial position'},
        }

        # Physical quantities
        self._quantities = {
            # Electrostatics
            'pot': {'symbol': r'\phi', 'unit': 'V', 'description': 'Electrostatic potential'},
            'phi': {'symbol': r'\phi', 'unit': 'V', 'description': 'Electrostatic potential'},
            'efield': {'symbol': r'E', 'unit': 'V/cm', 'description': 'Electric field'},
            'ex': {'symbol': r'E_x', 'unit': 'V/cm', 'description': 'Electric field (x)'},
            'ey': {'symbol': r'E_y', 'unit': 'V/cm', 'description': 'Electric field (y)'},
            'ez': {'symbol': r'E_z', 'unit': 'V/cm', 'description': 'Electric field (z)'},

            # Carrier densities
            'ndens': {'symbol': r'n', 'unit': r'cm^{-3}', 'description': 'Electron density'},
            'pdens': {'symbol': r'p', 'unit': r'cm^{-3}', 'description': 'Hole density'},
            'rho': {'symbol': r'\rho', 'unit': r'C/cm^3', 'description': 'Charge density'},

            # Currents and current densities
            'current': {'symbol': r'I', 'unit': 'A', 'description': 'Current'},
            'jn': {'symbol': r'J_n', 'unit': r'A/cm^2', 'description': 'Electron current density'},
            'jp': {'symbol': r'J_p', 'unit': r'A/cm^2', 'description': 'Hole current density'},
            'jnx': {'symbol': r'J_{n,x}', 'unit': r'A/cm^2', 'description': 'Electron current density (x)'},
            'jny': {'symbol': r'J_{n,y}', 'unit': r'A/cm^2', 'description': 'Electron current density (y)'},
            'jpx': {'symbol': r'J_{p,x}', 'unit': r'A/cm^2', 'description': 'Hole current density (x)'},
            'jpy': {'symbol': r'J_{p,y}', 'unit': r'A/cm^2', 'description': 'Hole current density (y)'},

            # Voltages and biases
            'vds': {'symbol': r'V_{DS}', 'unit': 'V', 'description': 'Drain-source voltage'},
            'vgs': {'symbol': r'V_{GS}', 'unit': 'V', 'description': 'Gate-source voltage'},
            'vbs': {'symbol': r'V_{BS}', 'unit': 'V', 'description': 'Bulk-source voltage'},
            'v_src': {'symbol': r'V_{SRC}', 'unit': 'V', 'description': 'Source voltage'},
            'v_drn': {'symbol': r'V_{DRN}', 'unit': 'V', 'description': 'Drain voltage'},
            'v_gat': {'symbol': r'V_{GAT}', 'unit': 'V', 'description': 'Gate voltage'},
            'v_blk': {'symbol': r'V_{BLK}', 'unit': 'V', 'description': 'Bulk voltage'},

            # Device currents
            'ids': {'symbol': r'I_{DS}', 'unit': 'A', 'description': 'Drain-source current'},
            'igs': {'symbol': r'I_{GS}', 'unit': 'A', 'description': 'Gate-source current'},
            'i_src': {'symbol': r'I_{SRC}', 'unit': 'A', 'description': 'Source current'},
            'i_drn': {'symbol': r'I_{DRN}', 'unit': 'A', 'description': 'Drain current'},
            'i_gat': {'symbol': r'I_{GAT}', 'unit': 'A', 'description': 'Gate current'},

            # Material properties
            'doping': {'symbol': r'N_D - N_A', 'unit': r'cm^{-3}', 'description': 'Net doping'},
            'nd': {'symbol': r'N_D', 'unit': r'cm^{-3}', 'description': 'Donor concentration'},
            'na': {'symbol': r'N_A', 'unit': r'cm^{-3}', 'description': 'Acceptor concentration'},
            'mobility': {'symbol': r'\mu', 'unit': r'cm^2/V/s', 'description': 'Mobility'},
            'mu_n': {'symbol': r'\mu_n', 'unit': r'cm^2/V/s', 'description': 'Electron mobility'},
            'mu_p': {'symbol': r'\mu_p', 'unit': r'cm^2/V/s', 'description': 'Hole mobility'},

            # Time and frequency
            'time': {'symbol': r't', 'unit': 's', 'description': 'Time'},
            'freq': {'symbol': r'f', 'unit': 'Hz', 'description': 'Frequency'},
            'omega': {'symbol': r'\omega', 'unit': 'rad/s', 'description': 'Angular frequency'},

            # Energy and bandstructure
            'energy': {'symbol': r'E', 'unit': 'eV', 'description': 'Energy'},
            'ec': {'symbol': r'E_c', 'unit': 'eV', 'description': 'Conduction band edge'},
            'ev': {'symbol': r'E_v', 'unit': 'eV', 'description': 'Valence band edge'},
            'ef': {'symbol': r'E_F', 'unit': 'eV', 'description': 'Fermi level'},
            'eg': {'symbol': r'E_g', 'unit': 'eV', 'description': 'Band gap'},

            # Convergence and iterations
            'iteration': {'symbol': r'n', 'unit': '', 'description': 'Iteration number'},
            'residual': {'symbol': r'||R||', 'unit': '', 'description': 'Residual norm'},
            'error': {'symbol': r'\epsilon', 'unit': '', 'description': 'Error'},

            # Temperature
            'temperature': {'symbol': r'T', 'unit': 'K', 'description': 'Temperature'},
            'temp': {'symbol': r'T', 'unit': 'K', 'description': 'Temperature'},
        }

        # Special cases and aliases
        self._aliases = {
            # Common aliases
            'pot': 'phi',
            'potential': 'phi',
            'electron_density': 'ndens',
            'hole_density': 'pdens',
            'n_density': 'ndens',
            'p_density': 'pdens',

            # Coordinate aliases
            'position': 'x',
            'coord': 'x',

            # Current aliases
            'current_density': 'jn',
            'j_electron': 'jn',
            'j_hole': 'jp',
        }

    def get_label(self, var_name: str, unit_override: Optional[str] = None) -> str:
        """
        Get LaTeX label for a variable

        Parameters:
        -----------
        var_name : str
            Variable name (e.g., 'x', 'ndens', 'pot')
        unit_override : str, optional
            Override default unit

        Returns:
        --------
        str
            LaTeX-formatted label with unit
        """
        # Handle aliases
        lookup_name = self._aliases.get(var_name.lower(), var_name.lower())

        # Check coordinates first
        if lookup_name in self._coordinates:
            info = self._coordinates[lookup_name]
        elif lookup_name in self._quantities:
            info = self._quantities[lookup_name]
        else:
            # Fallback for unknown variables
            return f"${var_name}$"

        # Use unit override if provided
        unit = unit_override if unit_override is not None else info['unit']

        # Format label for LaTeX rendering (proper SciencePlots format)
        symbol = info['symbol']
        if unit:
            # Use \mathrm{} inside math mode for upright text units (more reliable than \textup{})
            return r'$' + symbol + r'$ ($\mathrm{' + unit + r'}$)'
        else:
            return r'$' + symbol + r'$'

    def get_labels(self, x_var: str, y_var: str,
                   x_unit: Optional[str] = None, y_unit: Optional[str] = None) -> Tuple[str, str]:
        """
        Get xlabel and ylabel pair

        Parameters:
        -----------
        x_var, y_var : str
            Variable names for x and y axes
        x_unit, y_unit : str, optional
            Unit overrides

        Returns:
        --------
        tuple
            (xlabel, ylabel) tuple with LaTeX formatting
        """
        xlabel = self.get_label(x_var, x_unit)
        ylabel = self.get_label(y_var, y_unit)
        return xlabel, ylabel

    def get_title(self, var_name: str) -> str:
        """
        Get descriptive title for a variable

        Parameters:
        -----------
        var_name : str
            Variable name

        Returns:
        --------
        str
            Descriptive title
        """
        lookup_name = self._aliases.get(var_name.lower(), var_name.lower())

        if lookup_name in self._coordinates:
            return self._coordinates[lookup_name]['description']
        elif lookup_name in self._quantities:
            return self._quantities[lookup_name]['description']
        else:
            return var_name.replace('_', ' ').title()

    def add_custom_label(self, var_name: str, symbol: str, unit: str, description: str):
        """
        Add custom label to dictionary

        Parameters:
        -----------
        var_name : str
            Variable name key
        symbol : str
            LaTeX symbol (without $)
        unit : str
            Physical unit
        description : str
            Descriptive name
        """
        self._quantities[var_name.lower()] = {
            'symbol': symbol,
            'unit': unit,
            'description': description
        }

    def list_available(self) -> Dict[str, str]:
        """
        List all available variables and their symbols

        Returns:
        --------
        dict
            Dictionary mapping variable names to symbols
        """
        result = {}

        # Add coordinates
        for name, info in self._coordinates.items():
            result[name] = info['symbol']

        # Add quantities
        for name, info in self._quantities.items():
            result[name] = info['symbol']

        return result
