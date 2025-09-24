"""
Label dictionary for Cybear plots - clean LaTeX formatting
"""

# Physics field labels with proper LaTeX formatting
FIELD_LABELS = {
    # Potential and fields
    'pot': (r'$\phi$', 'V'),
    'efield': (r'$E$', 'V/cm'),
    'Ex': (r'$E_x$', 'V/cm'),
    'Ey': (r'$E_y$', 'V/cm'),
    'Ez': (r'$E_z$', 'V/cm'),

    # Carrier densities
    'ndens': (r'$n$', r'cm$^{-3}$'),
    'pdens': (r'$p$', r'cm$^{-3}$'),

    # Current densities
    'current': (r'$J$', r'A/cm$^2$'),
    'jn': (r'$J_n$', r'A/cm$^2$'),
    'jp': (r'$J_p$', r'A/cm$^2$'),

    # Material properties
    'mobility': (r'$\mu$', r'cm$^2$/V·s'),
    'temperature': (r'$T$', 'K'),
    'recomb': (r'$R$', r'cm$^{-3}$s$^{-1}$'),

    # Coordinates
    'x': (r'$x$', 'nm'),
    'y': (r'$y$', 'nm'),
    'z': (r'$z$', 'nm'),
}

def get_label(field_name):
    """Get formatted label and unit for a field"""
    return FIELD_LABELS.get(field_name, (field_name, ''))

def format_axis_label(field_name):
    """Get full axis label string"""
    label, unit = get_label(field_name)
    if unit:
        return f"{label} ({unit})"
    return label