"""
Planet system configurations for the impact-reinflation model.

This module defines various planet-star systems that can be used across
different simulation scripts. Each system configuration includes:
- Planet properties (name, mass, radius, semi-major axis)
- Stellar properties (luminosity, radius, saturation timescale)
- Additional parameters (escape efficiency, etc.)

Usage:
    from planet_systems import TEST_PLANET, get_planet_system
    
    # Use predefined test planet
    star, planet = get_planet_system('TEST_PLANET')
    
    # Or use a real planet
    star, planet = get_planet_system('GJ3929b')
"""

from impactmodel import Planet, Star
from typing import Tuple, Dict, Any


# ==============================
# PLANET SYSTEM CONFIGURATIONS
# ==============================

TEST_PLANET = {
    'name': 'TestPlanet',
    'mass': 1.0,  # Earth masses
    'radius': 1.0,  # Earth radii
    'distance': 0.05,  # AU
    'star': {
        'luminosity': 0.01,  # L_sun
        'radius': 0.25,  # R_sun
        't_sat': 1.3e9  # years
    },
    'escape_efficiency': 1e-2
}

EARTH_LIKE = {
    'name': 'Earth-like',
    'mass': 1.0,  # Earth masses
    'radius': 1.0,  # Earth radii
    'distance': 0.1,  # AU
    'star': {
        'luminosity': 0.01,  # L_sun
        'radius': 0.25,  # R_sun
        't_sat': 1.3e9  # years
    },
    'escape_efficiency': 1e-2
}

# Known exoplanet systems from observations
GJ3929B = {
    'name': 'GJ3929b',
    'mass': 1.75,  # Earth masses
    'radius': 1.09,  # Earth radii
    'distance': 0.0252,  # AU
    'star': {
        'radius': 0.32,  # R_sun
        'luminosity': 10**(-1.96),  # L_sun
        't_sat': 1.1e9  # years
    },
    'escape_efficiency': 1e-2
}

LTT1445AC = {
    'name': 'LTT1445Ac',
    'mass': 1.37,  # Earth masses
    'radius': 1.07,  # Earth radii
    'distance': 0.02659,  # AU
    'star': {
        'radius': 0.271,  # R_sun
        'luminosity': 10**(-2.04),  # L_sun
        't_sat': 1.9e9  # years
    },
    'escape_efficiency': 1e-2
}

LTT1445AB = {
    'name': 'LTT1445Ab',
    'mass': 2.73,  # Earth masses
    'radius': 1.34,  # Earth radii
    'distance': 0.03810,  # AU
    'star': {
        'radius': 0.271,  # R_sun
        'luminosity': 10**(-2.04),  # L_sun
        't_sat': 1.9e9  # years
    },
    'escape_efficiency': 1e-2
}

# Dictionary of all available systems
PLANET_SYSTEMS = {
    'TEST_PLANET': TEST_PLANET,
    'EARTH_LIKE': EARTH_LIKE,
    'GJ3929b': GJ3929B,
    'LTT1445Ac': LTT1445AC,
    'LTT1445Ab': LTT1445AB,
}


# ==============================
# HELPER FUNCTIONS
# ==============================

def get_planet_system(system_name: str) -> Tuple[Star, Planet]:
    """
    Create Star and Planet objects from a system configuration.
    
    Args:
        system_name: Name of the planet system (e.g., 'TEST_PLANET', 'GJ3929b')
    
    Returns:
        Tuple of (Star, Planet) objects
    
    Raises:
        KeyError: If system_name is not found in PLANET_SYSTEMS
    """
    if system_name not in PLANET_SYSTEMS:
        available = ', '.join(PLANET_SYSTEMS.keys())
        raise KeyError(
            f"Unknown planet system '{system_name}'. "
            f"Available systems: {available}"
        )
    
    config = PLANET_SYSTEMS[system_name]
    star_config = config['star']
    
    # Create star
    star = Star(
        L_bol=star_config['luminosity'],
        Rs=star_config['radius'],
        t_sat=star_config['t_sat']
    )
    
    # Create planet
    planet = Planet(
        name=config['name'],
        mass=config['mass'],
        radius=config['radius'],
        star=star,
        a=config['distance'],
        escape_efficiency=config['escape_efficiency']
    )
    
    return star, planet


def get_planet_config(system_name: str) -> Dict[str, Any]:
    """
    Get the configuration dictionary for a planet system.
    
    Args:
        system_name: Name of the planet system
    
    Returns:
        Configuration dictionary
    """
    if system_name not in PLANET_SYSTEMS:
        available = ', '.join(PLANET_SYSTEMS.keys())
        raise KeyError(
            f"Unknown planet system '{system_name}'. "
            f"Available systems: {available}"
        )
    
    return PLANET_SYSTEMS[system_name].copy()


def list_available_systems() -> list:
    """
    List all available planet system names.
    
    Returns:
        List of system names
    """
    return list(PLANET_SYSTEMS.keys())


if __name__ == "__main__":
    # Print available systems
    print("Available planet systems:")
    for name, config in PLANET_SYSTEMS.items():
        print(f"  {name}: {config['name']} "
              f"({config['mass']} M_earth, {config['radius']} R_earth, "
              f"a={config['distance']} AU)")
    
    # Test creating a system
    print("\nTesting system creation:")
    star, planet = get_planet_system('TEST_PLANET')
    print(f"  Star: L={star.L_bol:.4f} L_sun, Rs={star.Rs:.3f} R_sun")
    print(f"  Planet: {planet.name}, M={planet.mass_earth:.2f} M_earth, "
          f"R={planet.radius_earth:.2f} R_earth")

