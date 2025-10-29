"""
Physical constants used throughout the impact-reinflation model.

This module contains fundamental physical constants including astronomical 
constants (Earth parameters, stellar parameters), material properties (CO2 
sublimation properties), and unit conversions.
"""

import numpy as np

# Gravitational constant
G = 6.67e-11  # m^3 kg^-1 s^-2

# Earth parameters
M_EARTH = 5.972e24  # kg
R_EARTH = 6.371e6  # m
M_EARTH_KG = M_EARTH  # Alias for clarity
R_EARTH_M = R_EARTH  # Alias for clarity
g_EARTH = G * M_EARTH / R_EARTH**2  # m/s^2 (surface gravity)

# Stellar parameters (Sun)
M_SUN = 1.989e30  # kg
R_SUN = 6.957e8  # m
L_SUN = 3.828e26  # W (bolometric luminosity)

# CO2 properties
LATENT_HEAT_CO2 = 571e3  # J/kg (latent heat of sublimation)
CP_CO2 = 820  # J/(kg·K) (specific heat capacity)
T_SUBLIMATION_CO2 = 195  # K (sublimation temperature)
KAPPA_CO2 = 1.6e-4  # m²/kg (CO2 opacity, see Wordsworth et al. (2010b))
CO2_TRIPLE_POINT_PRESSURE = 5.185e5  # Pa (CO2 triple point pressure, see Wordsworth et al. (2010b))
CO2_TRIPLE_POINT_TEMP = 216.592  # K (CO2 triple point temperature, see Wordsworth et al. (2010b))

# Unit conversions
AU_TO_M = 1.496e11  # meters per astronomical unit
SECONDS_PER_YEAR = 365.25 * 24 * 60 * 60  # seconds per year

# Earth reference values
CURRENT_EARTH_OUTGASSING_RATE = 6.5e-4  # kg/m²/yr (middle ground value from Wallman & Aloisi (2012))
EARTH_INSOLATION = 1361  # W/m²

# Physical constants
SIGMA_BOLTZMANN = 5.670374419e-8  # W/m²/K⁴ (Stefan-Boltzmann constant)
KB = 1.380649e-23  # J/K (Boltzmann constant)
K_BOLTZMANN = KB  # Alias for consistency
