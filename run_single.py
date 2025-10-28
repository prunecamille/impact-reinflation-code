"""
Single evolution run script for the impact-reinflation model.

This script demonstrates how to set up and run a single evolution simulation
for a planet with specified stellar properties, outgassing rate, and impact rate.

User-configurable parameters:
- Stellar properties (L_bol, Rs, t_sat)
- Planet properties (mass, radius, semi-major axis)
- Outgassing rate
- Impact rate
- Escape efficiency
- Evolution time and time step

The script saves results to a numpy npz file for later analysis and plotting.
"""

import numpy as np
import sys
from impactmodel import Planet, Star

# ==============================
# USER-CONFIGURABLE PARAMETERS
# ==============================

# Stellar properties
L_bol = 0.01  # Bolometric luminosity in solar units
Rs = 0.25  # Stellar radius in solar radii
t_sat = 1.3e9  # Saturation time in years

# Planet properties
planet_name = "Earth-like"
planet_mass = 1.0  # Earth masses
planet_radius = 1.0  # Earth radii
semi_major_axis = 0.1  # AU

# Atmospheric evolution parameters
outgassing_rate = 6.5e-4  # kg/m²/yr (default is Earth-like)
escape_efficiency = 1e-2  # Escape efficiency parameter

# Impact parameters
impact_rate = 1e-9  # Impacts per year

# Evolution parameters
age = 5e9  # years
dt = 1e6  # years per time step
method = 'euler'  # Integration method: 'euler' or 'rk4'

# Output parameters
output_file = 'results_single.npz'

# ==============================
# RUN SIMULATION
# ==============================

def main():
    """Run a single evolution simulation."""
    
    print("=" * 60)
    print("Impact-Reinflation Model: Single Run")
    print("=" * 60)
    
    # Create star
    print(f"\nCreating star with:")
    print(f"  Bolometric luminosity: {L_bol} L_sun")
    print(f"  Stellar radius: {Rs} R_sun")
    print(f"  Saturation time: {t_sat:.2e} years")
    star = Star(L_bol=L_bol, Rs=Rs, t_sat=t_sat)
    
    # Create planet
    print(f"\nCreating planet '{planet_name}' with:")
    print(f"  Mass: {planet_mass} M_earth")
    print(f"  Radius: {planet_radius} R_earth")
    print(f"  Semi-major axis: {semi_major_axis} AU")
    print(f"  Outgassing rate: {outgassing_rate:.2e} kg/m²/yr")
    print(f"  Escape efficiency: {escape_efficiency}")
    planet = Planet(
        name=planet_name,
        mass=planet_mass,
        radius=planet_radius,
        star=star,
        a=semi_major_axis,
        outgassing_rate=outgassing_rate,
        escape_efficiency=escape_efficiency
    )
    
    # Add atmosphere (starting with no atmosphere)
    print("\nAdding atmosphere...")
    planet.add_atmosphere(pressure0=0, condensed_volatile_thickness0=0)
    critical_pressure = planet.atmosphere.critical_pressure
    print(f"  Critical pressure: {critical_pressure*1e-5:.4f} bar")
    
    # Add impactor
    print(f"\nAdding impactor with impact rate: {impact_rate:.2e} yr^-1")
    impactor_size = planet.add_default_impactor(impact_rate)
    print(f"  Impactor size: {impactor_size:.2f} km")
    
    # Evolve the system
    print(f"\nEvolving system for {age:.2e} years...")
    print(f"  Time step: {dt:.2e} years")
    print(f"  Method: {method}")
    
    results = planet.evolve_system(age=age, dt=dt, method=method)
    
    # Print summary
    print("\n" + "=" * 60)
    print("EVOLUTION SUMMARY")
    print("=" * 60)
    print(f"Number of impacts: {results['num_impacts']}")
    print(f"Time in inflated state: {results['inflation_percentage']:.2f}%")
    print(f"Final pressure: {results['pressures'][-1]:.4e} bar")
    print(f"Final atmospheric mass: {results['atm_masses'][-1]:.2e} kg")
    print(f"Final condensed mass: {results['vol_masses'][-1]:.2e} kg")
    
    # Save results
    print(f"\nSaving results to: {output_file}")
    np.savez(
        output_file,
        times=results['times'],
        atm_masses=results['atm_masses'],
        vol_masses=results['vol_masses'],
        pressures=results['pressures'],
        states=results['states'],
        impact_times=results['impact_times'],
        inflation_percentage=results['inflation_percentage'],
        num_impacts=results['num_impacts'],
        # Store simulation parameters
        planet_name=planet_name,
        planet_mass=planet_mass,
        planet_radius=planet_radius,
        semi_major_axis=semi_major_axis,
        outgassing_rate=outgassing_rate,
        escape_efficiency=escape_efficiency,
        impact_rate=impact_rate,
        age=age,
        dt=dt,
        critical_pressure=critical_pressure*1e-5,  # Convert to bar
        L_bol=L_bol,
        Rs=Rs,
        t_sat=t_sat
    )
    
    print("\n" + "=" * 60)
    print("Simulation complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()
