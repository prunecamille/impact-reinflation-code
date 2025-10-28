"""
Monte Carlo evolution script for the impact-reinflation model.

This script runs multiple simulations over a parameter space defined by ranges
of outgassing rates and impact rates, allowing exploration of atmospheric 
evolution under different conditions.

User-configurable parameters:
- Parameter ranges (outgassing rate, impact rate, planet properties)
- Number of samples for each parameter
- Stellar properties
- Output parameters

The script saves results to a numpy npz file and a CSV file for each simulation.
"""

import numpy as np
import sys
import os
from impactmodel import Planet, Star
from impactmodel.constants import CURRENT_EARTH_OUTGASSING_RATE

# ==============================
# USER-CONFIGURABLE PARAMETERS
# ==============================

# Stellar properties (fixed for all runs)
L_bol = 0.01  # Bolometric luminosity in solar units
Rs = 0.25  # Stellar radius in solar radii
t_sat = 1.3e9  # Saturation time in years

# Planet properties (fixed for all runs)
planet_name = "TestPlanet"
planet_mass = 1.0  # Earth masses
planet_radius = 1.0  # Earth radii
semi_major_axis = 0.1  # AU

# Parameter ranges for Monte Carlo sampling
# These define the grid of parameters to explore
outgassing_rate_min = -1  # log10 relative to Earth
outgassing_rate_max = 1  # log10 relative to Earth
outgassing_rate_num = 5  # Number of samples

impact_rate_min = 1e-11  # yr^-1
impact_rate_max = 1e-8  # yr^-1
impact_rate_num = 5  # Number of samples

# Additional parameters
escape_efficiency = 1e-2  # Escape efficiency parameter

# Evolution parameters
age = 5e9  # years
dt = 1e6  # years per time step
method = 'euler'  # Integration method

# Output parameters
output_dir = 'montecarlo_results'
output_prefix = 'montecarlo'

# ==============================
# RUN MONTE CARLO SIMULATION
# ==============================

def main():
    """Run Monte Carlo simulations."""
    
    print("=" * 60)
    print("Impact-Reinflation Model: Monte Carlo Simulation")
    print("=" * 60)
    
    # Create parameter grids
    outgassing_rates = np.logspace(
        outgassing_rate_min, 
        outgassing_rate_max, 
        outgassing_rate_num
    ) * CURRENT_EARTH_OUTGASSING_RATE
    
    impact_rates = np.logspace(
        np.log10(impact_rate_min),
        np.log10(impact_rate_max),
        impact_rate_num
    )
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Print simulation parameters
    print(f"\nStellar properties:")
    print(f"  L_bol: {L_bol} L_sun")
    print(f"  Rs: {Rs} R_sun")
    print(f"  t_sat: {t_sat:.2e} years")
    
    print(f"\nPlanet properties:")
    print(f"  Mass: {planet_mass} M_earth")
    print(f"  Radius: {planet_radius} R_earth")
    print(f"  Semi-major axis: {semi_major_axis} AU")
    
    print(f"\nParameter ranges:")
    print(f"  Outgassing rates: {outgassing_rate_min} to {outgassing_rate_max} (log10)")
    print(f"    Samples: {outgassing_rate_num}")
    print(f"  Impact rates: {impact_rate_min:.2e} to {impact_rate_max:.2e} yr^-1")
    print(f"    Samples: {impact_rate_num}")
    
    total_runs = len(outgassing_rates) * len(impact_rates)
    print(f"\nTotal runs: {total_runs}")
    
    # Initialize arrays to store results
    all_results = []
    
    # Create star
    star = Star(L_bol=L_bol, Rs=Rs, t_sat=t_sat)
    
    # Run simulations
    run_count = 0
    for i, outgassing_rate in enumerate(outgassing_rates):
        for j, impact_rate in enumerate(impact_rates):
            run_count += 1
            print(f"\n[{run_count}/{total_runs}] "
                  f"Outgassing: {outgassing_rate/CURRENT_EARTH_OUTGASSING_RATE:.2e}x Earth, "
                  f"Impact rate: {impact_rate:.2e} yr^-1")
            
            # Create planet
            planet = Planet(
                name=f"{planet_name}_{i}_{j}",
                mass=planet_mass,
                radius=planet_radius,
                star=star,
                a=semi_major_axis,
                outgassing_rate=outgassing_rate,
                escape_efficiency=escape_efficiency
            )
            
            # Add atmosphere (starting with no atmosphere)
            planet.add_atmosphere(pressure0=0, condensed_volatile_thickness0=0)
            critical_pressure = planet.atmosphere.critical_pressure
            
            # Add impactor
            impactor_size = planet.add_default_impactor(impact_rate)
            
            # Evolve the system
            results = planet.evolve_system(age=age, dt=dt, method=method)
            
            # Store results
            all_results.append({
                'outgassing_rate': outgassing_rate,
                'impact_rate': impact_rate,
                'impactor_size': impactor_size,
                'num_impacts': results['num_impacts'],
                'inflation_percentage': results['inflation_percentage'],
                'final_pressure': results['pressures'][-1],
                'final_atm_mass': results['atm_masses'][-1],
                'final_vol_mass': results['vol_masses'][-1],
                'critical_pressure': critical_pressure * 1e-5,  # Convert to bar
                'times': results['times'],
                'atm_masses': results['atm_masses'],
                'vol_masses': results['vol_masses'],
                'pressures': results['pressures'],
                'states': results['states'],
                'impact_times': results['impact_times']
            })
            
            print(f"  -> Impacts: {results['num_impacts']}, "
                  f"Inflation: {results['inflation_percentage']:.1f}%")
    
    # Save results
    print(f"\n{'=' * 60}")
    print("Saving results...")
    
    output_file = f"{output_dir}/{output_prefix}.npz"
    np.savez(output_file, *[r for r in all_results])
    
    # Also save a summary CSV
    csv_file = f"{output_dir}/{output_prefix}_summary.csv"
    import csv
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'outgassing_rate', 'impact_rate', 'impactor_size', 'num_impacts',
            'inflation_percentage', 'final_pressure', 'final_atm_mass', 
            'final_vol_mass', 'critical_pressure'
        ])
        for r in all_results:
            writer.writerow([
                r['outgassing_rate'], r['impact_rate'], r['impactor_size'],
                r['num_impacts'], r['inflation_percentage'], r['final_pressure'],
                r['final_atm_mass'], r['final_vol_mass'], r['critical_pressure']
            ])
    
    print(f"Results saved to:")
    print(f"  {output_file}")
    print(f"  {csv_file}")
    
    print("\n" + "=" * 60)
    print("Monte Carlo simulation complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()
