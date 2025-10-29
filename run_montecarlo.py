"""
Monte Carlo evolution script for the impact-reinflation model.

This script runs multiple simulations over a parameter space defined by ranges
of outgassing rates and impact rates, using random Monte Carlo sampling to 
explore atmospheric evolution under different conditions.

User-configurable parameters:
- Parameter ranges (outgassing rate, impact rate, planet properties)
- Number of random Monte Carlo samples
- Stellar properties
- Output parameters

The script saves results to a numpy npz file and a CSV file with all simulation results.
"""

import numpy as np
import sys
import os
from impactmodel.constants import CURRENT_EARTH_OUTGASSING_RATE, M_EARTH, R_EARTH, AU_TO_M
from planet_systems import get_planet_system, get_planet_config

# ==============================
# USER-CONFIGURABLE PARAMETERS
# ==============================

# Planet system to use (see planet_systems.py for available systems)
PLANET_SYSTEM = 'LTT1445Ac'  # Options: 'TEST_PLANET', 'EARTH_LIKE', 'GJ3929b', etc.

# Parameter ranges for Monte Carlo sampling
# These define the ranges to sample from randomly (matching parameter_space_plot.py)
outgassing_rate_min = -1.5  # log10 relative to Earth
outgassing_rate_max = 0.5   # log10 relative to Earth

impact_rate_min = -12  # log10 impacts/year
impact_rate_max = -6   # log10 impacts/year

# Number of random Monte Carlo samples
n_samples = 1000  # Number of random samples to draw


# Evolution parameters (matching parameter_space_plot.py)
age = 12e9  # years (12 billion years)
dt = 1e7    # years per time step
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
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get planet system - we need both the objects and the config
    star, _ = get_planet_system(PLANET_SYSTEM)  # Star object for reuse
    planet_config = get_planet_config(PLANET_SYSTEM)  # Config dict with Earth units
    
    # Print simulation parameters
    print(f"\nStellar properties:")
    print(f"  L_bol: {star.L_bol:.4f} L_sun")
    print(f"  Rs: {star.Rs:.3f} R_sun")
    print(f"  t_sat: {star.t_sat:.2e} years")
    
    print(f"\nPlanet system: {PLANET_SYSTEM}")
    print(f"  Planet: {planet_config['name']}")
    print(f"  Mass: {planet_config['mass']:.2f} M_earth")
    print(f"  Radius: {planet_config['radius']:.2f} R_earth")
    print(f"  Semi-major axis: {planet_config['distance']:.4f} AU")
    
    print(f"\nParameter ranges:")
    print(f"  Outgassing rates relative to Earth: {10**outgassing_rate_min} to {10**outgassing_rate_max} (x Earth)")
    print(f"  Impact rates in impacts/Gyr: {10**impact_rate_min*1e9} to {10**impact_rate_max*1e9} (impacts/Gyr)")
    print(f"\nNumber of Monte Carlo samples: {n_samples}")
    
    # Initialize arrays to store results
    all_results = []
    
    # Generate random Monte Carlo samples
    np.random.seed(42)  # This ensures reproducibility, remove to have actual random samples
    # randomly sample in log space
    outgassing_logs = np.random.uniform(outgassing_rate_min, outgassing_rate_max, n_samples)
    impact_logs = np.random.uniform(impact_rate_min, impact_rate_max, n_samples)
    outgassing_factors = 10**outgassing_logs
    impact_rates = 10**impact_logs
    
    print(f"\nRunning {n_samples} Monte Carlo samples...")
    
    # Run simulations
    for i in range(n_samples):
        if (i+1) % 100 == 0:
            print(f"  Progress: {i+1}/{n_samples} ({100*(i+1)/n_samples:.1f}%)")
        
        # Get random parameters
        outgassing_factor = outgassing_factors[i]
        impact_rate = impact_rates[i]
        
        # Calculate actual outgassing rate
        outgassing_rate = outgassing_factor * CURRENT_EARTH_OUTGASSING_RATE # convert to kg/m²/yr
        
        # Create planet (use config values in Earth units, reuse star)
        from impactmodel import Planet
        planet = Planet(
            name=f"{planet_config['name']}_{i}",
            mass=planet_config['mass'],  # Earth masses
            radius=planet_config['radius'],  # Earth radii
            star=star,
            a=planet_config['distance'],  # AU
            outgassing_rate=outgassing_rate,
            escape_efficiency=planet_config['escape_efficiency']
        )
        
        # Add atmosphere (starting with no atmosphere)
        planet.add_atmosphere(pressure0=0, condensed_volatile_thickness0=0)
        critical_pressure = planet.atmosphere.critical_pressure
        
        # Add impactor
        impactor_size = planet.add_default_impactor(impact_rate)
        
        # Evolve the system
        results = planet.evolve_system(age=age, dt=dt, method=method)
        
        # Store results (maybe not the best way to do this, but it works for now)
        all_results.append({
            'outgassing_rate': outgassing_rate,
            'outgassing_factor': outgassing_factor,
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
    
    # Save results
    print(f"\n{'=' * 60}")
    print("Saving results...")
    
    # Prepare arrays for npz save (matching parameter_space_plot.py format)
    # Note: 'outgassing_rates' stores linear factors (x Earth), matching parameter_space_plot.py
    outgassing_factors_array = np.array([r['outgassing_factor'] for r in all_results])
    impact_rates_array = np.array([r['impact_rate'] for r in all_results])
    impactor_sizes_array = np.array([r['impactor_size'] for r in all_results])
    inflation_percentages_array = np.array([r['inflation_percentage'] for r in all_results])
    critical_pressures_array = np.array([r['critical_pressure'] for r in all_results])
    
    output_file = f"{output_dir}/{output_prefix}.npz"
    # Match parameter_space_plot.py format: 'outgassing_rates' stores linear factors (x Earth)
    np.savez(
        output_file,
        outgassing_rates=outgassing_factors_array,  # Linear factors (x Earth), not actual rates
        impact_rates=impact_rates_array,
        impactor_sizes=impactor_sizes_array,
        inflation_percentages=inflation_percentages_array,
        critical_pressures=critical_pressures_array
    )
    
    # Also save a summary CSV (matching parameter_space_plot.py format)
    csv_file = f"{output_dir}/{output_prefix}_summary.csv"
    import csv
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'outgassing_rate_factor', 'impact_rate', 'inflation_percentage',
            'impactor_size', 'critical_pressure'
        ])
        for r in all_results:
            writer.writerow([
                r['outgassing_factor'],  # Store linear factor (x Earth), not actual rate
                r['impact_rate'],
                r['inflation_percentage'],
                r['impactor_size'],
                r['critical_pressure']
            ])
    
    print(f"Results saved to:")
    print(f"  {output_file}")
    print(f"  {csv_file}")
    
    print("\n" + "=" * 60)
    print("Monte Carlo simulation complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()
