"""
Demo script to run three simulations of Earth-like planets around an M0-type star
at different orbital distances: 0.01, 0.05, and 0.1 AU.
"""

import numpy as np
import sys
import os

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEMO_DIR = SCRIPT_DIR  # The demo directory

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(SCRIPT_DIR, '..')))

from impactmodel import Planet, Star, constants
constants.CURRENT_EARTH_OUTGASSING_RATE = 6.5e-4  # kg/m²/yr
from plot_results import plot_single_evolution

def main():
    """Run three simulations and generate plots for different orbital distances."""
    
    print("=" * 80)
    print("Impact-Reinflation Model: M0-Type Star Demo")
    print("=" * 80)
    
    # Stellar properties for M0-type star
    L_bol = 0.072  # Bolometric luminosity in solar units
    Rs = 0.62  # Stellar radius in solar radii
    t_sat = 1e9  # Saturation time in years
    
    # Planetary properties (Earth-like)
    planet_name = "Earth-like"
    planet_mass = 1.0  # Earth masses
    planet_radius = 1.0  # Earth radii
    
    # Orbital distances to test
    semi_major_axes = [0.01, 0.05, 0.1]  # AU
    
    # Atmospheric evolution parameters (Earth-like)
    outgassing_rate = constants.CURRENT_EARTH_OUTGASSING_RATE  # kg/m²/yr
    escape_efficiency = 1e-2  # Escape efficiency parameter
    
    # Impact parameters
    impact_rate = 1e-9  # Impacts per year
    
    # Evolution parameters
    age = 5e9  # years
    dt = 1e6  # years per time step
    method = 'euler'
    
    # Create star
    print(f"\nCreating M0-type star with:")
    print(f"  Bolometric luminosity: {L_bol} L_sun")
    print(f"  Stellar radius: {Rs} R_sun")
    print(f"  Saturation time: {t_sat:.2e} years")
    star = Star(L_bol=L_bol, Rs=Rs, t_sat=t_sat)
    
    results_files = []
    
    # Run simulations for each orbital distance
    for a in semi_major_axes:
        print("\n" + "=" * 80)
        print(f"Running simulation for semi-major axis: {a} AU")
        print("=" * 80)
        
        # Create planet
        print(f"\nCreating planet '{planet_name}' with:")
        print(f"  Mass: {planet_mass} M_earth")
        print(f"  Radius: {planet_radius} R_earth")
        print(f"  Semi-major axis: {a} AU")
        print(f"  Outgassing rate: {outgassing_rate:.2e} kg/m²/yr")
        print(f"  Escape efficiency: {escape_efficiency}")
        
        planet = Planet(
            name=planet_name,
            mass=planet_mass,
            radius=planet_radius,
            star=star,
            a=a,
            outgassing_rate=outgassing_rate,
            escape_efficiency=escape_efficiency
        )
        
        # Add atmosphere (starting with no atmosphere)
        print("\nAdding atmosphere...")
        planet.add_atmosphere(pressure0=0, condensed_volatile_thickness0=0)
        critical_pressure = planet.atmosphere.critical_pressure
        print(f"  Critical pressure: {critical_pressure*1e-5:.4f} bar")
        print(f"  Critical mass: {planet.atmosphere.critical_mass:.2e} kg")
        
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
        print("\n" + "=" * 80)
        print("EVOLUTION SUMMARY")
        print("=" * 80)
        print(f"Number of impacts: {results['num_impacts']}")
        print(f"Time in inflated state: {results['inflation_percentage']:.2f}%")
        print(f"Final pressure: {results['pressures'][-1]:.4e} bar")
        print(f"Final atmospheric mass: {results['atm_masses'][-1]:.2e} kg")
        print(f"Final condensed mass: {results['vol_masses'][-1]:.2e} kg")
        
        # Save results
        # Format filename properly for small decimal AU values
        a_str = str(a).replace('.', 'p')  # 0.01 -> "0p01", 0.05 -> "0p05", 0.1 -> "0p1"
        output_file = os.path.join(DEMO_DIR, f'3earths/earth_M0_{a_str}AU.npz')
        results_files.append(output_file)
        
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
            semi_major_axis=a,
            outgassing_rate=outgassing_rate,
            escape_efficiency=escape_efficiency,
            impact_rate=impact_rate,
            age=age,
            dt=dt,
            critical_pressure=critical_pressure*1e-5,  # Convert to bar
            impactor_diameter=impactor_size,  # km
            L_bol=L_bol,
            Rs=Rs,
            t_sat=t_sat
        )
        
        # Generate plot
        plot_file = output_file.replace('.npz', '.png')
        print(f"Generating plot: {plot_file}")
        plot_single_evolution(output_file, save_path=plot_file)
    
    print("\n" + "=" * 80)
    print("All simulations complete!")
    print("=" * 80)
    print("\nResults saved:")
    for rf in results_files:
        print(f"  - {rf}")
        print(f"  - {rf.replace('.npz', '.png')}")
    
    # Create a comparison plot showing all three scenarios
    print("\nGenerating comparison plot...")
    create_comparison_plot(results_files)
    
    print("\n" + "=" * 80)
    print("Demo complete!")
    print("=" * 80)


def create_comparison_plot(results_files):
    """Create a comparison plot showing all three scenarios."""
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    
    colors = {
        'inflated': '#06c2ac',
        'collapsed': '#742E1E',
        'pressure': '#380835',
        'impacts': '#e91e63',  # Reddish-pink for impacts and critical pressure
        'critical': '#e91e63'
    }
    
    # Extract semi-major axes and get impactor info from filenames
    for idx, data_file in enumerate(results_files):
        ax = axes[idx]
        
        # Load data
        data = np.load(data_file)
        times = data['times'] / 1e9  # Convert to Gyr
        pressures = data['pressures']
        states = data['states']
        impact_times = data['impact_times'] / 1e9  # Convert to Gyr
        critical_pressure = data['critical_pressure']  # Already in bar
        semi_major_axis = data['semi_major_axis']
        
        # Get impactor diameter from saved data
        impactor_diameter_km = data.get('impactor_diameter', 0)
        
        # Format semi-major axis for label
        if '0p01' in data_file:
            label = '0.01 AU'
        elif '0p05' in data_file:
            label = '0.05 AU'
        elif '0p1' in data_file:
            label = '0.1 AU'
        else:
            label = f'{semi_major_axis:.2f} AU'
        
        # Plot pressure
        ax.semilogy(times, pressures, color=colors['pressure'], linewidth=2)
        
        # Add critical pressure line
        ax.axhline(y=critical_pressure, color=colors['critical'], linestyle='--', 
                  linewidth=1.5, alpha=0.7)
        
        # Add state background
        ymin, ymax = ax.get_ylim()
        x = times
        y1 = np.full_like(x, ymin)
        y2 = np.full_like(x, ymax)
        
        inflated_mask = states == 1
        collapsed_mask = states == 0
        
        if np.any(inflated_mask):
            ax.fill_between(x, y1, y2, where=inflated_mask,
                           color=colors['inflated'], alpha=0.2)
        
        if np.any(collapsed_mask):
            ax.fill_between(x, y1, y2, where=collapsed_mask,
                           color=colors['collapsed'], alpha=0.2)
        
        # Add impact times
        if len(impact_times) > 0:
            ax.vlines(impact_times, ymin=ymin, ymax=ymax,
                     colors=colors['impacts'], linestyles='dashed', 
                     linewidth=1, alpha=0.6)
        
        ax.set_ylabel('Pressure (bar)', fontsize=12)
        ax.grid(True, which="both", ls="-", alpha=0.2)
        
        # Add parameters text in top-left (cleaner style)
        param_text = (
            f"$d = {label}$\n"
            f"$p_c = {critical_pressure:.3f}$ bar\n"
            f"Impactor Ø: {impactor_diameter_km} km"
        )
        ax.text(0.02, 0.95, param_text, transform=ax.transAxes,
               va='top', ha='left', fontsize=11, family='monospace')
    
    axes[-1].set_xlabel('Time (Gyr)', fontsize=12)
    plt.tight_layout()
    
    comparison_file = os.path.join(DEMO_DIR, 'comparison_M0_star.png')
    plt.savefig(comparison_file, dpi=150)
    print(f"Comparison plot saved to: {comparison_file}")
    plt.close()


if __name__ == "__main__":
    main()
