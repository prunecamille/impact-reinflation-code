"""
Create parameter space plots for GJ 3929 b, LTT 1445 Ac, and LTT 1445 Ab.
Recreates the plot showing atmospheric inflation as a function of 
outgassing rate and impact rate.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import ndimage
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import sys
import os
import csv

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from impactmodel import Planet, Star
from impactmodel.constants import CURRENT_EARTH_OUTGASSING_RATE

def add_smoothed_contours(ax, outgassing_rates, impact_rates, inflation_data,
                         contour_levels=[5, 25, 50], sigma=8.0, color='black'):
    """Add smoothed contours to existing scatter plot."""
    # Create grid for interpolation
    xi = np.logspace(np.log10(outgassing_rates.min()), 
                     np.log10(outgassing_rates.max()), 200)
    yi = np.logspace(np.log10(impact_rates.min()), 
                     np.log10(impact_rates.max()), 200)
    xi_grid, yi_grid = np.meshgrid(xi, yi)
    
    # Interpolate and smooth
    points = np.column_stack((outgassing_rates, impact_rates))
    zi = griddata(points, inflation_data, (xi_grid, yi_grid), 
                  method='cubic', fill_value=0)
    zi = gaussian_filter(zi, sigma=sigma)
    
    # Add contours
    contours = ax.contour(xi_grid, yi_grid, zi, levels=contour_levels, 
                         colors=color, linewidths=2, alpha=1)
    ax.clabel(contours, inline=True, fontsize=18, fmt='$\\mathbf{%.0f\\%%}$',
              colors=color)

# Planet configurations
PLANETS = {
    'GJ3929b': {
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
    },
    'LTT1445Ac': {
        'name': 'LTT1445Ac',
        'mass': 1.37,
        'radius': 1.07,
        'distance': 0.02659,
        'star': {
            'radius': 0.271,
            'luminosity': 10**(-2.04),
            't_sat': 1.9e9
        },
        'escape_efficiency': 1e-2
    },
    'LTT1445Ab': {
        'name': 'LTT1445Ab',
        'mass': 2.73,
        'radius': 1.34,
        'distance': 0.03810,
        'star': {
            'radius': 0.271,
            'luminosity': 10**(-2.04),
            't_sat': 1.9e9
        },
        'escape_efficiency': 1e-2
    }
}

# Parameter ranges for Monte Carlo
# Outgassing rates: log10 relative to Earth's current outgassing rate
outgassing_log_min = -2.0  # 0.01x Earth
outgassing_log_max = 1.0   # ~3x Earth

# Impact rates: log10 impacts per year (matching notebook ranges)
impact_log_min = -12  # log10(impacts/year)
impact_log_max = -6   # log10(impacts/year)

# Number of samples for Monte Carlo
n_samples = 10000  # Random samples across parameter space

# Evolution parameters
age = 12e9  # years (12 billion years)
dt = 1e7    # years per time step
method = 'euler'

# Create output directory
output_dir = os.path.join(os.path.dirname(__file__), 'parameter_space')
os.makedirs(output_dir, exist_ok=True)

print("=" * 80)
print("Parameter Space Monte Carlo Simulation")
print("=" * 80)
print(f"Outgassing range: {outgassing_log_min} to {outgassing_log_max} (log10)")
print(f"Impact rate range: {impact_log_min} to {impact_log_max} (log10 impacts/year)")
print(f"Number of random samples per planet: {n_samples}")
print("=" * 80)

# Run simulations for each planet
all_results = {}

for planet_key, planet_config in PLANETS.items():
    print(f"\nRunning simulations for {planet_key}...")
    
    # Create star
    star_config = planet_config['star']
    star = Star(
        L_bol=star_config['luminosity'],
        Rs=star_config['radius'],
        t_sat=star_config['t_sat']
    )
    
    # Storage for results
    inflation_percentages = []
    impactor_sizes = []
    critical_pressures = []
    outgassing_samples = []
    impact_samples = []
    
    # Generate random Monte Carlo samples
    np.random.seed(42)  # For reproducibility
    outgassing_logs = np.random.uniform(outgassing_log_min, outgassing_log_max, n_samples)
    impact_logs = np.random.uniform(impact_log_min, impact_log_max, n_samples)
    
    outgassing_factors = 10**outgassing_logs  # Convert to linear scale
    impact_rates = 10**impact_logs  # Convert to impacts/year
    
    print(f"  Running {n_samples} Monte Carlo samples...")
    
    for i in range(n_samples):
        if (i+1) % 100 == 0:
            print(f"  Progress: {i+1}/{n_samples} ({100*(i+1)/n_samples:.1f}%)")
        
        # Get random parameters
        outgassing_factor = outgassing_factors[i]
        impact_rate = impact_rates[i]
        
        # Calculate actual outgassing rate
        outgassing_rate = outgassing_factor * CURRENT_EARTH_OUTGASSING_RATE
        
        # Create planet
        planet = Planet(
            name=planet_config['name'],
            mass=planet_config['mass'],
            radius=planet_config['radius'],
            star=star,
            a=planet_config['distance'],
            outgassing_rate=outgassing_rate,
            escape_efficiency=planet_config['escape_efficiency']
        )
        
        # Add atmosphere (starting with no atmosphere)
        planet.add_atmosphere(pressure0=0, condensed_volatile_thickness0=0)
        critical_pressure = planet.atmosphere.critical_pressure * 1e-5  # Convert to bar
        
        # Add impactor
        impactor_size = planet.add_default_impactor(impact_rate)
        
        # Evolve the system
        results = planet.evolve_system(age=age, dt=dt, method=method)
        
        # Store results
        inflation_percentages.append(results['inflation_percentage'])
        impactor_sizes.append(impactor_size)
        critical_pressures.append(critical_pressure)
        outgassing_samples.append(outgassing_factor)
        impact_samples.append(impact_rate)
    
    # Store results for this planet
    all_results[planet_key] = {
        'outgassing_rates': np.array(outgassing_samples),
        'impact_rates': np.array(impact_samples),
        'inflation_percentages': np.array(inflation_percentages),
        'impactor_sizes': np.array(impactor_sizes),
        'critical_pressures': np.array(critical_pressures),
        'planet': planet  # Keep last planet instance for calculations
    }
    
    print(f"  ✓ Completed {planet_key}")
    print(f"    Inflation range: {min(inflation_percentages):.2f}% - {max(inflation_percentages):.2f}%")
    print(f"    Mean inflation: {np.mean(inflation_percentages):.2f}%")

# Save Monte Carlo results to files
print("\nSaving Monte Carlo results...")
for planet_key, data in all_results.items():
    # Save npz file for each planet
    output_file = os.path.join(output_dir, f'{planet_key}_montecarlo.npz')
    np.savez(
        output_file,
        outgassing_rates=data['outgassing_rates'],
        impact_rates=data['impact_rates'],
        inflation_percentages=data['inflation_percentages'],
        impactor_sizes=data['impactor_sizes'],
        critical_pressures=data['critical_pressures']
    )
    print(f"  ✓ Saved {planet_key} results to: {output_file}")
    
    # Save CSV summary for each planet
    csv_file = os.path.join(output_dir, f'{planet_key}_montecarlo_summary.csv')
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'outgassing_rate_factor', 'impact_rate', 'inflation_percentage',
            'impactor_size', 'critical_pressure'
        ])
        for i in range(len(data['outgassing_rates'])):
            writer.writerow([
                data['outgassing_rates'][i],
                data['impact_rates'][i],
                data['inflation_percentages'][i],
                data['impactor_sizes'][i],
                data['critical_pressures'][i]
            ])
    print(f"  ✓ Saved {planet_key} summary to: {csv_file}")

# Create the plot
print("\nGenerating plot...")

fig, axes = plt.subplots(1, 3, figsize=(16, 7), sharex=True, sharey=True)

vmin = 1
vmax = 300
colormap = 'PuRd'

planet_order = ['GJ3929b', 'LTT1445Ac', 'LTT1445Ab']

for idx, planet_key in enumerate(planet_order):
    ax = axes[idx]
    data = all_results[planet_key]
    
    # Use the Monte Carlo samples directly for scatter plot
    outgassing_samples = data['outgassing_rates']
    impact_samples = data['impact_rates']
    
    # Planet name (bottom left)
    planet_name = data['planet'].name
    if planet_key == 'GJ3929b':
        ax.text(0.015, 0.002, f'GJ 3929 b', fontsize=18, ha='left', va='bottom')
    elif planet_key == 'LTT1445Ac':
        ax.text(0.015, 0.002, f'LTT 1445 Ac', fontsize=18, ha='left', va='bottom')
    else:
        ax.text(0.015, 0.002, f'LTT 1445 Ab', fontsize=18, ha='left', va='bottom')
    
    # Convert impact rates to impacts/Gyr (1 Gyr = 1e9 years)
    impact_samples_gyr = impact_samples * 1e9  # Convert from impacts/year to impacts/Gyr
    
    # Scatter plot with correct units
    sc = ax.scatter(
        outgassing_samples, impact_samples_gyr,
        c=data['inflation_percentages'],
        cmap=colormap,
        s=4,
        norm=LogNorm(vmin=vmin, vmax=vmax)
    )
    
    # Set logarithmic scales
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Labels with colors
    if idx == 0:
        ax.set_xlabel('CO$_2$ outgassing rate (x modern Earth)', fontsize=15, color='xkcd:scarlet')
        ax.set_ylabel('Impact rate (impacts/Gyr)', fontsize=15, color='xkcd:dark aqua')
    
    # Color the spines
    ax.spines['bottom'].set_color('xkcd:scarlet')
    ax.spines['left'].set_color('xkcd:dark aqua')
    
    # Ticks
    ax.set_xticks([0.05, 1, 3])
    ax.set_yticks([0.01, 0.1, 1, 10, 100, 1000])
    
    # Tick labels with colors
    ax.set_xticklabels(['0.05', '1', '3'], color='xkcd:scarlet')
    ax.set_yticklabels(['0.01', '0.1', '1', '10', '100', '1000'], color='xkcd:dark aqua')
    
    # Grid
    ax.grid(True, which="both", ls="-", alpha=0.2)
    
    # Limits
    ax.set_xlim(0.01, 3)
    ax.set_ylim(0.001, 1000)
    
    # Add smoothed contours (using impacts/Gyr for contour calculation)
    add_smoothed_contours(
        ax, outgassing_samples, impact_samples_gyr, data['inflation_percentages'],
        contour_levels=[5, 25, 50], sigma=8.0, color='black'
    )
    
    # Add top x-axis for cumulative outgassed CO2
    planet = data['planet']
    # Calculate cumulative outgassed CO2 inventory as percentage of planetary mass
    # outgassing_samples are already relative to Earth, convert to actual rate
    total_inventory = (outgassing_samples * CURRENT_EARTH_OUTGASSING_RATE * 
                      planet.surface_area * age / planet.mass)
    
    x_axis_top = ax.twiny()
    x_axis_top.set_xscale('log')
    x_axis_top.set_xlim(total_inventory.min(), total_inventory.max())
    
    # Change colors of the axes
    x_axis_top.spines['top'].set_color('xkcd:royal blue')
    x_axis_top.spines['bottom'].set_color('xkcd:scarlet')
    x_axis_top.spines['left'].set_color('xkcd:dark aqua')
    
    # Get ticks and set labels
    xtopticks = x_axis_top.get_xticks()
    x_axis_top.set_xticklabels([f'{tick*1e3:.2f}' if tick > 1e-5 else f'{tick*1e3:.1f}' for tick in xtopticks], 
                               color='xkcd:royal blue', fontsize=15)
    
    # Set title
    if idx == 0:
        x_axis_top.set_title('Cumulative outgassed CO$_2$ (% of planetary mass)', 
                            color='xkcd:royal blue', fontsize=15)

# Add colorbar at the top
cbar_ax = fig.add_axes([0.45, 0.85, 0.4, 0.035])
cbar = fig.colorbar(sc, cax=cbar_ax, orientation='horizontal')
cbar.set_label('Atmospheric inflation (%)', fontsize=15, labelpad=10)
cbar.ax.xaxis.set_label_position('top')
cbar.ax.xaxis.set_ticks_position('top')

# Adjust layout
plt.subplots_adjust(top=0.8, bottom=0.15, left=0.08, right=0.95, hspace=0.05, wspace=0.05)

# Save figure
output_file = os.path.join(output_dir, 'parameter_space_plot.pdf')
plt.savefig(output_file, dpi=300)
print(f"✓ Plot saved to: {output_file}")
plt.close()

print("\n" + "=" * 80)
print("Parameter space plot complete!")
print("=" * 80)

