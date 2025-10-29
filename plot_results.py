"""
Plotting utilities for the impact-reinflation model results.

This module provides functions to visualize the results of single runs
and Monte Carlo simulations from the impact-reinflation model.

Functions:
- plot_single_evolution: Plot evolution of a single simulation
- plot_montecarlo_summary: Plot summary statistics from Monte Carlo results

Usage:
    from plot_results import plot_single_evolution
    plot_single_evolution('results_single.npz')
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from typing import Optional, Union

# Set default plot style
rcParams['font.size'] = 12
rcParams['axes.labelsize'] = 12
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 10
rcParams['figure.figsize'] = [10, 6]


def plot_single_evolution(data_file: str, save_path: Optional[str] = None) -> None:
    """
    Plot the evolution of atmospheric properties from a single simulation.
    
    Creates a 2-panel plot showing:
    - Top panel: Atmospheric pressure over time with state shading
    - Bottom panel: Atmospheric and condensed volatile masses over time
    
    Args:
        data_file (str): Path to the .npz file containing the results
        save_path (str, optional): Path to save the figure. If None, displays the plot.
    """
    # Load data
    data = np.load(data_file)
    
    times = data['times'] / 1e6  # Convert to Myr
    pressures = data['pressures']
    atm_masses = data['atm_masses']
    vol_masses = data['vol_masses']
    states = data['states']
    impact_times = data['impact_times'] / 1e6  # Convert to Myr
    
    # Define colors
    colors = {
        'inflated': '#06c2ac',  # Teal
        'collapsed': '#742E1E',  # Chestnut
        'pressure': '#380835',  # Deep purple
        'atm_mass': '#be0119',  # Deep red
        'vol_mass': '#6a79f7',  # Soft blue
        'impacts': '#80013f'  # Deep purple
    }
    
    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), 
                                    sharex=True, 
                                    gridspec_kw={'height_ratios': [3, 1]})
    
    # Get critical pressure if available
    critical_pressure = data.get('critical_pressure', None)
    
    # Top panel: Pressure with state background
    ax1.semilogy(times, pressures, color=colors['pressure'], linewidth=2)
    
    # Add critical pressure line if available
    if critical_pressure is not None:
        ax1.axhline(y=critical_pressure, color='red', linestyle='--', 
                    linewidth=1.5, alpha=0.7, label=f'Crit P = {critical_pressure:.4f} bar')
    
    # Add state background
    ymin, ymax = ax1.get_ylim()
    x = times
    y1 = np.full_like(x, ymin)
    y2 = np.full_like(x, ymax)
    
    inflated_mask = states == 1
    collapsed_mask = states == 0
    
    if np.any(inflated_mask):
        ax1.fill_between(x, y1, y2, where=inflated_mask,
                         color=colors['inflated'], alpha=0.2, 
                         label='Inflated State')
    
    if np.any(collapsed_mask):
        ax1.fill_between(x, y1, y2, where=collapsed_mask,
                         color=colors['collapsed'], alpha=0.2, 
                         label='Collapsed State')
    
    ax1.set_ylabel('Pressure (bar)')
    ax1.grid(True, which="both", ls="-", alpha=0.2)
    
    # Add impact times
    if len(impact_times) > 0:
        ax1.vlines(impact_times, ymin=ymin, ymax=ymax,
                  colors=colors['impacts'], linestyles='dashed', 
                  alpha=0.5, label='Impacts')
    
    # Bottom panel: Masses
    ax2.semilogy(times, atm_masses, color=colors['atm_mass'], 
                 linewidth=2, label='Atmospheric Mass')
    ax2_twin = ax2.twinx()
    ax2_twin.semilogy(times, vol_masses, color=colors['vol_mass'], 
                      linewidth=2, label='Condensed Mass')
    
    ax2.set_xlabel('Time (Myr)')
    ax2.set_ylabel('Atmospheric Mass (kg)', color=colors['atm_mass'])
    ax2_twin.set_ylabel('Condensed Mass (kg)', color=colors['vol_mass'])
    
    ax2.tick_params(axis='y', labelcolor=colors['atm_mass'])
    ax2_twin.tick_params(axis='y', labelcolor=colors['vol_mass'])
    
    ax2.grid(True, which="both", ls="-", alpha=0.2)
    
    # Add impact times to bottom panel
    if len(impact_times) > 0:
        ax2.vlines(impact_times, 
                  ymin=min(min(atm_masses), min(vol_masses)),
                  ymax=max(max(atm_masses), max(vol_masses)),
                  colors=colors['impacts'], linestyles='dashed', alpha=0.5)
    
    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        mpatches.Patch(facecolor=colors['inflated'], alpha=0.2, 
                      label='Inflated State'),
        mpatches.Patch(facecolor=colors['collapsed'], alpha=0.2, 
                      label='Collapsed State'),
        Line2D([0], [0], color=colors['impacts'], linestyle='dashed', 
              label='Impact')
    ]
    
    # Add critical pressure to legend if available
    if critical_pressure is not None:
        legend_elements.append(Line2D([0], [0], color='red', linestyle='--', 
                                    label=f'Crit P = {critical_pressure:.4f} bar'))
    
    ax1.legend(handles=legend_elements, loc='upper left', fontsize=9)
    
    # Add summary text
    crit_text = f" | Crit P: {critical_pressure:.4f} bar" if critical_pressure is not None else ""
    summary_text = (
        f"Inflation: {data['inflation_percentage']:.1f}% | "
        f"Impacts: {data['num_impacts']} | "
        f"Final P: {pressures[-1]:.2e} bar{crit_text}"
    )
    ax1.text(0.02, 0.98, summary_text, transform=ax1.transAxes,
             va='top', ha='left', bbox=dict(boxstyle='round', 
             facecolor='wheat', alpha=0.5), fontsize=9)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Figure saved to: {save_path}")
    else:
        plt.show()
    
    plt.close()


def add_smoothed_contours(ax, outgassing_rates, impact_rates, inflation_data,
                         contour_levels=[1, 10, 30], sigma=5.0, color='black'):
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


def plot_montecarlo_summary(csv_file: str, save_path: Optional[str] = None,
                            planet: Optional[Union['Planet', dict]] = None,
                            age: Optional[float] = None) -> None:
    """
    Plot summary statistics from Monte Carlo results in notebook style.
    
    Creates a single-panel scatter plot matching the style from the impact paper notebook.
    
    Args:
        csv_file (str): Path to the CSV file with summary results
        save_path (str, optional): Path to save the figure. If None, displays the plot.
        planet (Planet or dict, optional): Planet object or dict with planet properties.
            If provided, adds top x-axis showing total CO2 inventory. Dict should have
            'name', 'mass', 'radius', and 'surface_area' (or will use Planet object methods).
        age (float, optional): Age of simulation in years. Needed for top x-axis calculation.
    """
    import pandas as pd
    from impactmodel.constants import CURRENT_EARTH_OUTGASSING_RATE
    
    # Load data
    df = pd.read_csv(csv_file)
    
    # Get outgassing rates (linear factors x Earth) - matching parameter_space_plot.py format
    # CSV column is 'outgassing_rate_factor' (from run_montecarlo.py)
    if 'outgassing_rate_factor' in df.columns:
        outgassing_samples = df['outgassing_rate_factor'].values
    elif 'outgassing_factor' in df.columns:
        outgassing_samples = df['outgassing_factor'].values
    elif 'outgassing_rate' in df.columns:
        # Fallback: assume it's actual rate, convert to factor
        outgassing_samples = df['outgassing_rate'].values / CURRENT_EARTH_OUTGASSING_RATE
    else:
        raise ValueError("CSV must contain 'outgassing_rate_factor', 'outgassing_factor', or 'outgassing_rate' column")
    
    # Get impact rates and convert to impacts/Gyr
    impact_samples = df['impact_rate'].values
    impact_samples_gyr = impact_samples * 1e9  # Convert from impacts/year to impacts/Gyr
    
    # Get inflation percentages
    inflation_percentages = df['inflation_percentage'].values
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    
    # Parameters matching parameter_space_plot.py
    vmin = 1
    vmax = 300
    colormap = 'PuRd'
    
    # Planet name (bottom left)
    if planet is not None:
        if hasattr(planet, 'name'):
            planet_name = planet.name
        elif isinstance(planet, dict):
            planet_name = planet.get('name', 'Planet')
        else:
            planet_name = 'Planet'
        ax.text(0.015, 0.002, planet_name, fontsize=18, ha='left', va='bottom')
    
    # Scatter plot with correct units
    sc = ax.scatter(
        outgassing_samples, impact_samples_gyr,
        c=inflation_percentages,
        cmap=colormap, s=4,
        norm=LogNorm(vmin=vmin, vmax=vmax)
    )
    
    
    # Set logarithmic scales
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # Labels with colors
    ax.set_xlabel('CO$_2$ outgassing rate (x modern Earth)', fontsize=15, color='xkcd:scarlet', labelpad=1)
    ax.set_ylabel('Impact rate (impacts/Gyr)', fontsize=15, color='xkcd:dark aqua', labelpad=1)
    
    # Color the spines
    ax.spines['bottom'].set_color('xkcd:scarlet')
    ax.spines['left'].set_color('xkcd:dark aqua')
    
    # Ticks
    ax.set_xticks([0.05, 1, 3])
    ax.set_yticks([0.1, 1, 10, 100, 1000])
    
    # Tick labels with colors
    ax.set_xticklabels(['0.05', '1', '3'], color='xkcd:scarlet', fontsize=15)
    ax.set_yticklabels(['0.1', '1', '10', '100', '1000'], color='xkcd:dark aqua', fontsize=15)
    
    # Grid
    ax.grid(True, which="both", ls="-", alpha=0.2)
    
    # Limits (matching parameter_space_plot.py)
    ax.set_xlim(0.03, 3)
    ax.set_ylim(0.1, 1000)
    
    # Add smoothed contours (using impacts/Gyr for contour calculation)
    add_smoothed_contours(
        ax, outgassing_samples, impact_samples_gyr, inflation_percentages,
        contour_levels=[5, 25, 50], sigma=8.0, color='black'
    )
    
    # Add top x-axis for cumulative outgassed CO2
    # Use provided planet/age or defaults if not provided
    if planet is None:
        # Create default planet using typical values
        from impactmodel.planet import Planet as PlanetClass
        from impactmodel.star import Star
        default_star = Star(L_bol=0.01, Rs=0.25, t_sat=1.3e9)
        planet_obj = PlanetClass(
            name='Planet',
            mass=1.0,  # Earth masses
            radius=1.0,  # Earth radii
            star=default_star,
            a=0.05  # AU
        )
        planet_obj.add_atmosphere(pressure0=0, condensed_volatile_thickness0=0)
    elif isinstance(planet, dict):
        from impactmodel.planet import Planet as PlanetClass
        from impactmodel.star import Star
        temp_star = Star(L_bol=0.01, Rs=0.25, t_sat=1.3e9)
        planet_obj = PlanetClass(
            name=planet.get('name', 'Planet'),
            mass=planet.get('mass', 1.0),
            radius=planet.get('radius', 1.0),
            star=temp_star,
            a=planet.get('distance', 0.05)
        )
        planet_obj.add_atmosphere(pressure0=0, condensed_volatile_thickness0=0)
    else:
        planet_obj = planet
    
    # Use provided age or default
    if age is None:
        age = 5e9  # Default 5 billion years
    
    # Calculate cumulative outgassed CO2 inventory as percentage of planetary mass
    # outgassing_samples are already relative to Earth (linear factors), convert to actual rate
    total_inventory = (outgassing_samples * CURRENT_EARTH_OUTGASSING_RATE * 
                      planet_obj.surface_area * age / planet_obj.mass)
    
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
    x_axis_top.set_title('Cumulative outgassed CO$_2$ (% of planetary mass)', 
                        color='xkcd:royal blue', fontsize=15)
    
    # Add vertical colorbar on the right side
    cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.75])  # [left, bottom, width, height]
    cbar = fig.colorbar(sc, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Atmospheric inflation (%)', fontsize=15, labelpad=1)
    
    # Adjust layout to make room for colorbar on right and left labels
    plt.subplots_adjust(top=0.9, bottom=0.15, left=0.13, right=0.88)
    
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Figure saved to: {save_path}")
    else:
        plt.show()
    
    plt.close()


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python plot_results.py <data_file> [save_path]")
        print("  data_file: Path to .npz or .csv results file")
        print("  save_path: Optional path to save the figure")
        sys.exit(1)
    
    data_file = sys.argv[1]
    save_path = sys.argv[2] if len(sys.argv) > 2 else None
    
    if data_file.endswith('.npz'):
        plot_single_evolution(data_file, save_path)
    elif data_file.endswith('.csv'):
        plot_montecarlo_summary(data_file, save_path)
    else:
        print("Error: File must be .npz or .csv")
        sys.exit(1)
