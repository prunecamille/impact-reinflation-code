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
from typing import Optional

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


def plot_montecarlo_summary(csv_file: str, save_path: Optional[str] = None) -> None:
    """
    Plot summary statistics from Monte Carlo results.
    
    Creates a 2-panel plot showing:
    - Heat map of inflation percentage vs outgassing and impact rates
    - Scatter plot of final atmospheric mass vs parameters
    
    Args:
        csv_file (str): Path to the CSV file with summary results
        save_path (str, optional): Path to save the figure. If None, displays the plot.
    """
    import pandas as pd
    
    # Load data
    df = pd.read_csv(csv_file)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Panel 1: Heat map of inflation percentage
    # Reshape data for heatmap
    outgassing_unique = np.unique(df['outgassing_rate'])
    impact_unique = np.unique(df['impact_rate'])
    
    inflation_grid = np.zeros((len(outgassing_unique), len(impact_unique)))
    for i, og in enumerate(outgassing_unique):
        for j, ir in enumerate(impact_unique):
            mask = (df['outgassing_rate'] == og) & (df['impact_rate'] == ir)
            if mask.any():
                inflation_grid[i, j] = df[mask]['inflation_percentage'].values[0]
    
    im1 = ax1.imshow(inflation_grid, aspect='auto', cmap='viridis', 
                     origin='lower')
    ax1.set_xlabel('Impact Rate')
    ax1.set_ylabel('Outgassing Rate')
    ax1.set_title('Inflation Percentage (%)')
    
    # Set tick labels
    ax1.set_xticks(range(len(impact_unique)))
    ax1.set_xticklabels([f'{ir:.2e}' for ir in impact_unique], rotation=45)
    ax1.set_yticks(range(len(outgassing_unique)))
    ax1.set_yticklabels([f'{og/CURRENT_EARTH_OUTGASSING_RATE:.1f}x' 
                         for og in outgassing_unique])
    
    plt.colorbar(im1, ax=ax1)
    
    # Panel 2: Final atmospheric mass scatter plot
    scatter = ax2.scatter(np.log10(df['impact_rate']), 
                         np.log10(df['outgassing_rate']),
                         c=df['inflation_percentage'], 
                         s=50, cmap='viridis', alpha=0.7)
    ax2.set_xlabel('log10(Impact Rate)')
    ax2.set_ylabel('log10(Outgassing Rate / Earth rate)')
    ax2.set_title('Final Atmospheric State')
    plt.colorbar(scatter, ax=ax2, label='Inflation %')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150)
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
