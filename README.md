# Impact-Reinflation Model

A Python package for modeling planetary atmosphere evolution through impacts, atmospheric escape, and outgassing.

## Overview

This package simulates the evolution of planetary atmospheres under the combined effects of:
- **Atmospheric escape** driven by stellar XUV irradiation
- **Outgassing** from the planetary interior
- **Impact events** that deliver volatiles to the atmosphere

The model tracks the transition between **inflated** (gas phase) and **collapsed** (condensed phase) atmospheric states, where collapse occurs when the atmospheric pressure drops below a critical threshold determined by the CO₂ condensation temperature.

## Features

- Simple, clean API for setting up planetary and stellar systems
- Configurable parameters: erosion efficiency, outgassing rates, impact rates
- Single-run and Monte Carlo simulation capabilities
- Built-in plotting utilities for visualization
- No flare functionality (simplified for publication)

## Installation

```bash
# Clone or download the repository
cd impact-reinflation-code

# Install dependencies
pip install -r requirements.txt
```

The package uses only standard scientific Python libraries: numpy, scipy, and matplotlib.

## Quick Start

### Single Run

You can use predefined planet systems from `planet_systems.py`:

```python
from planet_systems import get_planet_system, get_planet_config
from impactmodel import Planet

# Get a predefined system (options: 'TEST_PLANET', 'EARTH_LIKE', 'GJ3929b', 'LTT1445Ac', 'LTT1445Ab')
star, _ = get_planet_system('EARTH_LIKE')
config = get_planet_config('EARTH_LIKE')
planet = Planet(
    name=config['name'],
    mass=config['mass'],  # Earth masses
    radius=config['radius'],  # Earth radii
    star=star,
    a=config['distance'],  # AU
    outgassing_rate=6.5e-4,  # kg/m²/yr
    escape_efficiency=config['escape_efficiency']
)
```

Or create custom planet/star objects:

```python
from impactmodel import Planet, Star

star = Star(L_bol=0.01, Rs=0.25, t_sat=1.3e9)
planet = Planet(
    name="Earth-like",
    mass=1.0,           # Earth masses
    radius=1.0,          # Earth radii
    star=star,
    a=0.1,              # AU
    outgassing_rate=6.5e-4,   # kg/m²/yr
    escape_efficiency=1e-2
)

# Add atmosphere
planet.add_atmosphere(pressure0=0, condensed_volatile_thickness0=0)

# Add impactor
planet.add_default_impactor(impact_rate=1e-9)  # yr^-1

# Evolve the system
results = planet.evolve_system(age=5e9, dt=1e6, method='euler')

print(f"Inflation percentage: {results['inflation_percentage']:.1f}%")
print(f"Number of impacts: {results['num_impacts']}")
```

Or use the provided script:

```bash
python run_single.py
```

Configure the planet system at the top of the script using `PLANET_SYSTEM` (see `planet_systems.py` for available systems).

### Monte Carlo Simulation

Run a parameter sweep over outgassing and impact rates:

```bash
python run_montecarlo.py
```

This will generate results in the `montecarlo_results/` directory. You can configure the planet system and parameter ranges at the top of `run_montecarlo.py`.

### Plotting Results

```bash
# Plot single run results
python plot_results.py results_single.npz output.png

# Plot Monte Carlo summary
python plot_results.py montecarlo_results/montecarlo_summary.csv output.png
```

Or use the plotting functions directly:

```python
from plot_results import plot_single_evolution

plot_single_evolution('results_single.npz', save_path='evolution.png')
```

## User-Configurable Parameters

### In `run_single.py`:

- **Planet system**: Choose from predefined systems in `planet_systems.py` via `PLANET_SYSTEM` variable
- **Atmospheric evolution**: `outgassing_rate`, `escape_efficiency`
- **Impacts**: `impact_rate`
- **Simulation**: `age`, `dt`, `method`

### In `run_montecarlo.py`:

- **Planet system**: Choose from predefined systems in `planet_systems.py` via `PLANET_SYSTEM` variable
- **Parameter ranges**: `outgassing_rate_min/max` (log10 relative to Earth), `impact_rate_min/max` (log10 impacts/year)
- **Number of samples**: `n_samples`
- **Evolution parameters**: `age`, `dt`, `method`
- **Output**: `output_dir`, `output_prefix`

## Model Physics

### Atmospheric States

The atmosphere can exist in two states:

1. **Inflated**: Gas phase, pressure above critical threshold
   - Subject to escape driven by XUV irradiation
   - Can collapse when pressure drops below critical value

2. **Collapsed**: Condensed phase, pressure below critical threshold
   - No escape occurs
   - Outgassing accumulates as condensed volatiles
   - Can be re-inflated by impacts

### Critical Pressure

The critical pressure is determined by solving for where the CO₂ condensation temperature equals the thin-radiator temperature (the equilibrium temperature of a thin CO₂ atmosphere). Below this pressure, CO₂ condenses.

### Impact Events

Impacts deliver kinetic energy that volatilizes condensed material:
- Uses 50% of the impactor's kinetic energy (0.5 efficiency)
- Energy goes into latent heat and heating condensed volatiles
- Only affects condensed volatiles (dayside outgassing)

### Atmospheric Escape

Escape flux is calculated using the stellar XUV emission and escape efficiency parameter:
```
dM/dt = η * π * R³ * F_XUV / (a² * G * M)
```

where `η` is the escape efficiency parameter.

### Outgassing

Constant outgassing rate from the planetary interior:
```
Outgassing = outgassing_rate * surface_area
```

Accumulates as condensed volatiles when atmosphere is collapsed.

## Code Structure

```
impact-reinflation-code/
├── impactmodel/
│   ├── __init__.py       # Package initialization
│   ├── constants.py      # Physical constants
│   ├── planet.py         # Planet class and evolution
│   ├── atmosphere.py     # Atmosphere state and evolution
│   ├── star.py          # Stellar properties
│   └── impactor.py      # Impact event modeling
├── run_single.py         # Single run script
├── run_montecarlo.py     # Monte Carlo script
├── plot_results.py       # Plotting utilities
├── planet_systems.py     # Predefined planet and star systems
├── requirements.txt      # Dependencies
├── demo/                 # Demo scripts for recreating paper plots (not required for use)
│   ├── parameter_space_plot.py  # Multi-planet parameter space exploration
│   └── 3earths_plot.py  # Three Earth-like planets at different distances
└── README.md            # This file
```

**Note:** The `demo/` folder contains scripts used to recreate plots from the publication and is not required for general use of the package. Output directories like `montecarlo_results/` are excluded from version control (see `.gitignore`).

## Examples

See the parameter configurations in:
- `run_single.py` - Single evolution run
- `run_montecarlo.py` - Parameter sweep over outgassing and impact rates

## Output Format

Simulations save results as numpy `.npz` files containing:
- `times`: Time array in years
- `pressures`: Atmospheric pressures in bar
- `atm_masses`: Atmospheric gas masses in kg
- `vol_masses`: Condensed volatile masses in kg
- `states`: State array (0=collapsed, 1=inflated)
- `impact_times`: Times when impacts occur in years
- `inflation_percentage`: Percentage of time in inflated state
- `num_impacts`: Total number of impacts
- Plus simulation parameters

## References

- Wordsworth et al. (2010): CO₂ condensation thermodynamics
- Escape efficiency parameterization based on XUV-driven atmospheric loss

## License

MIT License - see LICENSE file

## Author

Prune Camille August

## Acknowledgments

Part of this code was developed with the assistance of Cursor AI (an AI-powered code editor).
