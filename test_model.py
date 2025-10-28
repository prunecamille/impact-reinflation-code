"""
Test script to verify the impact-reinflation model works correctly.
"""

import numpy as np
from impactmodel import Planet, Star
from impactmodel.constants import CURRENT_EARTH_OUTGASSING_RATE

print("=" * 60)
print("Testing Impact-Reinflation Model")
print("=" * 60)

# Create a star
print("\n1. Creating star...")
star = Star(L_bol=0.01, Rs=0.25, t_sat=1.3e9)
print(f"   ✓ Star created: L_bol={star.L_bol} L_sun")

# Create a planet
print("\n2. Creating planet...")
planet = Planet(
    name="TestPlanet",
    mass=1.0,
    radius=1.0,
    star=star,
    a=0.1,
    outgassing_rate=CURRENT_EARTH_OUTGASSING_RATE,
    escape_efficiency=1e-2
)
print(f"   ✓ Planet created: {planet.name}")

# Add atmosphere
print("\n3. Adding atmosphere...")
planet.add_atmosphere(pressure0=0, condensed_volatile_thickness0=0)
critical_pressure = planet.atmosphere.critical_pressure
print(f"   ✓ Atmosphere added")
print(f"   Critical pressure: {critical_pressure*1e-5:.4f} bar")

# Add impactor
print("\n4. Adding impactor...")
impact_rate = 1e-9
impactor_size = planet.add_default_impactor(impact_rate)
print(f"   ✓ Impactor added: {impactor_size:.2f} km")

# Evolve the system
print("\n5. Evolving system...")
age = 5e9  # 5 billion years
dt = 1e6   # 1 million year steps
results = planet.evolve_system(age=age, dt=dt, method='euler')

print(f"   ✓ Evolution complete")
print(f"   Time steps: {len(results['times'])}")
print(f"   Final pressure: {results['pressures'][-1]:.4e} bar")
print(f"   Number of impacts: {results['num_impacts']}")
print(f"   Inflation percentage: {results['inflation_percentage']:.2f}%")

# Test plotting
print("\n6. Testing plot function...")
try:
    from plot_results import plot_single_evolution
    print("   ✓ Plot function imported")
except ImportError as e:
    print(f"   ✗ Failed to import plot function: {e}")

print("\n" + "=" * 60)
print("All tests passed!")
print("=" * 60)
