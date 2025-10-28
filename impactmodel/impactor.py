"""
Module for handling impactor properties and stochastic impact events.

This module provides the Impactor class, which models the physical properties
of impactors (size, mass, velocity) and determines when impacts occur based
on a given impact rate.
"""

import numpy as np


class Impactor:
    """
    Class representing an impactor object that can impact a planet.
    
    The impactor has physical properties (diameter, mass, velocity) and an
    impact rate that determines when impacts occur using Poisson statistics.
    
    Attributes:
        diameter (float): Impactor diameter in meters
        impact_rate (float): Rate of impacts in years^-1
        velocity (float): Impact velocity in m/s
        density (float): Impactor density in kg/m^3
        mass (float): Impactor mass in kg
        kinetic_energy (float): Kinetic energy in Joules
    """
    
    def __init__(self, diameter: float, impact_rate: float, velocity: float, 
                 density: float = 3000):
        """
        Initialize an Impactor instance.
        
        Args:
            diameter (float): Impactor diameter in kilometers
            impact_rate (float): Impact rate in yr^-1
            velocity (float): Impact velocity in m/s
            density (float, optional): Impactor density in kg/m^3. Defaults to 3000.
        """
        # Store properties
        self.impact_rate = impact_rate  # yr^-1
        self.diameter = diameter * 1e3  # Convert km to m
        self.velocity = velocity  # m/s
        self.density = density  # kg/m^3
        
        # Calculate derived properties
        self.mass = 4/3 * np.pi * (self.diameter / 2)**3 * self.density  # kg
        self.kinetic_energy = 0.5 * self.mass * self.velocity**2  # J
    
    def time_to_next_impact(self) -> float:
        """
        Calculate the time until the next impact event.
        
        Uses an exponential distribution based on the impact rate (Poisson process).
        
        Returns:
            float: Time in years until next impact
        """
        if self.impact_rate > 0:
            return -np.log(np.random.uniform(0, 1)) / self.impact_rate
        else:
            return np.inf
