"""
Module for handling stellar properties and XUV emission evolution.

This module provides the Star class, which models stellar properties with a 
focus on XUV emission and its evolution over time due to stellar activity decay.
"""

from typing import Optional
import numpy as np
from .constants import L_SUN, R_SUN, SIGMA_BOLTZMANN


class Star:
    """
    A class representing a star and its properties relevant for planetary atmosphere evolution.
    
    This class models stellar properties with a focus on XUV emission and its 
    evolution over time. The XUV emission follows a power-law decay after a 
    saturation time when stellar activity decreases.
    
    Attributes:
        L_bol (float): Bolometric luminosity in solar units
        Rs (float): Stellar radius in solar radii
        t_sat (float): Saturation time in years (when XUV emission begins to decay)
        L_XUV0 (float): Initial XUV luminosity in W
        planets (List[Planet]): List of associated planets
    """
    
    def __init__(self, L_bol: float, Rs: float, t_sat: float = 1e9):
        """
        Initialize a Star instance.
        
        Args:
            L_bol (float): Bolometric luminosity in solar units
            Rs (float): Stellar radius in solar radii
            t_sat (float, optional): Saturation time in years. Defaults to 1e9.
            
        Raises:
            ValueError: If L_bol or Rs are non-positive
        """
        if L_bol <= 0 or Rs <= 0:
            raise ValueError("L_bol and Rs must be positive")
            
        # Store properties
        self.L_bol = L_bol
        self.Rs = Rs
        self.t_sat = t_sat
        
        # Convert to SI units
        self.L_bol_SI = L_bol * L_SUN  # W
        self.Rs_SI = Rs * R_SUN  # m
        
        # Calculate initial XUV luminosity (10^-3.5 of bolometric luminosity)
        self.L_XUV0 = 10**(-3.5) * self.L_bol_SI  # W
        
        self.planets = []
    
    def assign_planet(self, planet: 'Planet') -> None:
        """
        Assign a planet to this star.
        
        Args:
            planet (Planet): The planet to assign
        """
        self.planets.append(planet)
    
    def XUV_emission(self, t: float) -> float:
        """
        Calculate the XUV emission at a given time.
        
        The XUV emission is constant before saturation time (t_sat) and follows
        a power law decay after saturation: L_XUV(t) = L_XUV0 * (t/t_sat)^(-1.23)
        
        Args:
            t (float): Time in years
            
        Returns:
            float: XUV luminosity in W
        """
        if t <= self.t_sat:
            return self.L_XUV0
        else:
            return self.L_XUV0 * (t/self.t_sat)**(-1.23)
    
    def F_XUV(self, t: float, planet: 'Planet') -> float:
        """
        Calculate the XUV flux at a specific planet's position at a given time.
        
        Args:
            t (float): Time in years
            planet (Planet): The planet to calculate flux for
            
        Returns:
            float: XUV flux in W/m²
        """
        L_XUV = self.XUV_emission(t)
        return L_XUV / (4 * np.pi * planet.a**2)
    
    def F(self, planet: 'Planet') -> float:
        """
        Calculate the bolometric flux at a specific planet's position.
        
        Args:
            planet (Planet): The planet to calculate flux for
            
        Returns:
            float: Bolometric flux in W/m²
        """
        return self.L_bol_SI / (4 * np.pi * planet.a**2)
    
    def integrated_XUV(self, t: float, planet: 'Planet') -> float:
        """
        Calculate the integrated XUV flux over a planet's lifetime.
        
        Args:
            t (float): Time in years
            planet (Planet): The planet to calculate integrated flux for
            
        Returns:
            float: Integrated XUV flux in J/m²
        """
        if t <= self.t_sat:
            integral = self.L_XUV0 * t * (365.25 * 24 * 60 * 60)  # Convert years to seconds
        else:
            # For t > t_sat, integrate the power law decay analytically
            # The integral of (t/t_sat)^(-1.23) from t_sat to t is evaluated
            seconds_per_year = 365.25 * 24 * 60 * 60
            integral = self.L_XUV0 * self.t_sat * seconds_per_year + \
                      self.L_XUV0 * self.t_sat * seconds_per_year * \
                      (1 - (t/self.t_sat)**-0.23) / 0.23
            
        return integral / (4 * np.pi * planet.a**2)
