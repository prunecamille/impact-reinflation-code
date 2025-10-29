"""
Module for handling atmospheric properties and evolution.

This module provides the Atmosphere class, which models atmospheric evolution 
with a focus on the interaction between gas and condensed phases of volatiles. 
It handles transitions between inflated and collapsed states.
"""

from typing import Optional
import numpy as np
from .constants import LATENT_HEAT_CO2, CP_CO2, T_SUBLIMATION_CO2


class Atmosphere:
    """
    A class representing a planet's atmosphere and its evolution.
    
    This class models the atmosphere's state and evolution, including the 
    interaction between gas and condensed phases of volatiles. It handles 
    transitions between inflated and collapsed states, and processes impacts 
    that affect the atmosphere.
    
    Attributes:
        planet (Planet): The planet this atmosphere belongs to
        pressure (float): Current atmospheric pressure in Pa
        condensed_volatile_thickness (float): Thickness of condensed volatiles in m
        condensed_volatile_mass (float): Mass of condensed volatiles in kg
        gas_mass (float): Mass of atmospheric gas in kg
        critical_pressure (float): Pressure threshold for collapse in Pa
        state (str): Current state ('inflated' or 'collapsed')
    """
    
    def __init__(self, planet: 'Planet', pressure0: float, 
                 condensed_volatile_thickness0: float, 
                 critical_pressure: Optional[float] = None):
        """
        Initialize an atmosphere with given properties.
        
        Args:
            planet: Planet object this atmosphere belongs to
            pressure0: Initial pressure in bar
            condensed_volatile_thickness0: Initial condensed volatile thickness in m
            (Like if you had a layer of ice on the nightside)
            critical_pressure: Critical pressure in bar. If None, will be calculated.
        """
        self.initialparams = {'pressure': pressure0, 
                             'condensed_volatile_thickness': condensed_volatile_thickness0}
        self._planet = planet
        self.set_volatile_properties()  # Set CO2 properties
        
        # Set condensed volatile properties
        self.condensed_volatile_thickness = condensed_volatile_thickness0  # m
        # Condensed volatiles only on nightside (divide by 2)
        self.condensed_volatile_mass = (condensed_volatile_thickness0 * 
                                        self._planet.surface_area / 2)  # kg
        
        # Set critical pressure for atmospheric collapse
        if critical_pressure is None:
            self.critical_pressure = self._planet.solve_for_critical_pressure()  # Pa
        else:
            self.critical_pressure = critical_pressure * 1e5  # Convert bar to Pa
        
        self.critical_mass = self.pressure_to_mass(self.critical_pressure)  # kg
        
        # Initialize pressure and gas mass
        self.pressure = pressure0 * 1e5  # Convert bar to Pa
        self.gas_mass = self.pressure_to_mass(self.pressure)  # kg
        
        # Set initial state based on pressure
        if self.pressure > self.critical_pressure:
            self.state = 'inflated'
            # If starting inflated, convert any initial condensed mass to gas
            self.gas_mass += self.condensed_volatile_mass
            self.condensed_volatile_mass = 0
        else:
            self.state = 'collapsed'
            # If starting collapsed, add gas mass to condensed mass
            self.condensed_volatile_mass += self.gas_mass
            self.gas_mass = 0
            self.pressure = 0

    @classmethod
    def from_gas_mass(cls, planet: 'Planet', gas_mass0: float, 
                      condensed_volatile_thickness0: float,
                      critical_pressure: Optional[float] = None) -> 'Atmosphere':
        """
        Alternative constructor to initialize atmosphere with gas mass instead of pressure.
        
        Args:
            planet: Planet object this atmosphere belongs to
            gas_mass0: Initial gas mass in kg
            condensed_volatile_thickness0: Initial condensed volatile thickness in m
            critical_pressure: Critical pressure in bar. If None, will be calculated.
            
        Returns:
            Atmosphere: A new Atmosphere instance initialized with the given gas mass
        """
        # Convert gas mass to pressure in bar for the standard constructor
        pressure0 = (planet.surface_gravity * gas_mass0 / 
                    (4 * np.pi * planet.radius**2))  # Pa
        pressure0 = pressure0 / 1e5  # Convert Pa to bar
        
        return cls(planet, pressure0, condensed_volatile_thickness0, critical_pressure)

    def update_state(self) -> None:
        """
        Update the state of the atmosphere based on current pressure.
        
        Checks if the atmosphere should transition between inflated and 
        collapsed states based on the current pressure relative to the 
        critical pressure.
        """
        self.pressure = self.mass_to_pressure(self.gas_mass)
        if self.state == 'inflated':
            if self.pressure <= self.critical_pressure:
                self.collapse()
        elif self.state == 'collapsed':
            if self.pressure > self.critical_pressure:
                self.inflate()
    
    def inflate(self) -> None:
        """
        Transition the atmosphere to the inflated state.
        
        In the inflated state, all condensed volatiles become gas.
        """
        # All condensed volatile becomes gas
        self.gas_mass += self.condensed_volatile_mass
        self.condensed_volatile_mass = 0
        self.pressure = self.mass_to_pressure(self.gas_mass)
        self.state = 'inflated'
    
    def collapse(self) -> None:
        """
        Transition the atmosphere to the collapsed state.
        
        In the collapsed state, all gas becomes condensed volatiles.
        """
        if self.gas_mass > 0:
            self.condensed_volatile_mass += self.gas_mass
            self.gas_mass = 0
            self.pressure = 0
            self.state = 'collapsed'
        else:
            self.pressure = 0
            self.gas_mass = 0
            self.state = 'collapsed'
    
    def process_impact(self, mass_volatilized: float) -> None:
        """
        Process an impact event, converting condensed volatiles to gas.
        
        Args:
            mass_volatilized (float): Mass of volatiles to convert in kg
        """
        # Convert condensed volatiles to gas
        mass_volatilized = min(mass_volatilized, self.condensed_volatile_mass)
        self.condensed_volatile_mass -= mass_volatilized
        self.gas_mass += mass_volatilized
        self.pressure = self.mass_to_pressure(self.gas_mass)
        # Check for state transition
        if self.pressure > self.critical_pressure:
            self.inflate()
        else:
            self.collapse()
    
    def reset(self) -> None:
        """
        Reset the atmosphere to its initial parameters.
        """
        self.pressure = self.initialparams['pressure'] * 1e5  # Convert bar to Pa
        self.condensed_volatile_thickness = self.initialparams['condensed_volatile_thickness']
        self.condensed_volatile_mass = (self.initialparams['condensed_volatile_thickness'] * 
                                       self._planet.surface_area / 2)
        self.gas_mass = self.pressure_to_mass(self.pressure)
        self.state = 'inflated' if self.pressure > self.critical_pressure else 'collapsed'
    
    def set_volatile_properties(self, volatile_latent_heat: float = LATENT_HEAT_CO2,
                              volatile_cp: float = CP_CO2,
                              volatile_sublimation_temp: float = T_SUBLIMATION_CO2) -> None:
        """
        Set the properties of the volatile (by default CO2).
        
        Args:
            volatile_latent_heat (float): Latent heat of sublimation in J/kg
            volatile_cp (float): Specific heat capacity in J/kg/K
            volatile_sublimation_temp (float): Sublimation temperature in K
        """
        self.volatile_latent_heat = volatile_latent_heat
        self.volatile_cp = volatile_cp
        self.volatile_sublimation_temp = volatile_sublimation_temp
    
    def pressure_to_mass(self, pressure: float) -> float:
        """
        Convert pressure to gas mass.
        
        Args:
            pressure (float): Pressure in Pa
            
        Returns:
            float: Gas mass in kg
        """
        return pressure * 4 * np.pi * self._planet.radius**2 / self._planet.surface_gravity
    
    def mass_to_pressure(self, mass: float) -> float:
        """
        Convert gas mass to pressure.
        
        Args:
            mass (float): Gas mass in kg
            
        Returns:
            float: Pressure in Pa
        """
        return mass * self._planet.surface_gravity / (4 * np.pi * self._planet.radius**2)
    
    def evolve_euler(self, t: float, dt: float) -> None:
        """
        Evolve the atmosphere over a given time step using Euler's method.
        
        In inflated state:
        - Gas mass changes due to outgassing and escape
        - Total mass decreases because escape removes mass from the system
        
        In collapsed state:
        - No escape occurs
        - Outgassing goes into condensed volatiles
        - Inflation can ONLY occur from impacts
        
        Args:
            t (float): Current time in years
            dt (float): Time step in years
        """
        if self.state == 'inflated':
            # Euler equation format: dM/dt = outgassing - escape
            self.gas_mass += dt * (self._planet.outgassed_mass_flux(t) - 
                                  self._planet.escape_mass_flux(t))
            # If gas mass is negative, set it to 0 and collapse
            if self.gas_mass < 0:
                self.gas_mass = 0
                self.collapse()
        elif self.state == 'collapsed':
            # Only outgassing occurs in collapsed state, goes to condensed volatiles
            dM = dt * self._planet.outgassed_mass_flux(t)
            self.condensed_volatile_mass += dM
        
        self.update_state()
    
    def evolve_rk4(self, t: float, dt: float) -> None:
        """
        Evolve the atmosphere over a given time step using the RK4 method.
        
        This method uses a 4th order Runge-Kutta integration to solve the 
        differential equations governing atmospheric evolution, providing better 
        accuracy for stiff systems with rapid escape rates.
        
        Args:
            t (float): Current time in years
            dt (float): Time step in years
        """
        if self.state == 'inflated':
            # Define the rate function for RK4
            def rate_function(t, mass):
                return (self._planet.outgassed_mass_flux(t) - 
                       self._planet.escape_mass_flux(t))
            
            # Calculate the four RK4 steps
            k1 = rate_function(t, self.gas_mass)
            k2 = rate_function(t + dt/2, self.gas_mass + dt/2 * k1)
            k3 = rate_function(t + dt/2, self.gas_mass + dt/2 * k2)
            k4 = rate_function(t + dt, self.gas_mass + dt * k3)
            
            # Update gas mass using RK4 formula
            self.gas_mass += dt/6 * (k1 + 2*k2 + 2*k3 + k4)
            
            # Handle negative mass case
            if self.gas_mass < 0:
                self.gas_mass = 0
                self.collapse()
        elif self.state == 'collapsed':
            # For collapsed state, we only have outgassing
            # Since this is a simple accumulation, we can use the average rate
            outgassing_rate = self._planet.outgassed_mass_flux(t)
            self.condensed_volatile_mass += dt * outgassing_rate
        
        self.update_state()
