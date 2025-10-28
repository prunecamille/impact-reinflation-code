"""
Module for handling planetary properties and evolution.

This module provides the Planet class, which models planetary properties and 
evolution, particularly focusing on atmospheric evolution under stellar irradiation 
and impacts.
"""

from typing import Optional, List, Dict, Union
import numpy as np
from .constants import (
    G, AU_TO_M, SECONDS_PER_YEAR, R_EARTH, M_EARTH,
    SIGMA_BOLTZMANN, KAPPA_CO2, CURRENT_EARTH_OUTGASSING_RATE, 
    CO2_TRIPLE_POINT_PRESSURE
)
from .impactor import Impactor
from .atmosphere import Atmosphere
from scipy.optimize import brentq


class Planet:
    """
    A class representing a planet and its properties relevant for atmosphere evolution.
    
    This class models planetary properties and evolution, including atmospheric 
    evolution under stellar irradiation and impacts. It handles the interaction 
    between the planet's atmosphere, impacts, and stellar irradiation.
    
    Attributes:
        name (str): Name of the planet
        mass (float): Planet mass in kg
        radius (float): Planet radius in m
        star (Star): Associated star
        a (float): Semi-major axis in m
        outgassing_rate (float): Outgassing rate in kg/m²/yr
        escape_efficiency (float): Efficiency of atmospheric escape
        atmosphere (Optional[Atmosphere]): Planet's atmosphere
        impactor (Optional[Impactor]): Impactor object
    """
    
    def __init__(self, name: str, mass: float, radius: float, star: 'Star', 
                 a: float, outgassing_rate: float = 1e-2, escape_efficiency: float = 0.1):
        """
        Initialize a Planet instance.
        
        Args:
            name (str): Name of the planet
            mass (float): Planet mass in Earth masses
            radius (float): Planet radius in Earth radii
            star (Star): Associated star
            a (float): Semi-major axis in AU
            outgassing_rate (float, optional): Outgassing rate in kg/m²/yr. Defaults to 1e-2.
            escape_efficiency (float, optional): Escape efficiency. Defaults to 0.1.
            
        Raises:
            ValueError: If mass, radius are non-positive
        """
        if mass <= 0 or radius <= 0:
            raise ValueError("mass, radius must be positive")
            
        self.name = name
        self.mass = mass * M_EARTH  # kg
        self.radius = radius * R_EARTH  # m
        self.star = star
        self.a = a * AU_TO_M  # m
        self.outgassing_rate = outgassing_rate  # kg/m²/yr
        self.escape_efficiency = escape_efficiency
        
        # Derived properties
        self.surface_gravity = G * self.mass / self.radius**2  # m/s^2
        self.surface_area = 4 * np.pi * self.radius**2  # m^2
        self.escape_velocity = np.sqrt(2 * G * self.mass / self.radius)  # m/s
        self.grav_pot = -G * self.mass / self.radius  # m²/s²
        self.density = self.mass / (4/3 * np.pi * self.radius**3)  # kg/m^3
        
        # Initialize components
        self.atmosphere: Optional[Atmosphere] = None
        self.impactor: Optional[Impactor] = None
        
        # Assign planet to star
        self.star.assign_planet(self)
    
    def add_atmosphere(self, pressure0: Optional[float] = None, 
                      condensed_volatile_thickness0: float = 0.0, 
                      critical_pressure: Optional[float] = None, 
                      gas_mass0: Optional[float] = None) -> None:
        """
        Add an atmosphere to the planet.
        
        Args:
            pressure0 (Optional[float]): Initial pressure in bar. Either pressure0 
                or gas_mass0 must be provided.
            condensed_volatile_thickness0 (float): Initial condensed volatile 
                thickness in m. Defaults to 0.0.
            critical_pressure (Optional[float]): Critical pressure in bar. If None, 
                will be calculated.
            gas_mass0 (Optional[float]): Initial gas mass in kg. Either pressure0 
                or gas_mass0 must be provided.
            
        Raises:
            ValueError: If neither pressure0 nor gas_mass0 is provided, or if both are provided.
        """
        if (pressure0 is None and gas_mass0 is None) or \
           (pressure0 is not None and gas_mass0 is not None):
            raise ValueError("Exactly one of pressure0 or gas_mass0 must be provided")
        
        if pressure0 is not None:
            self.atmosphere = Atmosphere(self, pressure0, condensed_volatile_thickness0, 
                                       critical_pressure)
        else:
            self.atmosphere = Atmosphere.from_gas_mass(self, gas_mass0, 
                                                      condensed_volatile_thickness0, 
                                                      critical_pressure)
    
    def add_impactor(self, diameter: float, impact_rate: float, velocity: float, 
                    density: float = 3000) -> None:
        """
        Add an impactor to the planet.
        
        Args:
            diameter (float): Impactor diameter in km
            impact_rate (float): Impact rate in yr^-1
            velocity (float): Impact velocity in m/s
            density (float, optional): Impactor density in kg/m^3. Defaults to 3000.
        """
        self.impactor = Impactor(diameter, impact_rate, velocity, density)
    
    def add_default_impactor(self, impact_rate: float) -> float:
        """
        Add a default impactor sized to reinflate the atmosphere.
        
        This method calculates the minimum impactor size needed to volatilize
        enough mass to bring the atmosphere to the critical pressure.
        
        Args:
            impact_rate (float): Impact rate in yr^-1
            
        Returns:
            float: Impact diameter in km
        """
        if self.atmosphere is None:
            raise ValueError("No atmosphere found")
        
        if impact_rate == 0:
            # No impacts will occur
            self.impactor = None
            return 0
        
        # Calculate the mass of volatiles needed to reach critical pressure
        mvol = self.atmosphere.pressure_to_mass(self.atmosphere.critical_pressure)
        
        # Energy needed to volatilize this mass
        E_vol = mvol * (self.atmosphere.volatile_latent_heat + 
                       self.atmosphere.volatile_cp * 
                       self.atmosphere.volatile_sublimation_temp)
        
        # Impact energy needed (efficiency of 0.5)
        E_impact = 2 * E_vol
        
        # Impact mass needed
        m_impact = 2 * E_impact / self.escape_velocity**2
        
        # Impactor radius
        r_impact = (3 * m_impact / (4 * np.pi * 3000))**(1/3)
        r_impact = np.ceil(r_impact)
        
        # Impactor diameter (convert from m to km and round up)
        d_impact = 2 * np.ceil(r_impact * 1e-3)
        
        # Add the impactor
        self.add_impactor(d_impact, impact_rate, self.escape_velocity)
        
        return d_impact
    
    def escape_mass_flux(self, time: float) -> float:
        """
        Calculate the atmospheric escape mass flux at a given time.
        
        The escape flux is calculated using the stellar XUV emission and the
        escape efficiency.
        
        Args:
            time (float): Time in years
            
        Returns:
            float: Escape mass flux in kg/yr
        """
        if self.atmosphere is None:
            return 0.0
        
        # Escape flux equation
        dMdt = (self.escape_efficiency * np.pi * self.radius**3 * 
                self.star.XUV_emission(time) / (4 * np.pi * self.a**2 * G * self.mass))
        
        return dMdt * SECONDS_PER_YEAR  # kg/yr
    
    def outgassed_mass_flux(self, time: float) -> float:
        """
        Calculate the outgassed mass flux at a given time.
        
        Args:
            time (float): Time in years
            
        Returns:
            float: Outgassed mass flux in kg/yr
        """
        return self.outgassing_rate * self.surface_area  # kg/yr
    
    def impact(self, impactor: Impactor) -> None:
        """
        Process an impact event.
        
        An impact volatilizes condensed volatiles into gas using 50% of the
        impactor's kinetic energy.
        
        Args:
            impactor (Impactor): The impacting object
        """
        if self.atmosphere is None:
            return
        
        # Calculate mass volatilized using 50% of kinetic energy
        impact_energy = 0.5 * impactor.kinetic_energy
        mass_volatilized = impact_energy / (
            self.atmosphere.volatile_latent_heat + 
            self.atmosphere.volatile_cp * 
            self.atmosphere.volatile_sublimation_temp
        )
        
        self.atmosphere.process_impact(mass_volatilized)
    
    def calculate_co2_condensation_temp(self, pressure: float) -> float:
        """
        Calculate the CO2 condensation temperature as a function of pressure.
        
        Based on Wordsworth et al. (2010b) CO2 condensation data.
        
        Args:
            pressure (float): Pressure in Pa
            
        Returns:
            float: Condensation temperature in K
        """
        if pressure < CO2_TRIPLE_POINT_PRESSURE:
            # Below triple point
            return -3167.8 / (np.log(0.01 * pressure) - 23.23)
        else:
            # Above triple point
            log_p = np.log(pressure)
            return 684.2 - 92.3 * log_p + 4.32 * log_p**2
    
    def calculate_thin_radiator_temp(self, pressure: float, 
                                     kappa: float = KAPPA_CO2, 
                                     albedo: float = 0.2) -> float:
        """
        Calculate thin radiator temperature for a given pressure.
        
        Args:
            pressure (float): Pressure in Pa
            kappa (float): CO2 opacity in m²/kg. Defaults to KAPPA_CO2.
            albedo (float): Planet albedo. Defaults to 0.2.
            
        Returns:
            float: Temperature in K
        """
        F = self.star.F(self)
        F_eff = F * (1 - albedo)  # Effective flux accounting for albedo
        Ttr = (F_eff * kappa * pressure / 
              (4 * SIGMA_BOLTZMANN * self.surface_gravity))**0.25
        return Ttr
    
    @property
    def equilibrium_temperature(self) -> float:
        """
        Calculate the equilibrium temperature using energy balance.
        
        Uses the simple energy balance: T_eq = (F*(1-A)/(4*sigma))^0.25
        where F is stellar flux, A is albedo, and sigma is Stefan-Boltzmann constant.
        
        Returns:
            float: Equilibrium temperature in K
        """
        albedo = 0.2  # Fixed albedo
        F = self.star.F(self)  # Stellar flux at planet's distance
        F_eff = F * (1 - albedo)  # Effective flux accounting for albedo
        
        # Energy balance: F_eff = 4 * sigma * T_eq^4
        T_eq = (F_eff / (4 * SIGMA_BOLTZMANN))**0.25
        
        return T_eq
    
    def evolve_system(self, age: float, dt: float = 1e6, 
                     method: str = 'euler') -> Dict[str, Union[np.ndarray, float]]:
        """
        Evolve the system forward in time.
        
        This method evolves the planet's atmosphere over time, accounting for:
        - Atmospheric escape due to XUV irradiation
        - Outgassing from the planetary interior
        - Impact events that volatilize condensed gases
        
        Args:
            age (float): Time to evolve for in years
            dt (float): Time step in years
            method (str): Integration method ('euler' or 'rk4'). Defaults to 'euler'.
            
        Returns:
            Dict containing evolution history:
            - times: Array of times in years
            - atm_masses: Array of atmospheric gas masses in kg
            - vol_masses: Array of condensed volatile masses in kg
            - pressures: Array of pressures in bar
            - states: Array of states (0 for collapsed, 1 for inflated)
            - impact_times: Array of times when impacts occur in years
            - inflation_percentage: Percentage of time spent in inflated state
            - num_impacts: Number of impacts
        """
        if self.atmosphere is None:
            print("Warning: Starting with no atmosphere")
            self.add_atmosphere(pressure0=0, condensed_volatile_thickness0=0)
        
        self.atmosphere.reset()

        # Time array
        times = np.arange(0, age, dt)

        # Initialize variables
        pressures = [self.atmosphere.pressure / 1e5]  # Convert to bar
        atm_masses = [self.atmosphere.gas_mass]
        vol_masses = [self.atmosphere.condensed_volatile_mass]
        impact_times = []
        states = [1 if self.atmosphere.state == 'inflated' else 0]
        
        # Initialize next impact time
        next_impact_time = self.impactor.time_to_next_impact() if self.impactor else np.inf

        # Evolve system
        for t in times:
            # Evolve atmosphere using specified method
            if method == 'euler':
                self.atmosphere.evolve_euler(t, dt)
            elif method == 'rk4':
                self.atmosphere.evolve_rk4(t, dt)
            else:
                print(f"Invalid method: {method}, using euler")
                self.atmosphere.evolve_euler(t, dt)
            
            # Store state
            atm_masses.append(self.atmosphere.gas_mass)
            vol_masses.append(self.atmosphere.condensed_volatile_mass)
            pressures.append(self.atmosphere.pressure / 1e5)
            states.append(1 if self.atmosphere.state == 'inflated' else 0)

            # Process impacts
            if t >= next_impact_time:
                self.impact(self.impactor)
                impact_times.append(t)
                next_impact_time = t + self.impactor.time_to_next_impact()

        num_impacts = len(impact_times)
        # Calculate inflation percentage
        inflation_percentage = 100 * np.sum(states) / len(states)

        # Add t=0 to the times array
        times = np.insert(times, 0, 0)
        
        return {
            'times': np.array(times),
            'atm_masses': np.array(atm_masses),
            'vol_masses': np.array(vol_masses),
            'pressures': np.array(pressures),
            'states': np.array(states),
            'impact_times': np.array(impact_times),
            'inflation_percentage': inflation_percentage,
            'num_impacts': num_impacts
        }
    
    def solve_for_critical_pressure(self) -> float:
        """
        Solve for the critical pressure where the condensation temperature of CO2 
        equals the equilibrium temperature of a thin CO2 atmosphere.
        
        The critical pressure is the threshold below which the atmosphere collapses
        from an inflated state to a collapsed state.
        
        Returns:
            float: Critical pressure in Pa
        """
        def difference(p):
            return (self.calculate_co2_condensation_temp(p) - 
                   self.calculate_thin_radiator_temp(p))

        # Choose pressure bounds (in Pascals)
        p_min = 0.5e2  # 0.0005 bar
        p_max = 1e8  # 1000 bar
        
        try:
            p_crit = brentq(difference, p_min, p_max)
            return p_crit
        except ValueError:
            # Depending on which boundary it hits, return the minimum or maximum pressure
            if difference(p_min) < difference(p_max):
                return p_max
            else:
                return p_min
