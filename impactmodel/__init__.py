"""
Impact-reinflation: A Python package for modeling planetary atmosphere evolution 
through impacts, atmospheric erosion, and outgassing.

This package simulates atmospheric evolution on exoplanets, tracking the 
interaction between stellar XUV irradiation, impact-driven volatile delivery, 
and internal outgassing processes.
"""

from .planet import Planet
from .star import Star
from .impactor import Impactor
from .atmosphere import Atmosphere
from .constants import *

__version__ = '1.0.0'
__all__ = ['Planet', 'Star', 'Impactor', 'Atmosphere']
