import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
import _utilities
from resonance import Resonance

# Global parameters
_energy_min = 10
_energy_max = 400
_energy_step = 0.01
# Input sample name or names as str, case sensitive
_layer_1 = 'U6O21Gd2'
_name = _layer_1
_thickness_1 = 3 # mm
_density_1 = 0.7875 # g/cm3 deviated due to porosity

o_reso = Resonance(energy_min=_energy_min, energy_max=_energy_max, energy_step=_energy_step)
o_reso.add_layer(formula=_layer_1, thickness=_thickness_1, density=_density_1)

o_reso.plot(all_elements=True, transmission=True, x_axis='time')
