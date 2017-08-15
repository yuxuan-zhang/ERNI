import pprint
from resonance import Resonance
import _utilities

# Global parameters
_energy_min = 1e-5
_energy_max = 300
_energy_step = 0.01
# Input sample name or names as str, case sensitive
_layer_1 = 'U3O8'
_thickness_1 = 0.5 # mm
_density_1 = 2 # g/cm3 deviated due to porosity

_layer_2 = 'Gd'
_thickness_2 = 0.075 # mm

o_reso = Resonance(energy_min=_energy_min, energy_max=_energy_max, energy_step=_energy_step)
o_reso.add_layer(formula=_layer_1, thickness=_thickness_1, density=_density_1)
o_reso.add_layer(formula=_layer_2, thickness=_thickness_2)
# new_list_ratio = [0., 0., 0.15, 0.85]
# o_reso.set_isotopic_ratio(compound='UO3', element='U', list_ratio=new_list_ratio)
# pprint.pprint(o_reso)
o_reso.plot(mixed=True, all_layers=True, transmission=False)
