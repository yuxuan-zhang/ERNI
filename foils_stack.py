import _plot_functions
import _functions
import matplotlib.pyplot as plt
import periodictable as pt
from periodictable import constants
import numpy as np
import pandas as pd


foils_thickmm_dict = {'Gd': 1,
                      'Cd': 1,
                      'U': 1,
                      'Ag': 1,
                      'Au': 1,
                      'Ta': 1,
                      'In': 1,
                      'W': 1,
                      'Co': 1,
                      'Hf': 1,
                      'Pb': 1,
                      'Eu': 1,
                      'B': 1
                      }
elements = list(dict.keys(foils_thickmm_dict))
thick = list(dict.values(foils_thickmm_dict))
elements_str = ''.join(elements)
print(elements_str)
# Parameters
_input = elements_str
_natural_ele = ['Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y']
print(len(_natural_ele) == len(elements))
thick_mm = 5  # mm
_input_density = 0.7875  # g/cm3  not needed if _input is single element
_input_ratios_dict = {'U': [0, 0, 0.15, 0.85]
                      # 'O': [1, 0, 0]
                      }  # at.%  not needed if _natural_ele is 'Y'
_database = 'ENDF_VIII'
energy_max = 21  # max incident energy in eV
energy_min = 0  # min incident energy in eV
energy_sub = 1000
sub_x = energy_sub * (energy_max - energy_min)  # subdivided new x-axis
# _multi_element = 'N'
_energy_x_axis = 'Y'  # 1 means plot x-axis as energy in eV
_trans_y_axis = 'N'  # 1 means plot y-axis as transmission
_plot_each_ele_contribution = 'Y'  # 1 means plot each element's contribution
_plot_each_iso_contribution = 'N'  # 1 means plot each isotope's contribution
_plot_mixed = 'N'  # 1 means plot mixed resonance

formula = _functions.input2formula(_input)  # Function called to parse input formula and return elements and ratios
elements = list(dict.keys(formula))
ratios = list(dict.values(formula))
natural_ele_dict, ratios_dict = _functions.formula_ratio_array(elements, _natural_ele, _input_ratios_dict)

unnatrual_ratio_array_dict = ratios_dict

if len(elements) == 1:
    sample_density = pt.elements.isotope(_input).density  # g/cm3  https://en.wikipedia.org/wiki/Cadmium
else:
    # _input_density = 0.7875  # g/cm3  need to input while the _input is multi-element mixture
    sample_density = _input_density

# if _natural_ele == 'Y':
#     for _each_ in elements:
#         unnatrual_ratio_array_dict[_each_] = []
# else:
#     _p = 0
#     for _each_ in elements:
#         unnatrual_ratio_array_dict[_each_] = [_p]
#         _p = _p + 1
# print('Unnatrual_ratio_array_dict: ', unnatrual_ratio_array_dict)


mass_iso_ele_dict = {}  # For number of atoms per cm3 calculation
y_i_iso_ele_dicts = {}  # For transmission calculation at isotope lever
y_i_iso_ele_sum_dict = {}  # For transmission calculation at element lever
df_raw_dict = {}  # Raw sigma data for elements and isotopes
isotopes_dict = {}  # List all isotopes for each element involved
abundance_dicts = {}  # List all natrual abundance for each isotope of each element involved
for _each_ in elements:
    _element = _each_
    ele_at_ratio = formula[_each_] / sum(ratios)
    # Get pre info (isotopes, abundance, mass, density) of each element from the formula
    isotopes_dict[_each_], iso_abundance, iso_density, iso_mass, abundance_dict, density_dict, mass_dict, file_names = \
        _plot_functions.get_pre_data(_database, _element)

    mass_iso_ele_dict[_each_] = _plot_functions.get_mass_iso_ele(iso_abundance, iso_mass, ele_at_ratio,
                                                                 natural_ele_dict[_each_],
                                                                 unnatrual_ratio_array_dict[_each_])

    x_energy, y_i_iso_ele_dict, y_i_iso_ele_sum, df_raw_dict[_each_] = \
        _plot_functions.get_xy(isotopes_dict[_each_], file_names, energy_min, energy_max, iso_abundance,
                               sub_x, ele_at_ratio, natural_ele_dict[_each_], unnatrual_ratio_array_dict[_each_])
    y_i_iso_ele_dicts[_each_] = y_i_iso_ele_dict
    y_i_iso_ele_sum_dict[_each_] = y_i_iso_ele_sum
    if natural_ele_dict[_each_] == 'N':
        abundance_dicts[_each_] = unnatrual_ratio_array_dict[_each_]
    else:
        abundance_dicts[_each_] = abundance_dict

print('Abundance_dicts: ', abundance_dicts)
# print(mass_iso_ele_dict)
# print(ele_at_ratio)

mass_iso_ele_list = list(dict.values(mass_iso_ele_dict))
mass_iso_ele_sum = sum(np.array(mass_iso_ele_list))
mixed_atoms_per_cm3 = sample_density * pt.constants.avogadro_number/mass_iso_ele_sum
# Use function: mixed_atoms_per_cm3 = _functions.atoms_per_cm3(density=sample_density, mass=mass_iso_ele_sum)

keys = list(dict.keys(y_i_iso_ele_sum_dict))
yi_values = list(dict.values(y_i_iso_ele_sum_dict))
yi_values_sum = sum(yi_values)
trans_sum = _functions.sig2trans_quick(thick_mm, mixed_atoms_per_cm3, yi_values_sum)
y_trans_tot = trans_sum
# print(y_i_iso_ele_sum_dict)


### Create the trans or absorb dict of ele for plotting if needed
if _plot_each_ele_contribution == 'Y':
    y_ele_dict = {}
    for _ele in elements:
        if _trans_y_axis == 'Y':
            y_ele_dict[_ele] = _functions.sig2trans_quick(thick_mm, mixed_atoms_per_cm3, y_i_iso_ele_sum_dict[_ele])
        else:
            y_ele_dict[_ele] = 1 - _functions.sig2trans_quick(thick_mm, mixed_atoms_per_cm3, y_i_iso_ele_sum_dict[_ele])

### Create the trans or absorb dict : y_iso_dicts of isotopes for plotting if needed
if _plot_each_iso_contribution == 'Y':
    y_iso_dicts = {}
    y_iso_dict = {}
    for _ele in elements:
        for _iso in isotopes_dict[_ele]:
            if _trans_y_axis == 'Y':
                y_iso_dict[_iso] = _functions.sig2trans_quick(thick_mm, mixed_atoms_per_cm3,
                                                              y_i_iso_ele_dicts[_ele][_iso])
            else:
                y_iso_dict[_iso] = 1 - _functions.sig2trans_quick(thick_mm, mixed_atoms_per_cm3,
                                                                  y_i_iso_ele_dicts[_ele][_iso])
        y_iso_dicts[_ele] = y_iso_dict
        y_iso_dict = {}  # Clear for following set of isotopes
    # print(y_iso_dicts)




### Determine x y axis types and captions
if _energy_x_axis == 'Y':
    _x_axis = x_energy
    _x_words = 'Energy (eV)'
else:
    _x_axis = _functions.ev2lamda(x_energy)
    _x_words = 'Wavelength (Å)'

if _trans_y_axis == 'Y':
    _y_words = 'Neutron transmission'
else:
    _y_words = 'Neutron attenuation'

### Determine x y axis values
if _plot_mixed == 'Y':
    if _trans_y_axis == 'Y':
        _y_axis = y_trans_tot
    else:
        _y_axis = 1 - y_trans_tot
    plt.plot(_x_axis, _y_axis, label=_input)

if _plot_each_ele_contribution == 'Y':
    for _ele in elements:
        _y_each_axis = y_ele_dict[_ele]
        plt.plot(_x_axis, _y_each_axis, label=_ele)

if _plot_each_iso_contribution == 'Y':
    for _ele in elements:
        for _iso in isotopes_dict[_ele]:
            _y_each_axis = y_iso_dicts[_ele][_iso]
            plt.plot(_x_axis, _y_each_axis, label=_iso)

plt.ylim(-0.01, 1.01)
plt.xlabel(_x_words)
plt.ylabel(_y_words)
plt.legend(loc='best')
plt.show()
