import _plot_functions
import _functions
import matplotlib.pyplot as plt
import periodictable as pt
from periodictable import constants
import numpy as np
import re


# Parameters
thick_mm = 0.26  # mm
# _element = 'U'
_input = 'UO3'
_input_parsed = re.findall(r'([A-Z][a-z]*)(\d*)', _input)
print(_input_parsed)
formula = {'U': 1, 'O': 3}
ratio_array = []
_natural_mix = 'Y'
sample_density = 0.7875  #2  # 0.7875  # pt.elements.isotope(element).density  # g/cm3  https://en.wikipedia.org/wiki/Cadmium
_database = 'ENDF_VIII'
energy_max = 300  # max incident energy in eV
energy_min = 0  # min incident energy in eV
energy_sub = 100
sub_x = energy_sub * (energy_max - energy_min)  # subdivided new x-axis
_multi_element = 'N'
_energy_x_axis = 'Y'  # 1 means plot x-axis as energy in eV
_trans_y_axis = 'N'  # 1 means plot y-axis as transmission
_plot_each_iso_contribution = 'N'  # 1 means plot each isotope's contribution
_plot_each_ele_contribution = 'Y'  # 1 means plot each element's contribution
_plot_mixed = 'Y'  # 1 means plot mixed resonance

elements = list(dict.keys(formula))
ratios = list(dict.values(formula))
mass_iso_ele_dict = {}
y_i_iso_ele_dicts = {}
y_i_iso_ele_sum_dict = {}
df_raw_dict = {}
for _each_ in elements:
    _element = _each_
    ele_at_ratio = formula[_each_] / sum(ratios)
    # Get pre info (isotopes, abundance, mass, density) of each element from the formula
    isotopes, iso_abundance, iso_density, iso_mass, abundance_dict, density_dict, mass_dict, file_names = \
        _plot_functions.get_pre_data(_database, _element)

    mass_iso_ele_dict[_each_] = _plot_functions.get_atom_per_cm3(iso_abundance, iso_mass, ele_at_ratio,
                                                                 _natural_mix, ratio_array)

    x_energy, y_i_iso_ele_dict, y_i_iso_ele_sum, df_raw_dict[_each_] = \
        _plot_functions.get_xy(isotopes, file_names, energy_min, energy_max, iso_abundance,
                               sub_x, ele_at_ratio, _natural_mix, ratio_array)
    y_i_iso_ele_dicts[_each_] = y_i_iso_ele_dict #list(dict.values(y_i_iso_ele_dict))
    y_i_iso_ele_sum_dict[_each_] = y_i_iso_ele_sum

# Get Number of atoms per unit volume (#/cm^3)
mass_iso_ele_list = list(dict.values(mass_iso_ele_dict))
mass_iso_ele_sum = sum(np.array(mass_iso_ele_list))
mixed_atoms_per_cm3 = sample_density * pt.constants.avogadro_number/mass_iso_ele_sum
# Use function: mixed_atoms_per_cm3 = _functions.atoms_per_cm3(density=sample_density, mass=mass_iso_ele_sum)

print(y_i_iso_ele_dicts)
print(y_i_iso_ele_sum_dict)
# print(df_raw_dict)
keys = list(dict.keys(y_i_iso_ele_sum_dict))
yi_values = list(dict.values(y_i_iso_ele_sum_dict))
yi_values_sum = sum(yi_values)
print(y_i_iso_ele_sum_dict)

trans_sum = _functions.sig2trans_quick(thick_mm, mixed_atoms_per_cm3, yi_values_sum)
y_trans_tot = trans_sum
_all = elements

_plot_functions.plot_xy(elements, _energy_x_axis, _trans_y_axis, _plot_each_ele_contribution, _plot_mixed,
                        x_energy, y_trans_tot, thick_mm, mixed_atoms_per_cm3, y_i_iso_ele_sum_dict)
# _x_axis = x_energy
# _y_axis = 1 - trans_sum
# plt.plot(_x_axis, _y_axis) #label=_element+' natural mixture')
# plt.ylim(-0.01, 1.01)
# plt.show()
# _plot_functions.plot_xy(_element, _energy_x_axis, _trans_y_axis, _plot_each_contribution, _plot_mixed,
#             x_energy, y_trans_tot, isotopes, trans_dict)



# _plot_functions.plot_xy(_element, _energy_x_axis, _trans_y_axis, _plot_each_contribution, _plot_mixed,
#                             x_energy, y_trans_tot, isotopes, trans_dict)

# mass_iso_ele_list = list(dict.values(mass_iso_ele_dict))
# mass_iso_ele_array = np.array(mass_iso_ele_list)
# mixed_atoms_per_cm3 = sample_density * pt.constants.avogadro_number/sum(mass_iso_ele_array)
# print('Number of atoms per unit volume (#/cm^3): {}'.format(mixed_atoms_per_cm3))
