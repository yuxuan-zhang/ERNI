import _plot_functions
import _functions
import matplotlib.pyplot as plt
import periodictable as pt
from periodictable import constants
import numpy as np


# Parameters
thick_mm = 0.26  # mm
# _element = 'U'
formula = {'U': 1, 'O': 3}
ratio_array = []
_natural_mix = 'Y'
sample_density = 2  # 0.7875  # pt.elements.isotope(element).density  # g/cm3  https://en.wikipedia.org/wiki/Cadmium
# ele_at_ratio = 0.25  # for single element, will be implanted for multiple elements compound
# element2 = '16-O'
_database = 'ENDF_VIII'
energy_max = 300  # max incident energy in eV
energy_min = 1  # min incident energy in eV
energy_sub = 100
sub_x = energy_sub * (energy_max - energy_min)  # subdivided new x-axis
_multi_element = 'N'
_energy_x_axis = 'Y'  # 1 means plot x-axis as energy in eV
_trans_y_axis = 'N'  # 1 means plot y-axis as transmission
_plot_each_iso_contribution = 'N'  # 1 means plot each isotope's contribution
_plot_each_ele_contribution = 'Y'  # 1 means plot each element's contribution
_plot_mixed = 'N'  # 1 means plot mixed resonance

elements = list(dict.keys(formula))
numbers = list(dict.values(formula))
mass_iso_ele_dict = {}
y_i_iso_ele_dicts = {}
y_i_iso_ele_sum_dict = {}
df_raw_dict = {}
for _each_ in elements:
    _element = _each_
    ele_at_ratio = formula[_each_] / sum(numbers)
    # Get pre info (isotopes, abundance, mass, density) of each element from the formula
    isotopes, iso_abundance, iso_density, iso_mass, abundance_dict, density_dict, mass_dict, file_names = \
        _plot_functions.get_pre_data(_database, _element)

    mass_iso_ele_dict[_each_] = _plot_functions.get_atom_per_cm3(iso_abundance, iso_mass, ele_at_ratio,
                                                                 _natural_mix, ratio_array)

    x_energy, y_i_iso_ele_dict, y_i_iso_ele_sum, df_raw = _plot_functions.get_xy(isotopes, file_names, energy_min,
                                                                                 energy_max, iso_abundance, sub_x, ele_at_ratio)
    y_i_iso_ele_dicts[_each_] = list(dict.values(y_i_iso_ele_dict))
    y_i_iso_ele_sum_dict[_each_] = y_i_iso_ele_sum

print(y_i_iso_ele_dicts)
print(y_i_iso_ele_sum_dict)
mass_iso_ele_list = list(dict.values(mass_iso_ele_dict))
mass_iso_ele_sum = sum(np.array(mass_iso_ele_list))

print(mass_iso_ele_sum)
mixed_atoms_per_cm3 = _functions.atoms_per_cm3(density=sample_density, mass=mass_iso_ele_sum)
# print('Number of atoms per unit volume (#/cm^3): {}'.format(mixed_atoms_per_cm3))


# _plot_functions.plot_xy(_element, _energy_x_axis, _trans_y_axis, _plot_each_contribution, _plot_mixed,
#                             x_energy, y_trans_tot, isotopes, trans_dict)

# mass_iso_ele_list = list(dict.values(mass_iso_ele_dict))
# mass_iso_ele_array = np.array(mass_iso_ele_list)
# mixed_atoms_per_cm3 = sample_density * pt.constants.avogadro_number/sum(mass_iso_ele_array)
# print('Number of atoms per unit volume (#/cm^3): {}'.format(mixed_atoms_per_cm3))
