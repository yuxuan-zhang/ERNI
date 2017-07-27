import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import _functions
import _fit_functions
import _plot_functions
from lmfit import minimize, Parameters
import os
import periodictable as pt
from periodictable.constants import avogadro_number
import peakutils as pku

# Parameters
source_to_detector_cm = 1610.569047826142  # cm
delay_ms = -12.112589054625005  # ms
delay_us = delay_ms * 1000
_slice = 220
_database = 'ENDF_VIII'
_input_formula = 'CoAg'
energy_min = 0
energy_max = 400
energy_sub = 100
sub_x = energy_sub * (energy_max - energy_min)  # subdivided new x-axis
_input_thick_mm = 0.025
_input_thick_cm = _input_thick_mm/10
time_lamda_ev_axis = 'eV'
_name = 'foil6'
data_path = 'data/' + _name + '.csv'
spectra_path = 'data/spectra.txt'


'''Describe your sample: '''
stacked_foil_boo = 'Y'  # Single element foil or stacked foils: Y/N?
# Thickness input:
special_thick_boo = 'N'
special_thick_element_str = str
special_thick_mm_list = []
special_thick_cm_list = np.array(special_thick_mm_list)/10
# Enriched isotope ratio input:
enrichment_boo = 'N'  # Isotopic enreiched or depleted: Y/N?
enriched_element_str = str
input_ratio_dict = {}
# Special density input:
special_density_boo = 'N'
special_density_element_str = str
special_density_gcm3_list = []


'''How you want the data to be plotted?'''
_plot_or_not = 'Y'
_energy_x_axis = 'Y'  # 1 means plot x-axis as energy in eV
_trans_y_axis = 'N'  # 1 means plot y-axis as transmission
_plot_each_ele_contribution = 'Y'  # 1 means plot each element's contribution
_plot_each_iso_contribution = 'N'  # 1 means plot each isotope's contribution
_plot_mixed = 'N'  # 1 means plot mixed resonance

formula_dict = _functions.input2formula(_input_formula)
elements = _functions.dict_key_list(formula_dict)
ratios = _functions.dict_value_list(formula_dict)
sum_ratios = sum(ratios)
isotope_dict = _functions.get_isotope_dicts(_database, elements)

# DICT 1: Thickness dict with option for modification
thick_cm_dict = _functions.repeat_value_dict(elements, _input_thick_cm)
if stacked_foil_boo == 'Y':
    if special_thick_boo == 'Y':
        thick_cm_dict = _plot_functions.modify_thick_cm_dict_by_input(thick_cm_dict, special_thick_element_str, special_thick_cm_list)

# DICT 2: Isotopic mass dict
iso_mass_dicts = _functions.get_iso_mass_dicts_quick(elements, isotope_dict)

# Dict 3: Molar mass dict
molar_mass_dict = _functions.get_molar_mass_dict(elements)

# DICT 4: Isotope at.% dict with option for modification
iso_ratio_dicts = _functions.get_iso_ratio_dicts_quick(elements, isotope_dict)

# DICT 5: Density dict
density_gcm3_dict = _functions.get_density_dict(elements)


# Update DICT 3 & 4 & 5: isotopic ratio changes lead to |Density| & |Molar mass| changes
if enrichment_boo == 'Y':
    # Update isotope at.% ratio dict
    iso_ratio_dicts, enriched_element = _plot_functions.modify_iso_ratio_dicts(elements, isotope_dict, enriched_element_str, input_ratio_dict)
    # Update molar mass dict
    molar_mass_dict = _plot_functions.modify_molar_mass_dict_by_enrichment(molar_mass_dict, enriched_element, isotope_dict, iso_ratio_dicts, iso_mass_dicts)
    # Update density dict
    density_gcm3_dict = _plot_functions.modify_density_dict_by_enrichment(density_gcm3_dict, enriched_element, isotope_dict, iso_ratio_dicts)

# Update DICT 5: Density dict, if special case considered
if stacked_foil_boo == 'Y':
    if special_density_boo == 'Y':
        # Stacked foils and would like to modify density for specific element
        density_gcm3_dict = _plot_functions.modify_density_dict_by_input(density_gcm3_dict, special_density_element_str, special_density_gcm3_list)
else:
    if special_density_boo == 'Y':
        # Not isolated elements or mixture or compound need density input currently
        input_tot_density = 0.7875


'''For plotting the database'''
avo_divi_mass_iso_ele_dict = {}  # For number of atoms per cm3 calculation
sigma_iso_ele_eleisodict = {}  # For transmission calculation at isotope level
sigma_iso_ele_sum_eledict = {}  # For transmission calculation at element level
sigma_iso_ele_sum_l_eledict = {}
sigma_iso_ele_l_eleisodict = {}
df_raw_dict = {}  # Raw sigma data for elements and isotopes
sample_density_dict = {}
atoms_per_cm3_dict = {}

for el in elements:
    # isotopes_list = list(dict.keys(iso_ratio_dicts[el]))
    iso_ratio_list = list(dict.values(iso_ratio_dicts[el]))
    iso_ratio_array = np.array(iso_ratio_list)
    iso_mass_list = list(dict.values(iso_mass_dicts[el]))
    iso_mass_array = np.array(iso_mass_list)
    ele_at_ratio = formula_dict[el] / sum_ratios

    # A part for getting atoms_per_cm3, this part is irrelevant to fitting parameters, and will be exported for fitting
    avo_divi_mass_iso_ele_dict[el] = avogadro_number / (sum(iso_ratio_array * iso_mass_array) * ele_at_ratio)
    if stacked_foil_boo == 'Y':
        # Multiple foils stacked
        # sample_density_dict[el] = density_gcm3_dict[el] * 1
        atoms_per_cm3_dict[el] = density_gcm3_dict[el] * avo_divi_mass_iso_ele_dict[el]
    else:
        if special_density_boo == 'Y':
            # Total density of mixture
            tot_density = input_tot_density

    # Get sigma related terms
    file_names = _functions.get_file_path(_database, el)
    x_energy, sigma_iso_ele_isodict, sigma_iso_ele_sum, df_raw_dict[el] \
        = _plot_functions.get_xy_from_database(iso_ratio_dicts[el],
                                               file_names,
                                               energy_min,
                                               energy_max,
                                               iso_ratio_list,
                                               sub_x,
                                               ele_at_ratio)
    # Two level dict of isotopic array of (L * sigma * iso_ratio * ele_ratio)
    # sigma_iso_ele_l_eleisodict[el] = sigma_iso_ele_l_isodict
    # One level dict of elemental array of (L * sigma * iso_ratio * ele_ratio)
    # sigma_iso_ele_sum_l_eledict[el] = sigma_iso_ele_sum * thick_cm_dict[el]

    # Two level dict of isotopic array of (sigma * iso_ratio * ele_ratio)
    sigma_iso_ele_eleisodict[el] = sigma_iso_ele_isodict
    # One level dict of elemental array of (sigma * iso_ratio * ele_ratio)
    sigma_iso_ele_sum_eledict[el] = sigma_iso_ele_sum


if stacked_foil_boo == 'Y':
    # Stacked foils
    mixed_l_n_avo = _plot_functions.l_x_n_multi_ele_stack(elements, thick_cm_dict, density_gcm3_dict, molar_mass_dict)
else:
    thick_cm_list = list(dict.values(thick_cm_dict))
    thick_cm = thick_cm_list[0]
    sample_density = tot_density
    avo_divi_mass_iso_ele_list = list(dict.values(avo_divi_mass_iso_ele_dict))
    avo_divi_mass_iso_ele_sum = sum(np.array(avo_divi_mass_iso_ele_list))
    mixed_atoms_per_cm3 = sample_density * avo_divi_mass_iso_ele_sum
    mixed_l_n_avo = thick_cm * mixed_atoms_per_cm3

# sum of (sigma * ele_ratio * iso_ratio * l)
yi_values_l = list(dict.values(sigma_iso_ele_sum_l_eledict))
yi_values_l_sum = sum(yi_values_l)
# sum of (sigma * ele_ratio * iso_ratio)
yi_values = list(dict.values(sigma_iso_ele_sum_eledict))
yi_values_sum = sum(yi_values)

trans_sum = _functions.sig_l_2trans_quick(mixed_l_n_avo, yi_values_sum)
y_trans_tot = trans_sum
y_attenu_tot = 1 - trans_sum

# Create the trans or absorb dict of ele for plotting if needed
if _plot_each_ele_contribution == 'Y':
    y_ele_dict = {}
    for _ele in elements:
        if _trans_y_axis == 'Y':
            y_ele_dict[_ele] = _functions.sig_l_2trans_quick(mixed_l_n_avo, sigma_iso_ele_sum_eledict[_ele])
        else:
            y_ele_dict[_ele] = 1 - _functions.sig_l_2trans_quick(mixed_l_n_avo, sigma_iso_ele_sum_eledict[_ele])

# Create the trans or absorb dict : y_iso_dicts of isotopes for plotting if needed
y_iso_dicts = {}
if _plot_each_iso_contribution == 'Y':
    for _ele in elements:
        y_iso_dict = {}
        for _iso in isotope_dict[_ele]:
            if _trans_y_axis == 'Y':
                y_iso_dict[_iso] = _functions.sig_l_2trans_quick(mixed_l_n_avo,
                                                                 sigma_iso_ele_eleisodict[_ele][_iso])
            else:
                y_iso_dict[_iso] = 1 - _functions.sig_l_2trans_quick(mixed_l_n_avo,
                                                                     sigma_iso_ele_eleisodict[_ele][_iso])
        y_iso_dicts[_ele] = y_iso_dict

# Plot the theoretical neutron resonance
if _plot_or_not == 'Y':
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
        plt.plot(_x_axis, _y_axis, label=_input_formula)

    if _plot_each_ele_contribution == 'Y':
        for _ele in elements:
            _y_each_axis = y_ele_dict[_ele]
            plt.plot(_x_axis, _y_each_axis, label=_ele)

    if _plot_each_iso_contribution == 'Y':
        for _ele in elements:
            for _iso in isotope_dict[_ele]:
                _y_each_axis = y_iso_dicts[_ele][_iso]
                plt.plot(_x_axis, _y_each_axis, label=_iso)


x_data_array = _functions.get_spectra_slice(spectra_path, time_lamda_ev_axis, delay_us,
                                                source_to_detector_cm, _slice)
print(x_data_array)
y_data_array = 1 - _functions.get_normalized_data_slice('data/'+_name+'.csv', _slice)/4.2
print(y_data_array)
plt.plot(x_data_array, y_data_array, 'r-', label=_name, alpha=0.7)

indexes = pku.indexes(y_data_array, thres=0.6, min_dist=50)

print(indexes)
peaks_x = pku.interpolate(x_data_array, y_data_array, ind=indexes)
print(peaks_x)
# plt.figure(figsize=(10, 6))
plt.plot(x_data_array[indexes], y_data_array[indexes], 'bx', label='peak')
plt.title('First estimation')

# paramas = Parameters()
# paramas.add('source_to_detector_cm', value=1610.9)
# paramas.add('delay_ms', value=4.5-16.6)
# paramas.add('density_gcm3', value=pt.elements.isotope(_input_element).density)

plt.ylim(-0.01, 1.01)
plt.xlim(energy_min, energy_max)
plt.xlabel(_x_words)
plt.ylabel(_y_words)
plt.legend(loc='best')
plt.show()