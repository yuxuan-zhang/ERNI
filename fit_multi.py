import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import _functions
import glob
import _fit_funtions
from lmfit import minimize, Parameters
import os
import periodictable as pt
from periodictable import constants
import peakutils as pku

# Parameters
source_to_detector_cm = 1610.9  # cm
delay_ms = 4.5 - 16.6127  # ms
delay_us = delay_ms * 1000
_slice = 220
_database = 'ENDF_VIII'
_input_element = 'Co'
energy_min = 0
energy_max = 400
energy_sub = 100
sub_x = energy_sub * (energy_max - energy_min)  # subdivided new x-axis
_thick_mm = 0.025
_thick_cm = _thick_mm/10
time_lamda_ev_axis = 'eV'
_name = 'foil6'
data_path = 'data/' + _name + '.csv'
spectra_path = 'data/spectra.txt'
# _fit_funtions.get_multi_data(file_name_signature, time_lamda_ev_axis, delay_us, source_to_detector_cm, _slice)


x_data_array = _functions.get_spectra_slice(spectra_path, time_lamda_ev_axis, delay_us,
                                                source_to_detector_cm, _slice)
print(x_data_array)
y_data_array = 1 - _functions.get_normalized_data_slice('data/'+_name+'.csv', _slice)/4.2
print(y_data_array)
plt.plot(x_data_array, y_data_array, 'r.', label=_name)

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

formula_dict = _functions.input2formula(_input_element)
elements = _functions.dict_key_list(formula_dict)
ratios = _functions.dict_value_list(formula_dict)

# To check whether the input are foils stacked
foils_stacked = ratios.count(ratios[0]) == len(ratios)
mass_iso_ele_dict = {}  # For number of atoms per cm3 calculation
sigma_iso_ele_eleisodict = {}  # For transmission calculation at isotope lever
sigma_iso_ele_sum_eledict = {}  # For transmission calculation at element lever
sigma_iso_ele_sum_l_eledict = {}
sigma_iso_ele_l_eleisodict = {}
df_raw_dict = {}  # Raw sigma data for elements and isotopes
isotopes_dict = {}  # List all isotopes for each element involved
abundance_dicts = {}  # List all natural abundance for each isotope of each element involved
all_ratio_dicts = {}
sample_density = .0
for _each_ in elements:
    _element = _each_
    if foils_stacked is False:
        if _mixture_or_ele_with_diff_density == 'Y':
            ele_at_ratio = formula_dict[_each_] / sum(ratios)
    else:
        if _mixture_or_ele_with_diff_density == 'Y':
            ele_at_ratio = formula_dict[_each_] / sum(ratios)
        else:
            ele_at_ratio = 1
    # Get pre info (isotopes, abundance, mass, density) of each element from the formula
    isotopes_dict[_each_], all_ratio_dicts[_each_], iso_abundance, iso_mass, file_names \
        = _plot_functions.get_pre_data(_database, _element)

    # Replace the at.% if isotope composition does not follow natural abundance in the main isotope ratio dict
    if _unnatural_ele_input == 'Y':
        for _ele_need_to_replace in unnatural_element:
            all_ratio_dicts[_ele_need_to_replace] = unnatural_ratio_dicts[_ele_need_to_replace]

    # A part for getting atoms_per_cm3, this part is irrelevant to fitting parameters, and will be exported for fitting
    mass_iso_ele_dict[_each_] = _plot_functions.get_mass_iso_ele(iso_abundance,
                                                                 iso_mass,
                                                                 ele_at_ratio,
                                                                 all_ele_boo_dict[_each_],
                                                                 all_ratio_dicts[_each_])
    # Total density calculation of mixture mixed by ele_at_ratio
    sample_density = sample_density + density_gcm3_dict[_each_] * ele_at_ratio

    thick_cm = thick_cm_dict[_each_]
    x_energy, sigma_iso_ele_isodict, sigma_iso_ele_l_isodict, sigma_iso_ele_sum, df_raw_dict[_each_] \
        = _plot_functions.get_xy(isotopes_dict[_each_],
                                 thick_cm,
                                 file_names,
                                 energy_min,
                                 energy_max,
                                 iso_abundance,
                                 sub_x,
                                 ele_at_ratio,
                                 all_ele_boo_dict[_each_],
                                 all_ratio_dicts[_each_])
    # Two level dict of isotopic array of (L * sigma * iso_ratio * ele_ratio)
    sigma_iso_ele_l_eleisodict[_each_] = sigma_iso_ele_l_isodict
    # One level dict of elemental array of (L * sigma * iso_ratio * ele_ratio)
    sigma_iso_ele_sum_l_eledict[_each_] = sigma_iso_ele_sum * thick_cm

    # Two level dict of isotopic array of (sigma * iso_ratio * ele_ratio)
    sigma_iso_ele_eleisodict[_each_] = sigma_iso_ele_isodict
    # One level dict of elemental array of (sigma * iso_ratio * ele_ratio)
    sigma_iso_ele_sum_eledict[_each_] = sigma_iso_ele_sum

print('Abundance_dicts: ', all_ratio_dicts)

if _mixture_or_ele_with_diff_density == 'Y':
    if _compound == 'Y':
        sample_density = float(input('Mixture or compound density of {} in g/cm3: '.format(_input_formula)))

mass_iso_ele_list = list(dict.values(mass_iso_ele_dict))
mass_iso_ele_sum = sum(np.array(mass_iso_ele_list))
avo_divided = pt.constants.avogadro_number/mass_iso_ele_sum
mixed_atoms_per_cm3 = sample_density * avo_divided
# Use function: mixed_atoms_per_cm3 = _functions.atoms_per_cm3(density=sample_density, mass=mass_iso_ele_sum)

# sum of (sigma * ele_ratio * iso_ratio * l)
yi_values_l = list(dict.values(sigma_iso_ele_sum_l_eledict))
yi_values_l_sum = sum(yi_values_l)
# sum of (sigma * ele_ratio * iso_ratio)
yi_values = list(dict.values(sigma_iso_ele_sum_eledict))
yi_values_sum = sum(yi_values)


trans_sum = _functions.sigl2trans_quick(mixed_atoms_per_cm3, yi_values_l_sum)
y_trans_tot = trans_sum


def get_residual():


    model = _functions.sig2trans_quick(_thick_cm, mixed_atoms_per_cm3, sigma_iso_ele_sum)
    resudual = data - model
    return

trans_sum = _functions.sig2trans_quick(_thick_cm, mixed_atoms_per_cm3, sigma_iso_ele_sum)
attenu_sum = 1 - trans_sum

plt.plot(x_energy, attenu_sum, 'b-', label=_input_element+'-'+_database)

plt.xlim(0, 300)
# plt.ylim(-0.01, 1.01)
# plt.xlabel(_x_words)
# plt.ylabel(_y_words)
plt.legend(loc='best')
plt.show()
# Transmission calculation of summed and separated contributions by each isotopes


# def get_distance_delay(paramas, ):





# ### Create the trans or absorb dict of ele for plotting if needed
# if _plot_each_ele_contribution == 'Y':
#     y_ele_dict = {}
#     for _ele in elements:
#         if _trans_y_axis == 'Y':
#             y_ele_dict[_ele] = _functions.sigl2trans_quick(mixed_atoms_per_cm3, sigma_iso_ele_sum_l_eledict[_ele])
#         else:
#             y_ele_dict[_ele] = 1 - _functions.sigl2trans_quick(mixed_atoms_per_cm3, sigma_iso_ele_sum_l_eledict[_ele])
#
#
# ### Create the trans or absorb dict : y_iso_dicts of isotopes for plotting if needed
# if _plot_each_iso_contribution == 'Y':
#     y_iso_dicts = {}
#     y_iso_dict = {}
#     for _ele in elements:
#         for _iso in isotopes_dict[_ele]:
#             if _trans_y_axis == 'Y':
#                 y_iso_dict[_iso] = _functions.sigl2trans_quick(mixed_atoms_per_cm3,
#                                                                sigma_iso_ele_l_eleisodict[_ele][_iso])
#             else:
#                 y_iso_dict[_iso] = 1 - _functions.sigl2trans_quick(mixed_atoms_per_cm3,
#                                                                    sigma_iso_ele_l_eleisodict[_ele][_iso])
#         y_iso_dicts[_ele] = y_iso_dict
#         y_iso_dict = {}  # Clear for following set of isotopes
#     # print(y_iso_dicts)
#
#
# # ### Export to clipboard for density and thickness manipulations with Excel or DataGraph
# # _name = _input_formula
# # df_yi_tot = pd.DataFrame(data=x_energy, index=None)
# # df_yi_tot.rename(columns={0: 'eV'+_name}, inplace=True)
# # df_yi_tot['lamda-'+_name] = _functions.ev2lamda(x_energy)
# # df_yi_tot['sample_density-'+_name] = sample_density
# # df_yi_tot['avo_divided-'+_name] = avo_divided
# # df_yi_tot['sigma-'+_name] = yi_values_sum
# #
# #
# # # print(y_i_iso_ele_sum_dict)
# # for _each_ in elements:
# #     _ele_str = str(_each_)
# #     df_yi_tot['sigma-'+_ele_str] = sigma_iso_ele_sum_eledict[_each_]
# #     df_test = pd.DataFrame(sigma_iso_ele_eleisodict[_each_])
# #     df_yi_tot = pd.concat([df_yi_tot, df_test], axis=1)
# #
# # print(df_yi_tot.head())
# # # # Export to clipboard
# # # df_yi_tot.to_clipboard(excel=True)
#
#
# ### Plot the theoretical neutron resonance
# if _plot_or_not == 'Y':
#     ### Determine x y axis types and captions
#     if _energy_x_axis == 'Y':
#         _x_axis = x_energy
#         _x_words = 'Energy (eV)'
#     else:
#         _x_axis = _functions.ev2lamda(x_energy)
#         _x_words = 'Wavelength (Å)'
#
#     if _trans_y_axis == 'Y':
#         _y_words = 'Neutron transmission'
#     else:
#         _y_words = 'Neutron attenuation'
#
#     ### Determine x y axis values
#     if _plot_mixed == 'Y':
#         if _trans_y_axis == 'Y':
#             _y_axis = y_trans_tot
#         else:
#             _y_axis = 1 - y_trans_tot
#         plt.plot(_x_axis, _y_axis, label=_input_formula)
#
#     if _plot_each_ele_contribution == 'Y':
#         for _ele in elements:
#             _y_each_axis = y_ele_dict[_ele]
#             plt.plot(_x_axis, _y_each_axis, label=_ele)
#
#     if _plot_each_iso_contribution == 'Y':
#         for _ele in elements:
#             for _iso in isotopes_dict[_ele]:
#                 _y_each_axis = y_iso_dicts[_ele][_iso]
#                 plt.plot(_x_axis, _y_each_axis, label=_iso)
#
#     plt.ylim(-0.01, 1.01)
#     plt.xlabel(_x_words)
#     plt.ylabel(_y_words)
#     plt.legend(loc='best')
#     plt.show()


