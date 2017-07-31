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
# from scipy.optimize import leastsq
import scipy.optimize
import pprint

# Input sample name or names as str, case sensitive
_input_ele_str = 'CoAg'  # input('Please input the chemicals? ')
_input_thick_mm = 0.05  # float(input('Please input the thickness or majority thickness of stacked foils in mm : '))
_input_thick_mm_list = [0.025, 0.025]
_database = 'ENDF_VIII'
energy_max = 300  # max incident energy in eV
energy_min = 0  # min incident energy in eV
energy_sub = 100
sub_x = (energy_max - energy_min) * energy_sub

thick_cm_list = np.array(_input_thick_mm_list) / 10
# Ideal
formula_dict = _functions.input2formula(_input_ele_str)
elements = _functions.dict_key_list(formula_dict)
ratios = _functions.dict_value_list(formula_dict)
sum_ratios = sum(ratios)
thick_cm_dict = _functions.dict_value_by_key(elements, thick_cm_list)
# stoichiometric_ratio = _functions.ele_ratio_dict(elements, thick_cm_dict, density_gcm3_dict, molar_mass_dict)
isotope_dict = _functions.get_isotope_dicts(_database, elements)
sigma_dicts = {}
for ele in elements:
    file_names = _functions.get_file_path(_database, ele)
    x_energy, sigma_dicts[ele] = _functions.get_sigma(isotope_dict[ele],
                                                      file_names,
                                                      energy_min,
                                                      energy_max,
                                                      energy_sub)
pprint.pprint(sigma_dicts)
print(x_energy)

params = Parameters()
# params.add('thick_cm_dict', value=thick_cm_dict)
# params.add('density_gcm3_dict', value=density_gcm3_dict)


# thick_cm_dict = _functions.dict_value_by_key(elements, thick_cm_list)
# print(elements)
# print(isotope_dict)
# print(sigma_iso_ele_eleisodict)
# print(sigma_iso_ele_sum_eledict)

# thick_cm_list = np.array(_input_thick_mm_list) / 10
# mixed_l_n_avo = _plot_functions.l_x_n_multi_ele_stack(elements,
#                                                       thick_cm_dict,
#                                                       density_gcm3_dict,
#                                                       molar_mass_dict)
#
# yi_values = list(dict.values(sigma_iso_ele_sum_eledict))
# yi_values_sum = sum(yi_values)
# # sum of (sigma * ele_ratio * iso_ratio)
# # print(yi_values)
# trans_sum = _functions.sig_l_2trans_quick(mixed_l_n_avo, yi_values_sum)
# y_trans_tot = trans_sum
#
# y_attenu_tot = 1 - y_trans_tot
#
# ideal_y_index = pku.indexes(y_attenu_tot, thres=0.15, min_dist=10)#, thres=0.1, min_dist=50)
# ideal_x_index = pku.interpolate(x_energy, y_attenu_tot, ind=ideal_y_index)
# print('x_ideal_peak: ', ideal_x_index)
# #
# plt.plot(x_energy, y_attenu_tot, 'b-', label=_input_ele_str+'_ideal')
# plt.plot(x_energy[ideal_y_index], y_attenu_tot[ideal_y_index], 'bo', label='peak_ideal')
#

# # Experiment
# source_to_detector_cm = 1610.1796635603498  # cm
# delay_ms = -12.112641168682274#-12.145 #4.5 - 16.61295379  # ms
# delay_us = delay_ms * 1000
# range_min = 500
# range_max = 2771
# _slice = range_min
# energy_min = 0
# time_lamda_ev_axis = 'eV'
# _name = 'foil6'
# data_path = 'data/' + _name + '.csv'
# spectra_path = 'data/spectra.txt'
# x_data_array = _functions.get_spectra_range(spectra_path, delay_us,
#                                             source_to_detector_cm, range_min, range_max)
# # print('x_exp: ', x_data_array)
# y_data_array = 1 - _functions.get_normalized_data_range(data_path, range_min, range_max)/4.25
# # print('y_exp: ', y_data_array)
# exp_y_index = pku.indexes(y_data_array, thres=0.12/max(y_data_array), min_dist=7)
# exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
# print('x_exp_peak: ', exp_x_index)
# print('Equal size: ', len(ideal_x_index) == len(exp_x_index))
# # baseline = pku.baseline(y_data_array)
# # print(baseline)
# # print('y_exp_peak: ', exp_y_index)
# #
# # df = pd.DataFrame()
# # df['Exp_x'] = x_data_array
# # df['Exp_y'] = y_data_array
# # df2 = pd.DataFrame()
# # df2['Ideal_x'] = x_energy
# # df2['Ideal_y'] = y_attenu_tot
# #
# # params = Parameters()
# # params.add('source_to_detector_cm', value=source_to_detector_cm)
# # params.add('delay_us', value=delay_us)
#
# # x_gap = _fit_funtions.peak_x_gap(params, ideal_x_index, y_data_array)
# # print('x_gap:', x_gap)
#
# # out = minimize(_fit_funtions.peak_x_gap, params, method='leastsq', args=(ideal_x_index, y_data_array))
# # out = scipy.optimize.minimize(_fit_funtions.peak_x_gap_scipy, delay_us, method='leastsq', args=(ideal_x_index, y_data_array))
# # print(out.residual)
#
# plt.plot(x_data_array, y_data_array, 'r-', label=_name)
# plt.plot(x_data_array[exp_y_index], y_data_array[exp_y_index], 'go', label='peak_exp')
# plt.title('Peak estimation')
#
# plt.ylim(-0.01, 1.01)
# plt.xlim(0, energy_max)
# plt.legend(loc='best')
# plt.show()
