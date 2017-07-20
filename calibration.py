import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import _functions
import _fit_funtions
import _plot_functions
from lmfit import minimize, Parameters
import os
import periodictable as pt
from periodictable.constants import avogadro_number
import peakutils as pku
from scipy.optimize import leastsq
from scipy.optimize import minimize


# Input sample name or names as str, case sensitive
_input_ele_str = 'AgCo'  # input('Please input the chemicals? ')
_input_thick_mm = 0.025  # float(input('Please input the thickness or majority thickness of stacked foils in mm : '))
_database = 'ENDF_VIII'
energy_max = 300  # max incident energy in eV
energy_min = 0  # min incident energy in eV
energy_sub = 100

# Ideal
x_energy, y_trans_tot = _plot_functions.get_tot_trans_for_single_ele(_input_ele_str, _input_thick_mm, energy_max, energy_min, energy_sub)
y_attenu_tot = 1 - y_trans_tot
print('x_ideal: ', x_energy)
print('y_ideal: ', y_attenu_tot)
plt.plot(x_energy, y_attenu_tot, 'b-', label=_input_ele_str+'_ideal')

ideal_y_index = pku.indexes(y_attenu_tot, thres=0.6, min_dist=50)
print(ideal_y_index)
ideal_x_index = pku.interpolate(x_energy, y_attenu_tot, ind=ideal_y_index)
print(ideal_x_index)
plt.plot(x_energy[ideal_y_index], y_attenu_tot[ideal_y_index], 'go', label='peak_ideal')


# Experiment
source_to_detector_cm = 1610.9  # cm
delay_ms = 4.5 - 16.612953  # ms
delay_us = delay_ms * 1000
_slice = 220
energy_min = 0
time_lamda_ev_axis = 'eV'
_name = 'foil6'
data_path = 'data/' + _name + '.csv'
spectra_path = 'data/spectra.txt'
x_data_array = _functions.get_spectra_slice(spectra_path, time_lamda_ev_axis, delay_us,
                                            source_to_detector_cm, _slice)
print('x_exp: ', x_data_array)
y_data_array = 1 - _functions.get_normalized_data_slice(data_path, _slice)/4.2
print('y_exp: ', y_data_array)
plt.plot(x_data_array, y_data_array, 'r-', label=_name)
exp_y_index = pku.indexes(y_data_array, thres=0.6, min_dist=50)

print(exp_y_index)
exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
print(exp_x_index)
# plt.figure(figsize=(10, 6))
varis = [delay_us, source_to_detector_cm]
gap = _fit_funtions.peak_position_gap(varis, ideal_x_index, y_data_array)
print('gap:', gap)

# minimize(_fit_funtions.peak_position_gap, delay_us, method='SLSQP')


plt.plot(x_data_array[exp_y_index], y_data_array[exp_y_index], 'x', label='peak_exp')
plt.title('First estimation')

plt.ylim(-0.01, 1.01)
plt.xlim(energy_min, energy_max)
plt.legend(loc='best')
plt.show()
