import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import peakutils as pku
# from scipy.optimize import leastsq
import scipy.signal
from lmfit import Parameters

import _functions
import _plot_functions

plt.style.use('bmh')

# Input sample name or names as str, case sensitive
_input_ele_str = 'Ta'  # input('Please input the chemicals? ')
_input_thick_mm = 0.03  # float(input('Please input the thickness or majority thickness of stacked foils in mm : '))
_database = 'ENDF_VIII'
energy_max = 220  # max incident energy in eV
energy_min = 5  # min incident energy in eV
energy_sub = 100

# Ideal
x_energy, y_trans_tot = _plot_functions.get_tot_trans_for_single_ele(_input_ele_str, _input_thick_mm, energy_max, energy_min, energy_sub)
y_attenu_tot = 1 - y_trans_tot
# print('x_ideal: ', x_energy)
# print('y_ideal: ', y_attenu_tot)
ideal_y_index = pku.indexes(y_attenu_tot, thres=0.05, min_dist=10)#, thres=0.1, min_dist=50)
ideal_x_index = pku.interpolate(x_energy, y_attenu_tot, ind=ideal_y_index)
print('x_ideal_peak: ', ideal_x_index)
# peaks_ind = pku.peak.indexes(y_attenu_tot, min_dist=50)
# print(peaks_ind)
# print('y_ideal_peak: ', ideal_y_index)
plt.plot(x_energy, y_attenu_tot, 'b-', label='Ideal', linewidth=1)
plt.plot(x_energy[ideal_y_index], y_attenu_tot[ideal_y_index], 'bx')#, label='peak_ideal')


# Experiment
source_to_detector_cm = 1612.5  # cm
delay_ms = 0.00299#-12.145 #4.5 - 16.61295379  # ms
delay_us = delay_ms * 1000
range_min = 0
range_max = 1992
_slice = range_min
energy_min = 0
time_lamda_ev_axis = 'eV'
_name = 'Exp. data'
data_path = 'data/July_2017/Values1.csv'
ob_path = 'data/July_2017/Values_ob.csv'
# data_path = 'data/' + _name + '.csv'
spectra_path = 'data/July_2017/Image052_Spectra.txt'
x_data_array = _functions.get_spectra_range(spectra_path, delay_us,
                                            source_to_detector_cm, range_min, range_max)
# print('x_exp: ', x_data_array)
df = pd.read_csv(data_path, header=None, skiprows=1)
df2 = pd.read_csv(ob_path, header=None, skiprows=1)

# y_data_array = np.array(df[1]) / np.array(df2[1])
y_data_array = np.array(df[1])
y_data_array = y_data_array[::-1]  # Flip array from descending to normal
print(y_data_array)
y_data_array = scipy.signal.detrend(y_data_array)

y_data_array = -y_data_array/1.25
y_data_array = y_data_array[range_min:range_max]
# y_data_array = 1 - _functions.get_normalized_data_range(data_path, range_min, range_max)/4.25
# print('y_exp: ', y_data_array)
exp_y_index = pku.indexes(y_data_array, thres=0.08/max(y_data_array), min_dist=10)
exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
# print('x_exp_peak: ', exp_x_index)
# print('Equal size: ', len(ideal_x_index) == len(exp_x_index))
# baseline = pku.baseline(y_data_array)
# print(baseline)
# print('y_exp_peak: ', exp_y_index)


df1 = pd.DataFrame()
df2 = pd.DataFrame()
time_array = _functions.get_spectra_range(spectra_path, delay_us,
                                          source_to_detector_cm, range_min, range_max, time_lamda_ev_axis='time')
df1['Exp_x'] = time_array
df1['Exp_y'] = y_data_array[::-1]
df2['Ideal_x'] = x_energy
df2['Ideal_y'] = y_attenu_tot

df1.to_clipboard(excel=True)

params = Parameters()
params.add('source_to_detector_cm', value=source_to_detector_cm)
params.add('delay_us', value=delay_us)

# x_gap = _fit_funtions.peak_x_gap(params, ideal_x_index, y_data_array)
# print('x_gap:', x_gap)

# out = minimize(_fit_funtions.peak_x_gap, params, method='leastsq', args=(ideal_x_index, y_data_array))
# out = scipy.optimize.minimize(_fit_funtions.peak_x_gap_scipy, delay_us, method='leastsq', args=(ideal_x_index, y_data_array))
# print(out.__dict__)

plt.plot(x_data_array, y_data_array, 'r-', label='Exp.', markersize=4, alpha=0.7)
plt.plot(x_data_array[exp_y_index], y_data_array[exp_y_index], 'go', label='peak_exp', markersize=3)

plt.title('Neutron resonances of Ta foil')
plt.ylim(-0.01, 0.6)
plt.xlim(0, energy_max)
plt.xlabel('Energy (eV)')
plt.ylabel('Neutron attenuation')
plt.legend(loc='best')
plt.show()
