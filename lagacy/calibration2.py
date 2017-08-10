import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters

import _functions
import _plot_functions

# from scipy.optimize import leastsq
# Input sample name or names as str, case sensitive
_input_ele_str = 'Hf'  # input('Please input the chemicals? ')
_input_thick_mm = 0.025  # float(input('Please input the thickness or majority thickness of stacked foils in mm : '))
_database = 'ENDF_VIII'
energy_max = 170  # max incident energy in eV
energy_min = 6  # min incident energy in eV
energy_sub = 100

# Ideal
x_energy, y_trans_tot = _plot_functions.get_tot_trans_for_single_ele(_input_ele_str, _input_thick_mm, energy_max, energy_min, energy_sub)
y_attenu_tot = 1 - y_trans_tot
# print('x_ideal: ', x_energy)
# print('y_ideal: ', y_attenu_tot)
ideal_y_index = pku.indexes(y_attenu_tot, thres=0.15, min_dist=10)#, thres=0.1, min_dist=50)
ideal_x_index = pku.interpolate(x_energy, y_attenu_tot, ind=ideal_y_index)
print('x_ideal_peak: ', ideal_x_index)
# peaks_ind = pku.peak.indexes(y_attenu_tot, min_dist=50)
# print(peaks_ind)
# print('y_ideal_peak: ', ideal_y_index)
plt.plot(x_energy, y_attenu_tot, 'b-', label=_input_ele_str+'_ideal')
plt.plot(x_energy[ideal_y_index], y_attenu_tot[ideal_y_index], 'bo', label='peak_ideal')


# Experiment
source_to_detector_cm = 1610.225604799349  # cm
delay_ms = -12.112632941196785#-12.145 #4.5 - 16.61295379  # ms
delay_us = delay_ms * 1000
range_min = 600
range_max = 2772
_slice = range_min
energy_min = 0
time_lamda_ev_axis = 'eV'
_name = 'foil8'
data_path = 'data/' + _name + '.csv'
spectra_path = 'data/spectra.txt'
x_data_array = _functions.get_spectra_range(spectra_path, time_lamda_ev_axis, delay_us,
                                            source_to_detector_cm, range_min, range_max)
x_data_array = x_data_array[::-1]
# print('x_exp: ', x_data_array)
y_data_array = 1 - _functions.get_normalized_data_range(data_path, range_min, range_max) / 4.25
y_data_array = y_data_array[::-1]
# print('y_exp: ', y_data_array)
exp_y_index = pku.indexes(y_data_array, thres=0.12/max(y_data_array), min_dist=7)
exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
print('x_exp_peak: ', exp_x_index)
print('Equal size: ', len(ideal_x_index) == len(exp_x_index))
# baseline = pku.baseline(y_data_array)
# print(baseline)
# print('y_exp_peak: ', exp_y_index)
#
# df = pd.DataFrame()
# df['Exp_x'] = x_data_array
# df['Exp_y'] = y_data_array
# df2 = pd.DataFrame()
# df2['Ideal_x'] = x_energy
# df2['Ideal_y'] = y_attenu_tot

params = Parameters()
params.add('source_to_detector_cm', value=source_to_detector_cm)
params.add('delay_us', value=delay_us)

# x_gap = _fit_funtions.peak_x_gap(params, ideal_x_index, y_data_array)
# print('x_gap:', x_gap)

# out = minimize(_fit_funtions.peak_x_gap, params, method='leastsq', args=(ideal_x_index, y_data_array))
# out = scipy.optimize.minimize(_fit_funtions.peak_x_gap_scipy, delay_us, method='leastsq', args=(ideal_x_index, y_data_array))
# print(out.__dict__)

plt.plot(x_data_array, y_data_array, 'r-', label=_name)
plt.plot(x_data_array[exp_y_index], y_data_array[exp_y_index], 'go', label='peak_exp')
plt.title('Peak estimation')

plt.ylim(-0.01, 1.01)
plt.xlim(0, energy_max)
plt.legend(loc='best')
plt.show()
