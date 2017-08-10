import matplotlib.pyplot as plt
import peakutils as pku

import _functions
import _plot_functions
# from scipy.optimize import leastsq
import detect_peaks

# Input sample name or names as str, case sensitive
_input_ele_str = 'Ag'  # input('Please input the chemicals? ')
_input_thick_mm = 0.025  # float(input('Please input the thickness or majority thickness of stacked foils in mm : '))
_database = 'ENDF_VIII'
energy_max = 200  # max incident energy in eV
energy_min = 0  # min incident energy in eV
energy_sub = 100

# Ideal
x_energy, y_trans_tot = _plot_functions.get_tot_trans_for_single_ele(_input_ele_str, _input_thick_mm, energy_max, energy_min, energy_sub)
y_attenu_tot = 1 - y_trans_tot
# print('x_ideal: ', x_energy)
# print('y_ideal: ', y_attenu_tot)
ideal_y_index = pku.indexes(y_attenu_tot, thres=0.05, min_dist=50)#, thres=0.1, min_dist=50)
ideal_x_index = pku.interpolate(x_energy, y_attenu_tot, ind=ideal_y_index)
print('x_ideal_peak: ', ideal_x_index)
# peaks_ind = pku.peak.indexes(y_attenu_tot, min_dist=50)
# print(peaks_ind)
# print('y_ideal_peak: ', ideal_y_index)
plt.plot(x_energy, y_attenu_tot, 'b-', label=_input_ele_str+'_ideal')
plt.plot(x_energy[ideal_y_index], y_attenu_tot[ideal_y_index], 'bo', label='peak_ideal')


# Experiment
source_to_detector_cm = 1607.1  # cm
delay_ms = -12.112952089537648#-12.145 #4.5 - 16.61295379  # ms
delay_us = delay_ms * 1000
_slice = 220
energy_min = 0
time_lamda_ev_axis = 'eV'
_name = 'foil7'
data_path = 'data/' + _name + '.csv'
spectra_path = 'data/spectra.txt'
x_data_array = _functions.get_spectra_slice(spectra_path, time_lamda_ev_axis, delay_us,
                                            source_to_detector_cm, _slice)
# print('x_exp: ', x_data_array)
y_data_array = 1 - _functions.get_normalized_data_slice(data_path, _slice) / 4.25
# print('y_exp: ', y_data_array)

exp_y_index = detect_peaks.detect_peaks(y_data_array, mph=0.07, mpd=10)
# exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
# exp_y_index = pku.indexes(y_data_array, thres=0.17, min_dist=50)
# exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
print('x_exp_peak: ', exp_y_index)
plt.plot(x_data_array, y_data_array, 'r-', label=_name)
plt.plot(x_data_array[exp_y_index], y_data_array[exp_y_index], 'go')
# baseline = pku.baseline(y_data_array)
# print(baseline)
# print('y_exp_peak: ', exp_y_index)


# params = Parameters()
# params.add('source_to_detector_cm', value=source_to_detector_cm)
# params.add('delay_us', value=delay_us)
#
# x_gap = _fit_funtions.peak_x_gap(params, ideal_x_index, y_data_array)
# print('x_gap:', x_gap)

# out = minimize(_fit_funtions.peak_x_gap, params, method='leastsq', args=(ideal_x_index, y_data_array))
# out = scipy.optimize.minimize(_fit_funtions.peak_x_gap_scipy, delay_us, method='leastsq', args=(ideal_x_index, y_data_array))
# print(out.__dict__)

# plt.plot(x_data_array, scipy.signal.detrend(y_data_array), 'g-', label='detrend')
# plt.plot(x_data_array[exp_y_index], y_data_array[exp_y_index], 'gx', label='peak_exp')

plt.title('Peak estimation')
plt.ylim(-0.01, 1.01)
plt.xlim(energy_min, energy_max)
plt.legend(loc='best')
plt.show()
