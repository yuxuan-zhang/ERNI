import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
import _utilities
from resonance import Resonance

# Global parameters
_energy_min = 8
_energy_max = 170
_energy_step = 0.01
# Input sample name or names as str, case sensitive
_layer_1 = 'Ag'
_thickness_1 = 0.025 # mm
# _density_1 = 8 # g/cm3 deviated due to porosity

o_reso = Resonance(energy_min=_energy_min, energy_max=_energy_max, energy_step=_energy_step)
o_reso.add_layer(formula=_layer_1, thickness=_thickness_1)

# Ideal
x_energy = o_reso.stack_sigma['Ag']['Ag']['energy_eV']
sigma_b_ = o_reso.stack_sigma['Ag']['Ag']['sigma_b']
y_attenu_tot = o_reso.stack_signal['Ag']['attenuation']
# print('x_ideal: ', x_energy)
# print('y_ideal: ', y_attenu_tot)
ideal_y_index = pku.indexes(y_attenu_tot, thres=0.15, min_dist=10)#, thres=0.1, min_dist=50)
ideal_x_index = pku.interpolate(x_energy, y_attenu_tot, ind=ideal_y_index)
print('x_ideal_peak: ', ideal_x_index)
# peaks_ind = pku.peak.indexes(y_attenu_tot, min_dist=50)
# print(peaks_ind)
# print('y_ideal_peak: ', ideal_y_index)
plt.plot(x_energy, y_attenu_tot, 'b-', label=_layer_1+'_ideal')
plt.plot(x_energy[ideal_y_index], y_attenu_tot[ideal_y_index], 'bo', label='peak_ideal')


# Experiment
source_to_detector_cm = 1612.3278721983177  # cm
delay_ms = -12.112494119089204#-12.145 #4.5 - 16.61295379  # ms
delay_us = delay_ms * 1000
range_min = 500
range_max = 2000
_slice = range_min
energy_min = 0
time_lamda_ev_axis = 'eV'
_name = 'foil7'
data_path = 'data/' + _name + '.csv'
spectra_path = 'data/spectra.txt'
x_data_array = _utilities.get_spectra_range(spectra_path, delay_us,
                                            source_to_detector_cm, range_min, range_max)
# print('x_exp: ', x_data_array)
y_data_array = 1 - _utilities.get_normalized_data_range(data_path, range_min, range_max) / 4.25
# print('y_exp: ', y_data_array)
exp_y_index = pku.indexes(y_data_array, thres=0.12/max(y_data_array), min_dist=7)
exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
print('x_exp_peak: ', exp_x_index)
equal_size_boo = len(ideal_x_index) == len(exp_x_index)
print('Equal size: ', equal_size_boo)


# Fitting the peak positions

params = Parameters()
params.add('source_to_detector_cm', value=source_to_detector_cm)
params.add('delay_us', value=delay_us)

plt.plot(x_data_array, y_data_array, 'r-', label=_name)
plt.plot(x_data_array[exp_y_index], y_data_array[exp_y_index], 'go', label='peak_exp')
plt.title('Peak estimation')

plt.ylim(-0.01, 1.01)
plt.xlim(0, _energy_max)
plt.legend(loc='best')
plt.show()

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
# x_gap = _fit_functions.peak_x_gap(params, ideal_x_index, y_data_array)
# print('x_gap:', x_gap)

# out = minimize(_fit_functions.peak_x_gap, params, method='leastsq', args=(ideal_x_index, y_data_array))
# out = scipy.optimize.minimize(_fit_funtions.peak_x_gap_scipy, delay_us, method='leastsq', args=(ideal_x_index, y_data_array))
# print(out.__dict__)