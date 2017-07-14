import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import _functions
import glob

# Parameters
source_to_detector_cm = 1610.9  # cm
delay_ms = 4.5 - 16.611  # ms
delay_us = delay_ms * 1000
_slice = 220
time_lamda_ev_axis = 'eV'
_name = 'foil'
path = 'data/*' + _name + '*.csv'
file_names = glob.glob(path)
x_axis_array = _functions.get_spectra_slice('data/data_spectra_20.txt', time_lamda_ev_axis, delay_us, source_to_detector_cm, _slice)
print(x_axis_array)
df = pd.DataFrame()
key_dict = {}
for i in range(len(file_names)):
    key_dict[i] = _name + str(i+1)
# df[time_lamda_ev_axis] = x_axis_array
key_list = list(dict.values(key_dict))
plot_boo_dict = _functions.boo_dict(key_list, 'N')
plot_foil_num = [1, 3, 4, 5, 6, 7, 8]
plot_foil = {}
for num in plot_foil_num:
    plot_foil[num] = _name+str(num)
plot_foil_list = list(dict.values(plot_foil))
print(plot_foil_list)
plot_boo_dict = _functions.boo_dict_invert_by_key(plot_foil_list, plot_boo_dict)
print(plot_boo_dict)
for foil in key_list:
    df[foil] = -1*_functions.get_normalized_data_slice('data/'+foil+'.csv', _slice)
    if plot_boo_dict[foil] == 'Y':
        plt.plot(x_axis_array, df[foil], '.', label=foil)
print(df.head())
print(df.tail())

plt.xlim(0, 300)
# plt.ylim(-0.01, 1.01)ß
plt.legend(loc=5)
plt.xlabel(time_lamda_ev_axis)
plt.show()


