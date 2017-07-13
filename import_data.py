import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import _functions
import glob

# Parameters
source_to_detector_cm = 1610.9  # cm
delay_ms = 4.5 - 16.6  # ms
delay_us = delay_ms * 1000
_slice = 220
time_lamda_ev_axis = 'eV'
_name = 'foil'
path = 'data/' + _name + '*.csv'
file_names = glob.glob(path)
x_axis_array = _functions.get_spectra_slice('data/data_spectra_20.txt', time_lamda_ev_axis, delay_us, source_to_detector_cm, _slice)
df = pd.DataFrame()
key_dict = {}
for i in range(len(file_names)):
    key_dict[i] = _name + str(i+1)
# df[time_lamda_ev_axis] = x_axis_array
key_list = list(dict.values(key_dict))
plot_boo_dict = _functions.boo_dict(key_list, 'Y')

plot_foil = [_name+'3']
plot_boo_dict = _functions.boo_dict_invert_by_key(plot_foil, plot_boo_dict)
print(plot_boo_dict)
for foil in key_list:
    df[foil] = _functions.get_normalized_data_slice('data/'+foil+'.csv', _slice)
    if plot_boo_dict[foil] == 'Y':
        plt.plot(x_axis_array, df[foil], '--')
print(df.head())
print(df.tail())

# df_all.plot()
# df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil1')
# df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil2')
# df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil3')
# df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil4')
# df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil5')
# df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil6')
# df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil7')
# df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil8')

# plt.xlim(0.00453, 0.00492)
plt.xlim(0, 300)
# plt.ylim(-0.01, 1.01)
# plt.xlabel(_x_words)
# plt.ylabel(_y_words)
#
# plt.legend(loc=2)
# plt.xlabel(time_lamda_ev_axis)
plt.show()


