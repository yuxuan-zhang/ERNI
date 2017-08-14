import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lagacy import _functions

# Parameters
source_to_detector_cm = 1612.3278721983177  # cm
delay_ms = -12.112494119089204
# source_to_detector_cm = 1612.3278721983177  # cm
# delay_ms = 0.00299  # ms
delay_us = delay_ms * 1000

data_20 = pd.read_csv('data/data_spectra_20.txt', sep='\t', header=None)
data_40 = pd.read_csv('data/data_spectra_40.txt', sep='\t', header=None)
data_3749 = pd.read_csv('data/data_spectra_3749.txt', sep='\t', header=None)
ob_18 = pd.read_csv('data/OB_image18_pectra.txt', sep='\t', header=None)
ob_delme = pd.read_csv('data/OB_delme_spectra.txt', sep='\t', header=None)
ob_Cd = pd.read_csv('data/OB_Cd_spectra.txt', sep='\t', header=None)
df_all = pd.DataFrame()
df_all['eV'] = _functions.time2ev(data_20[0], delay_ms, source_to_detector_cm)
df_all['lamda'] = _functions.time2lamda(data_20[0], delay_ms, source_to_detector_cm)
df_all['time'] = data_20[0]
df_all['20'] = data_20[1]
df_all['40'] = data_40[1]
df_all['3749'] = data_3749[1]
df_all['ob_18'] = ob_18[1]*4
df_all['ob_delme'] = ob_delme[1]*4
df_all['ob_Cd'] = ob_Cd[1]*4
df_all.set_index(df_all['eV'], inplace=True)

_name = 'foil3'
data_path = 'data/' + _name + '.csv'
spectra_path = 'data/spectra.txt'
x_data_array = _functions.get_spectra(spectra_path, delay_us, source_to_detector_cm, time_lamda_ev_axis='time')
df = pd.read_csv(data_path, header=None, skiprows=1)
data_array = np.array(df[1])
data = data_array[:int(len(data_array)/2)] * 100000
ob = data_array[int(len(data_array)/2):] * 300000

# plt.plot(x_data_array, ob)
# plt.xlim(-0.01, 500)
df_all['foil6_ob'] = ob
df_all['foil6_data'] = data


df_all.plot.line()
plt.xlim(0, 200)
plt.legend(loc='best')
plt.show()


# df_april = pd.DataFrame()
# df_22 = pd.read_csv('data/Image022_Spectra.txt', sep='\t', header=None)
# df_23 = pd.read_csv('data/Image023_Spectra.txt', sep='\t', header=None)
# df_april['eV'] = _functions.time2ev(df_22[0], delay_ms, source_to_detector_cm)
# df_april['lamda'] = _functions.time2lamda(df_22[0], delay_ms, source_to_detector_cm)
# df_april['time'] = df_22[0]
# df_april['OB_22'] = df_22[1]
# df_april['OB_23'] = df_23[1]
# df_april.set_index(df_april['lamda'], inplace=True)
# print(df_april.head())
# df_april.plot.line()
# # plt.xlim(-0.01, 300)
# # plt.ylim(-0.01, 150000)
# plt.legend(loc='best')
# plt.show()

