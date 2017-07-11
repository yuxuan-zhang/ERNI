import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import _functions

# Parameters
source_to_detector_cm = 1610.9  # cm
delay_ms = 4.5  # ms
data_20 = pd.read_csv('data_spectra_20.txt', sep='\t', header=None)
data_40 = pd.read_csv('data_spectra_40.txt', sep='\t', header=None)
data_3749 = pd.read_csv('data_spectra_3749.txt', sep='\t', header=None)
ob_18 = pd.read_csv('OB_image18_pectra.txt', sep='\t', header=None)
ob_delme = pd.read_csv('OB_delme_spectra.txt', sep='\t', header=None)
ob_Cd = pd.read_csv('OB_Cd_spectra.txt', sep='\t', header=None)
df_all = pd.DataFrame()
df_all['eV'] = _functions.time2ev(data_20[0], delay_ms, source_to_detector_cm)
df_all['20'] = data_20[1]
df_all['40'] = data_40[1]
df_all['3749'] = data_3749[1]
df_all['ob_18'] = ob_18[1]
df_all['ob_delme'] = ob_delme[1]
df_all['ob_Cd'] = ob_Cd[1]
df_all.set_index(df_all['eV'], inplace=True)
# print(len(data_20))
# print(len(data_40))
# print(len(data_3749))
#
# print(len(ob18))
# print(len(ob_delme))
# print(len(ob_Cd))
# print(data_20[0] == data_3749[0])
print(df_all)

# plt.plot(data_20[0], data_20[1], label='20')
# plt.plot(data_40[0], data_40[1], label='40')
# plt.plot(data_3749[0], data_3749[1], label='40')
# plt.plot(ob18[1], label='OB18')
# plt.plot(ob_Cd[1], label='OB_Cd')
# plt.plot(ob_delme[1], label='OB_delme')

# plt.plot(df_all['eV'], data_20[1], label='20')
# plt.plot(df_all['eV'], data_40[1], label='40')
# plt.plot(df_all['eV'], data_3749[1], label='40')
# plt.plot(df_all['eV'], ob_18[1], label='OB_18')
# plt.plot(df_all['eV'], ob_Cd[1], label='OB_Cd')
# plt.plot(df_all['eV'], ob_delme[1], label='OB_delme')
#
# plt.legend(loc='best')
# plt.show()

df_all.plot()
plt.show()
