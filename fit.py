import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import _functions
import glob

# Parameters
source_to_detector_cm = 1610.9  # cm
delay_ms = 4.5 - 15  # ms
delay_us = delay_ms * 1000
time_lamda_ev_axis = 'lamda'
_name = 'foil'
path = 'data/' + _name + '*.csv'
file_names = glob.glob(path)
df_spectra = pd.read_csv('data/data_spectra_40.txt', sep='\t', header=None)
# df_spectra = pd.read_csv('data/*11_Spectra.txt', sep='\t', header=None)
print(file_names)
print(len(df_spectra))
df_all = pd.DataFrame()
time_array = (np.array(df_spectra[0]))
if time_lamda_ev_axis == 't':
    df_all['time'] = time_array
    df_all.set_index(df_all['time'], inplace=True)
    # del df_all['time']
if time_lamda_ev_axis == 'lamda':
    lamda_array = _functions.time2lamda(time_array, delay_us, source_to_detector_cm)
    df_all['lamda'] = lamda_array
    df_all.set_index(df_all['lamda'], inplace=True)
    # del df_all['lamda']
if time_lamda_ev_axis == 'eV':
    ev_array = _functions.time2ev(time_array, delay_us, source_to_detector_cm)
    df_all['eV'] = ev_array
    df_all.set_index(df_all['eV'], inplace=True)
    # del df_all['eV']
p = 1
for _files in file_names:
    df = pd.read_csv(_files, header=None, skiprows=1)
    # df1 = df[::len(df)/2]  # drop rows beyond range
    # df2 = df[len(df)/2::-1]
    data_array = np.array(df[1])
    data = data_array[:int(len(data_array)/2)]
    ob = data_array[int(len(data_array)/2):]
    normalized = data/ob
    # OB at the end of 2773
    df_all[_name+str(p)] = -1 * normalized
    p = p+1
print(df_all.head())
print(df_all.tail())

# df_all.plot()
df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil1')
df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil2')
df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil3')
df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil4')
df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil5')
df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil6')
df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil7')
df_all.plot(kind='scatter', x=time_lamda_ev_axis, y='foil8')

# plt.xlim(0.00453, 0.00492)
# plt.xlim(0, 300)
# plt.ylim(-0.01, 1.01)
# plt.xlabel(_x_words)
# plt.ylabel(_y_words)

plt.legend(loc=2)
plt.xlabel(time_lamda_ev_axis)
plt.show()


