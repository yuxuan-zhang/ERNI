import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import _functions

# Parameters
source_to_detector_cm = 1612.3278721983177  # cm
delay_ms = -12.112494119089204#-12.145 #4.5 - 16.61295379  # ms
data_20 = pd.read_csv('data/data_spectra_20.txt', sep='\t', header=None)
data_40 = pd.read_csv('data/data_spectra_40.txt', sep='\t', header=None)
data_3749 = pd.read_csv('data/data_spectra_3749.txt', sep='\t', header=None)
ob_18 = pd.read_csv('data/OB_image18_pectra.txt', sep='\t', header=None)
ob_delme = pd.read_csv('data/OB_delme_spectra.txt', sep='\t', header=None)
ob_Cd = pd.read_csv('data/OB_Cd_spectra.txt', sep='\t', header=None)
df_all = pd.DataFrame()
df_all['eV'] = _functions.time2ev(data_20[0], delay_ms, source_to_detector_cm)
df_all['lamda'] = _functions.time2lamda(data_20[0], delay_ms, source_to_detector_cm)
df_all['20'] = data_20[1]
df_all['40'] = data_40[1]
df_all['3749'] = data_3749[1]
df_all['ob_18'] = ob_18[1]
df_all['ob_delme'] = ob_delme[1]
df_all['ob_Cd'] = ob_Cd[1]
df_all.set_index(df_all['eV'], inplace=True)

plt.plot(df_all)
# plt.xlim(-0.01, 1.01)
plt.show()

# df_a = pd.read_csv('Analysis_07_11_2017.csv')

print(df_all.head())


