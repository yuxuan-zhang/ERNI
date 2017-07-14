import numpy as np
import periodictable as pt
from periodictable import constants
import re
import os
import glob
import pandas as pd
import _functions
import matplotlib.pyplot as plt


def get_multi_data(file_name_signature, time_lamda_ev_axis, delay_us, source_to_detector_cm, _slice):
    # time_lamda_ev_axis = 'eV'
    # _file_name_signature = 'foil*'
    path = 'data/*' + file_name_signature + '*.csv'
    spectra_path = 'data/spectra.txt'
    file_names = glob.glob(path)
    x_axis_array = _functions.get_spectra_slice(spectra_path, time_lamda_ev_axis, delay_us,
                                                source_to_detector_cm, _slice)
    df = pd.DataFrame()
    key_dict = {}
    for i in range(len(file_names)):
        key_dict[i] = file_name_signature + str(i + 1)
    # df[time_lamda_ev_axis] = x_axis_array
    key_list = list(dict.values(key_dict))
    plot_boo_dict = _functions.boo_dict(key_list, 'N')
    plot_foil_num = [1, 2, 4]
    plot_foil = {}
    for num in plot_foil_num:
        plot_foil[num] = file_name_signature + str(num)
    plot_foil_list = list(dict.values(plot_foil))
    print(plot_foil_list)
    plot_boo_dict = _functions.boo_dict_invert_by_key(plot_foil_list, plot_boo_dict)
    print(plot_boo_dict)
    for foil in key_list:
        df[foil] = _functions.get_normalized_data_slice('data/' + foil + '.csv', _slice)
        if plot_boo_dict[foil] == 'Y':
            plt.plot(x_axis_array, df[foil], '.', label=foil)
    print(df.head())
    print(df.tail())

    # plt.xlim(0, 300)
    # plt.ylim(-0.01, 1.01)ß
    plt.legend(loc='best')
    plt.xlabel(time_lamda_ev_axis)
    plt.show()


# def get_ele_trans():
