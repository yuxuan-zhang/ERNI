import numpy as np
import periodictable as pt
from periodictable import constants
import re
import os
import glob
import pandas as pd
import _functions
from scipy.interpolate import *
import peakutils as pku


def peak_x_gap(params, ideal_x_index, y_data_array):
    # Unpack Parameters:
    parvals = params.valuesdict()
    source_to_detector_cm = parvals['source_to_detector_cm']
    delay_us = parvals['delay_us']
    # Model:
    time_lamda_ev_axis = 'eV'
    spectra_path = 'data/spectra.txt'
    _slice = 220
    x_data_array = _functions.get_spectra_slice(spectra_path, time_lamda_ev_axis, delay_us,
                                                source_to_detector_cm, _slice)
    exp_y_index = pku.indexes(y_data_array, thres=0.5, min_dist=100)
    exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
    gap = (exp_x_index[0] - ideal_x_index) ** 2
    return gap


def peak_x_gap_scipy(delay_us, ideal_x_index, y_data_array):
    # Unpack Parameters:
    source_to_detector_cm = 1610.9  # cm
    # Model:
    time_lamda_ev_axis = 'eV'
    spectra_path = 'data/spectra.txt'
    _slice = 220
    x_data_array = _functions.get_spectra_slice(spectra_path, time_lamda_ev_axis, delay_us,
                                                source_to_detector_cm, _slice)
    exp_y_index = pku.indexes(y_data_array, thres=0.6, min_dist=50)
    exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
    gap = (exp_x_index[0] - ideal_x_index) ** 2
    return gap


def peak_gap(params, ideal_x_index):
    # Unpack Parameters:
    parvals = params.valuesdict()
    source_to_detector_cm = parvals['source_to_detector_cm']
    time_us = parvals['time_us']
    # Model:
    energy_miliev = 81.787 / (0.3956 * time_us / source_to_detector_cm) ** 2
    energy_ev = energy_miliev / 1000
    return (energy_ev - ideal_x_index) ** 2

# def get_multi_data(file_name_signature, time_lamda_ev_axis, delay_us, source_to_detector_cm, _slice):
#     # time_lamda_ev_axis = 'eV'
#     # _file_name_signature = 'foil*'
#     path = 'data/*' + file_name_signature + '*.csv'
#     spectra_path = 'data/spectra.txt'
#     file_names = glob.glob(path)
#     x_axis_array = _functions.get_spectra_slice(spectra_path, time_lamda_ev_axis, delay_us,
#                                                 source_to_detector_cm, _slice)
#     df = pd.DataFrame()
#     key_dict = {}
#     for i in range(len(file_names)):
#         key_dict[i] = file_name_signature + str(i + 1)
#     # df[time_lamda_ev_axis] = x_axis_array
#     key_list = list(dict.values(key_dict))
#     plot_boo_dict = _functions.boo_dict(key_list, 'N')
#     plot_foil_num = [1, 2, 4]
#     plot_foil = {}
#     for num in plot_foil_num:
#         plot_foil[num] = file_name_signature + str(num)
#     plot_foil_list = list(dict.values(plot_foil))
#     print(plot_foil_list)
#     plot_boo_dict = _functions.boo_dict_invert_by_key(plot_foil_list, plot_boo_dict)
#     print(plot_boo_dict)
#     for foil in key_list:
#         df[foil] = _functions.get_normalized_data_slice('data/' + foil + '.csv', _slice)
#         if plot_boo_dict[foil] == 'Y':
#             plt.plot(x_axis_array, df[foil], '.', label=foil)
#     print(df.head())
#     print(df.tail())
#
#     # plt.xlim(0, 300)
#     # plt.ylim(-0.01, 1.01)ß
#     plt.legend(loc='best')
#     plt.xlabel(time_lamda_ev_axis)
#     plt.show()


def get_pre_data_to_fit(_database, _element):
    # main_dir = os.path.dirname(os.path.abspath(__file__))
    # path = main_dir + '/data_web/' + _database + '/' + _element + '*.csv'
    path = 'data_web/' + _database + '/' + _element + '*.csv'
    file_names = glob.glob(path)
    abundance_dict = {}
    mass_dict = {}
    for _i, file in enumerate(file_names):
        # Obtain element, z number from the basename
        _basename = os.path.basename(file)
        _name_number_csv = _basename.split('.')
        _name_number = _name_number_csv[0]
        _name = _name_number.split('-')
        _symbol = _name[1] + '-' + _name[0]
        abundance_dict[str(_symbol)] = pt.elements.isotope(_symbol).abundance / 100
        mass_dict[str(_symbol)] = pt.elements.isotope(_symbol).mass
    isotopes = list(dict.keys(abundance_dict))  # List of isotopes such as '238-U', ''235-U
    iso_abundance = list(dict.values(abundance_dict))  # List of isotopic abundance
    # iso_density = list(dict.values(density_dict))  # List of isotopic density
    iso_mass = list(dict.values(mass_dict))  # List of isotopic molar mass
    return isotopes, abundance_dict, iso_abundance, iso_mass, file_names


def get_mass_iso_ele_to_fit(iso_abundance, iso_mass, ele_at_ratio):
    # Calculate the number of atoms per unit volume
    abundance_array = np.array(iso_abundance)
    mass_array = np.array(iso_mass)
    mass_abundance_multiplied = mass_array * abundance_array
    mass_iso_ele = sum(mass_abundance_multiplied) * ele_at_ratio
    return mass_iso_ele


def get_xy_to_fit(isotopes, file_names, energy_min, energy_max, iso_abundance, sub_x, ele_at_ratio):
    # Transmission calculation of summed and separated contributions by each isotopes
    df = pd.DataFrame()
    df_raw = pd.DataFrame()
    sigma_iso_ele_isodict = {}
    # sigma_iso_ele_l_isodict = {}
    # sigma_iso_ele_l_isodict = {}
    # thick_cm = thick_mm/10
    sigma_iso_ele_sum = 0.
    # sigma_iso_ele_l_sum = 0.
    # if _natural_mix == 'Y':
    #     iso_at_ratio = iso_abundance
    # else:
    #     ratio_list = list(dict.values(_unnatural_ratio_dict))
    #     ratio_array = np.array(ratio_list)
    iso_at_ratio = iso_abundance
    for i, _isotope in enumerate(isotopes):
        # Read database .csv file
        df = pd.read_csv(file_names[i], header=1)
        # Drop rows beyond range
        df = df.drop(df[df.E_eV < energy_min].index)  # drop rows beyond range
        df = df.drop(df[df.E_eV > energy_max].index)  # drop rows beyond range
        df = df.reset_index(drop=True)  # Reset index after dropping values
        # print(df.head())
        # print(df.tail())
        '''
        Attention!!!
        The drop here not works perfect since all the data not at the same intervals.
        df1 ends at 4999 and df2 might end at 5000. 
        This will affect the accuracy of the summation performed later. 
        '''
        # Spline x-axis and y-axis for transmission calculation
        x_energy = np.linspace(df['E_eV'].min(), df['E_eV'].max(), sub_x)
        y_spline = interpolate.interp1d(x=df['E_eV'], y=df['Sig_b'], kind='linear')
        y_i = y_spline(x_energy)
        sigma_b = y_i
        # y_i_sum = y_i_sum + y_i * iso_abundance[i] * ele_at_ratio
        sigma_iso_ele_isodict[_isotope] = sigma_b * iso_at_ratio[i] * ele_at_ratio
        # sigma_iso_ele_l_isodict[_isotope] = sigma_iso_ele_isodict[_isotope] * thick_cm
        sigma_iso_ele_sum = sigma_iso_ele_sum + sigma_b * iso_at_ratio[i] * ele_at_ratio
        # sigma_iso_ele_l_sum = sigma_iso_ele_l_sum + sigma_b * iso_at_ratio[i] * ele_at_ratio * thick_cm

        """
        Attention:
        The following part is for producing df_raw of all isotopes for future reference
        """
        # Create a new DataFrame including all isotopic data
        # within the selected energy range
        first_col = _isotope + ', E_eV'
        second_col = _isotope + ', Sig_b'
        df.rename(columns={'E_eV': first_col, 'Sig_b': second_col}, inplace=True)
        df_raw = pd.concat([df_raw, df], axis=1)

    return x_energy, sigma_iso_ele_isodict, sigma_iso_ele_sum, df_raw

