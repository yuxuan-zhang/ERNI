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
import _plot_functions
from lmfit import minimize, Parameters


def def_params_from_dict(_dict, _dict_name_str):
    dict_key = _functions.dict_key_list(_dict)
    dict_value = _functions.dict_value_list(_dict)
    params = Parameters()
    for i in range(len(dict_key)):
        param_name = _dict_name_str + dict_key[i]
        params.add(param_name, value=dict_value[i])
    return params


def add_params_from_doct(params, _dict, _dict_name_str):
    dict_key = _functions.dict_key_list(_dict)
    dict_value = _functions.dict_value_list(_dict)
    for i in range(len(dict_key)):
        param_name = _dict_name_str + dict_key[i]
        params.add(param_name, value=dict_value[i])
    return params


def peak_x_gap(params, ideal_x_index, y_data_array):
    # Unpack Parameters:
    parvals = params.valuesdict()
    source_to_detector_cm = parvals['source_to_detector_cm']
    delay_us = parvals['delay_us']
    # Model:
    spectra_path = 'data/spectra.txt'
    range_min = 500
    range_max = 2000
    x_data_array = _functions.get_spectra_range(spectra_path, delay_us,
                                                source_to_detector_cm, range_min, range_max)
    exp_y_index = pku.indexes(y_data_array, thres=0.12/max(y_data_array), min_dist=7)
    exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
    # gap = (exp_x_index[0] - ideal_x_index) ** 2
    # print(exp_x_index)
    # print(ideal_x_index)
    gap = (exp_x_index - ideal_x_index) ** 2
    return gap


def peak_y_gap(params, ideal_x_index, y_data_array):
    # Unpack Parameters:
    parvals = params.valuesdict()
    thick_cm_dict = parvals['thick_cm_dict']
    density_gcm3_dict = parvals['density_gcm3_dict']
    iso_ratio_dicts = parvals['iso_ratio_dicts']

    # Model:


# def peak_x_gap_scipy(delay_us, ideal_x_index, y_data_array):
#     # Unpack Parameters:
#     source_to_detector_cm = 1610.9  # cm
#     # Model:
#     time_lamda_ev_axis = 'eV'
#     spectra_path = 'data/spectra.txt'
#     _slice = 220
#     x_data_array = _functions.get_spectra_slice(spectra_path, time_lamda_ev_axis, delay_us,
#                                                 source_to_detector_cm, _slice)
#     exp_y_index = pku.indexes(y_data_array, thres=0.6, min_dist=50)
#     exp_x_index = pku.interpolate(x_data_array, y_data_array, ind=exp_y_index)
#     gap = (exp_x_index[0] - ideal_x_index) ** 2
#     return gap


# def peak_gap(params, ideal_x_index):
#     # Unpack Parameters:
#     parvals = params.valuesdict()
#     source_to_detector_cm = parvals['source_to_detector_cm']
#     time_us = parvals['time_us']
#     # Model:
#     energy_miliev = 81.787 / (0.3956 * time_us / source_to_detector_cm) ** 2
#     energy_ev = energy_miliev / 1000
#     return (energy_ev - ideal_x_index) ** 2

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


def get_sigma(isotopes, file_names, energy_min, energy_max, sub_x):
    # Transmission calculation of summed and separated contributions by each isotopes
    df_raw = pd.DataFrame()
    df_inter = pd.DataFrame()
    # sigma_iso_ele_isodict = {}
    # sigma_iso_ele_l_isodict = {}
    # thick_cm = thick_mm/10
    # sigma_iso_ele_sum = 0.
    # sigma_iso_ele_l_sum = 0.
    # iso_at_ratio = iso_ratio_list
    sigma_dict = {}
    for i, iso in enumerate(isotopes):
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
        df_inter['E_eV'] = x_energy
        """
        Attention:
        The following part is for producing df_raw of all isotopes for future reference
        """
        # Create a new DataFrame including all isotopic data
        # within the selected energy range
        first_col = iso + ', E_eV'
        second_col = iso + ', Sig_b'
        sigma_dict[iso] = sigma_b
        df.rename(columns={'E_eV': first_col, 'Sig_b': second_col}, inplace=True)
        df_raw = pd.concat([df_raw, df], axis=1)
        df_inter[second_col] = sigma_b

    return x_energy, sigma_dict

# def get_y_dict():
# def get_sigma_term(_input_ele_str, energy_max, energy_min, energy_sub):
#     # Input sample name or names as str, case sensitive
#     # _input_formula = 'AgCo'  # input('Please input the chemicals? ')
#     # _input_thick_mm = 0.025  # float(input('Please input the thickness or majority thickness of stacked foils in mm : '))
#     # _input_thick_cm = _input_thick_mm / 10
#     _database = 'ENDF_VIII'
#     sub_x = energy_sub * (energy_max - energy_min)  # steps used to interpolate database
#
#     # '''Input for dict modification in certain cases: '''
#     # # Thickness input:
#     # special_thick_boo = 'N'
#     # special_thick_element_str = str
#     # special_thick_mm_list = []
#     # special_thick_cm_list = np.array(special_thick_mm_list) / 10
#     # # Enriched isotope ratio input:
#     # enrichment_boo = 'N'  # Isotopic enriched or depleted: Y/N?
#     # enriched_element_str = 'U'
#     # input_ratio_dict = {'U': [0., 0., .15, .85]}
#     # # 'O': [1., 0., 0.]}  #{'233-U': 0., '234-U': 0., '235-U': 0.15, '238-U': 0.85}}
#     # # Special density input:
#     # special_density_boo = 'N'
#     # special_density_element_str = str
#     # special_density_gcm3_list = []
#
#     ''' Parse input formula str and return:
#     (1) elements list, elemental ratio list
#     (2) isotopes dict in the form of {element1: [iso11, iso12, iso13, ...],
#                                       element2: [iso21, iso22, iso23, ...],
#                                       element3: [iso31, iso32, iso33, ...],
#                                       ...}
#     (3) isotopic ratio dict in the form of {element1: {iso11: iso_ratio11, iso12: iso_ratio12, iso13: iso_ratio13, ...},
#                                             element2: {iso21: iso_ratio21, iso22: iso_ratio12, iso23: iso_ratio23, ...},
#                                             element3: {iso31: iso_ratio31, iso32: iso_ratio12, iso33: iso_ratio33, ...},
#                                             ...}
#     '''
#     formula_dict = _functions.input2formula(_input_ele_str)
#     elements = _functions.dict_key_list(formula_dict)
#     ratios = _functions.dict_value_list(formula_dict)
#     sum_ratios = sum(ratios)
#     isotope_dict = _functions.get_isotope_dicts(_database, elements)
#
#     # DICT 1: Thickness dict with option for modification
#     # thick_cm_dict = _functions.repeat_value_dict(elements, _input_thick_cm)
#     # if compound_boo == 'N':
#     #     if special_thick_boo == 'Y':
#     #         thick_cm_dict = _plot_functions.modify_thick_cm_dict_by_input(thick_cm_dict, special_thick_element_str,
#     #                                                                       special_thick_cm_list)
#
#     # # DICT 2: Isotopic mass dict
#     # iso_mass_dicts = _functions.get_iso_mass_dicts_quick(elements, isotope_dict)
#
#     # Dict 3: Molar mass dict
#     molar_mass_dict = _functions.get_molar_mass_dict(elements)
#
#     # DICT 4: Isotope at.% dict with option for modification
#     iso_ratio_dicts = _functions.get_iso_ratio_dicts_quick(elements, isotope_dict)
#
#     # DICT 5: Density dict
#     density_gcm3_dict = _functions.get_density_dict(elements)
#
#     # DICT 6: Stoichiometric ratio
#     # if compound_boo == 'Y':
#     #     # If input is compound, input formula follows the stoichimetric ratios
#     #     ele_at_ratio_dict = {}
#     #     for el in elements:
#     #         ele_at_ratio_dict[el] = formula_dict[el] / sum_ratios
#     # else:
#         # If input is NOT compound, so the input are stack of elements,
#         # stoichimetric ratios need to be calculated based on density and thickness
#     # ele_at_ratio_dict = _functions.ele_ratio_dict(elements, thick_cm_dict, density_gcm3_dict, molar_mass_dict)
#
#     # print('Thickness (cm): ', thick_cm_dict)
#     print('Density (g/cm^3): ', density_gcm3_dict)
#     print('Molar weight (g/mol): ', molar_mass_dict)
#     print('Isotopic ratio (at.%): ', iso_ratio_dicts)
#
#     '''For plotting the database'''
#     sigma_iso_ele_eleisodict = {}  # For transmission calculation at isotope level
#     sigma_iso_ele_sum_eledict = {}  # For transmission calculation at element level
#     # sigma_iso_ele_sum_l_eledict = {}
#     # sigma_iso_ele_l_eleisodict = {}
#     df_raw_dict = {}  # Raw sigma data for elements and isotopes
#
#     for el in elements:
#         # isotopes_list = list(dict.keys(iso_ratio_dicts[el]))
#         iso_ratio_list = list(dict.values(iso_ratio_dicts[el]))
#         # ele_at_ratio = formula_dict[el] / sum_ratios
#
#         # Get sigma related terms
#         file_names = _functions.get_file_path(_database, el)
#         x_energy, sigma_iso_ele_isodict, sigma_iso_ele_sum, df_raw_dict[el] \
#             = _plot_functions.get_xy_from_database(iso_ratio_dicts[el],
#                                                    file_names,
#                                                    energy_min,
#                                                    energy_max,
#                                                    iso_ratio_list,
#                                                    sub_x,
#                                                    )
#         # Two level dict of isotopic array of (L * sigma * iso_ratio * ele_ratio)
#         # sigma_iso_ele_l_eleisodict[el] = sigma_iso_ele_l_isodict
#         # One level dict of elemental array of (L * sigma * iso_ratio * ele_ratio)
#         # sigma_iso_ele_sum_l_eledict[el] = sigma_iso_ele_sum * thick_cm_dict[el]
#
#         # Two level dict of isotopic array of (sigma * iso_ratio * ele_ratio)
#         sigma_iso_ele_eleisodict[el] = sigma_iso_ele_isodict
#         # One level dict of elemental array of (sigma * iso_ratio * ele_ratio)
#         sigma_iso_ele_sum_eledict[el] = sigma_iso_ele_sum
#
#     return elements, isotope_dict, density_gcm3_dict, molar_mass_dict, x_energy, sigma_iso_ele_eleisodict, sigma_iso_ele_sum_eledict
