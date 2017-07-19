import _functions
import numpy as np
import pandas as pd
import periodictable as pt
import matplotlib.pyplot as plt
from matplotlib import style
from scipy.interpolate import *
style.use('ggplot')


def get_xy_from_database(isotopes, thick_cm, file_names, energy_min, energy_max, iso_ratio_list, sub_x, ele_at_ratio):
    # Transmission calculation of summed and separated contributions by each isotopes
    df = pd.DataFrame()
    df_raw = pd.DataFrame()
    sigma_iso_ele_isodict = {}
    sigma_iso_ele_l_isodict = {}
    # thick_cm = thick_mm/10
    sigma_iso_ele_sum = 0.
    # sigma_iso_ele_l_sum = 0.
    iso_at_ratio = iso_ratio_list
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
        # y_i_sum = y_i_sum + y_i * iso_abundance[i] * ele_at_ratio
        sigma_iso_ele_isodict[iso] = sigma_b * iso_at_ratio[i] * ele_at_ratio
        sigma_iso_ele_l_isodict[iso] = sigma_iso_ele_isodict[iso] * thick_cm
        sigma_iso_ele_sum = sigma_iso_ele_sum + sigma_b * iso_at_ratio[i] * ele_at_ratio
        # sigma_iso_ele_l_sum = sigma_iso_ele_l_sum + sigma_b * iso_at_ratio[i] * ele_at_ratio * thick_cm

        """
        Attention:
        The following part is for producing df_raw of all isotopes for future reference
        """
        # Create a new DataFrame including all isotopic data
        # within the selected energy range
        first_col = iso + ', E_eV'
        second_col = iso + ', Sig_b'
        df.rename(columns={'E_eV': first_col, 'Sig_b': second_col}, inplace=True)
        df_raw = pd.concat([df_raw, df], axis=1)

    return x_energy, sigma_iso_ele_isodict, sigma_iso_ele_l_isodict, sigma_iso_ele_sum, df_raw


def l_x_n_multi_ele_stack(elements, thick_cm_dict, density_gcm3_dict, molar_mass_dict):
    l_x_n = 0.
    for ele in elements:
        l_x_n = l_x_n + density_gcm3_dict[ele] * thick_cm_dict[ele] / molar_mass_dict[ele]
    l_n_avo = l_x_n * pt.constants.avogadro_number
    print('Thickness(l) x atoms_per_cm^3(N) (g/cm^3): ', l_x_n)
    return l_n_avo


def plot_database(_energy_x_axis, _trans_y_axis, _plot_mixed, _plot_each_ele_contribution, _plot_each_iso_contribution,
            elements, isotope_dict, x_energy, y_trans_tot, y_ele_dict, y_iso_dicts, _input_formula):
    ### Plot the theoretical neutron resonance
    ### Determine x y axis types and captions
    if _energy_x_axis == 'Y':
        _x_axis = x_energy
        _x_words = 'Energy (eV)'
    else:
        _x_axis = _functions.ev2lamda(x_energy)
        _x_words = 'Wavelength (Å)'

    if _trans_y_axis == 'Y':
        _y_words = 'Neutron transmission'
    else:
        _y_words = 'Neutron attenuation'

    ### Determine x y axis values
    if _plot_mixed == 'Y':
        if _trans_y_axis == 'Y':
            _y_axis = y_trans_tot
        else:
            _y_axis = 1 - y_trans_tot
        plt.plot(_x_axis, _y_axis, label=_input_formula)

    if _plot_each_ele_contribution == 'Y':
        for _ele in elements:
            _y_each_axis = y_ele_dict[_ele]
            plt.plot(_x_axis, _y_each_axis, label=_ele)

    if _plot_each_iso_contribution == 'Y':
        for _ele in elements:
            for _iso in isotope_dict[_ele]:
                _y_each_axis = y_iso_dicts[_ele][_iso]
                plt.plot(_x_axis, _y_each_axis, label=_iso)

    plt.ylim(-0.01, 1.01)
    plt.xlabel(_x_words)
    plt.ylabel(_y_words)
    plt.legend(loc='best')
    plt.show()


### Functions to modify the dicts in special case

def modify_thick_cm_dict_by_input(thick_cm_dict, special_thick_element_str, special_thick_cm_list):
    # For elements with various thickness:
    special_element = special_thick_element_str.split(' ')
    thick_cm_dict = _functions.dict_replace_value_by_key(thick_cm_dict, special_element, special_thick_cm_list)
    print('Modified thickness by input (cm): ', thick_cm_dict)
    return thick_cm_dict


def modify_density_dict_by_input(density_gcm3_dict, special_element_str, special_density_gcm3_list):
    # For elements with special density:
    special_element = special_element_str.split(' ')
    density_gcm3_dict = _functions.dict_replace_value_by_key(density_gcm3_dict, special_element, special_density_gcm3_list)
    print('Modified density by input (g/cm^3): ', density_gcm3_dict)
    return density_gcm3_dict


def modify_iso_ratio_dicts(elements, isotope_dict, enriched_element_str, input_ratio_dict):
    iso_ratio_dicts = _functions.get_iso_ratio_dicts_quick(elements, isotope_dict)
    enriched_element = enriched_element_str.split(' ')
    for ele in enriched_element:
        isotopes = isotope_dict[ele]
        iso_ratio_dicts[ele] = _functions.dict_replace_value_by_key(iso_ratio_dicts[ele], isotopes, input_ratio_dict[ele])
    return iso_ratio_dicts, enriched_element


def modify_molar_mass_dict_by_enrichment(molar_mass_dict, enriched_element, isotope_dict, enriched_iso_ratio_dicts, iso_mass_dicts):
    for ele in enriched_element:
        molar_mass = 0.
        for iso in isotope_dict[ele]:
            molar_mass = molar_mass + enriched_iso_ratio_dicts[ele][iso] * iso_mass_dicts[ele][iso]
        molar_mass_dict[ele] = molar_mass
    print('Modified molar mass by enrichment (g/mol): ', molar_mass_dict)
    return molar_mass_dict


def modify_density_dict_by_enrichment(density_gcm3_dict, enriched_element, isotope_dict, enriched_iso_ratio_dicts):
    for ele in enriched_element:
        density_gcm3 = 0.
        for iso in isotope_dict[ele]:
            density_gcm3 = density_gcm3 + enriched_iso_ratio_dicts[ele][iso] * pt.elements.isotope(iso).density
        density_gcm3_dict[ele] = density_gcm3
    print('Modified density by enrichment (g/cm^3): ', density_gcm3_dict)
    return density_gcm3_dict

