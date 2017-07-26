import _functions
import numpy as np
import pandas as pd
import periodictable as pt
from periodictable.constants import avogadro_number
import matplotlib.pyplot as plt
from matplotlib import style
from scipy.interpolate import *
style.use('ggplot')


def get_xy_from_database(isotopes, file_names, energy_min, energy_max, iso_ratio_list, sub_x, ele_at_ratio):
    # Transmission calculation of summed and separated contributions by each isotopes
    df = pd.DataFrame()
    df_raw = pd.DataFrame()
    sigma_iso_ele_isodict = {}
    # sigma_iso_ele_l_isodict = {}
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
        # sigma_iso_ele_l_isodict[iso] = sigma_iso_ele_isodict[iso] * thick_cm
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

    return x_energy, sigma_iso_ele_isodict, sigma_iso_ele_sum, df_raw


def atoms_per_cm3(elements, thick_cm_dict, density_gcm3_dict, molar_mass_dict):
    n = 0.
    for ele in elements:
        l_x_n = l_x_n + thick_cm_dict[ele] * density_gcm3_dict[ele] / molar_mass_dict[ele]
    l_n_avo = l_x_n * pt.constants.avogadro_number
    return atoms_per_cm3


def l_x_n_multi_ele_stack(elements, thick_cm_dict, density_gcm3_dict, molar_mass_dict):
    l_x_n = 0.
    for ele in elements:
        l_x_n = l_x_n + thick_cm_dict[ele] * density_gcm3_dict[ele] / molar_mass_dict[ele]
    l_n_avo = l_x_n * pt.constants.avogadro_number
    print('Thickness(l) x atoms_per_cm^3(N) : ', l_n_avo)
    return l_n_avo


def l_x_n_compound(elements, thick_cm, compound_density, molar_mass_dict, formula_dict, sum_ratios):
    molar_mass_sum = 0.
    for ele in elements:
        molar_mass_sum = molar_mass_sum + molar_mass_dict[ele] * formula_dict[ele] / sum_ratios
    l_x_n = thick_cm * compound_density / molar_mass_sum
    l_n_avo = l_x_n * pt.constants.avogadro_number
    print('Thickness(l) x atoms_per_cm^3(N) : ', l_n_avo)
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


def get_tot_trans_for_single_ele(_input_ele_str, _input_thick_mm, energy_max, energy_min, energy_sub):
    # Input sample name or names as str, case sensitive
    # _input_formula = 'AgCo'  # input('Please input the chemicals? ')
    # _input_thick_mm = 0.025  # float(input('Please input the thickness or majority thickness of stacked foils in mm : '))
    _input_thick_cm = _input_thick_mm / 10
    _database = 'ENDF_VIII'
    # energy_max = 300  # max incident energy in eV
    # energy_min = 0  # min incident energy in eV
    # energy_sub = 100  # steps used to interpolate database
    sub_x = energy_sub * (energy_max - energy_min)  # steps used to interpolate database
    # compound_boo = 'N'  # Compound or single/multi elements foil/stacked foils: Y/N?

    # '''Input for dict modification in certain cases: '''
    # # Thickness input:
    # special_thick_boo = 'N'
    # special_thick_element_str = str
    # special_thick_mm_list = []
    # special_thick_cm_list = np.array(special_thick_mm_list) / 10
    # # Enriched isotope ratio input:
    # enrichment_boo = 'N'  # Isotopic enriched or depleted: Y/N?
    # enriched_element_str = 'U'
    # input_ratio_dict = {'U': [0., 0., .15, .85]}
    # # 'O': [1., 0., 0.]}  #{'233-U': 0., '234-U': 0., '235-U': 0.15, '238-U': 0.85}}
    # # Special density input:
    # special_density_boo = 'N'
    # special_density_element_str = str
    # special_density_gcm3_list = []

    # '''How you want the data to be plotted?'''
    # _plot_or_not = 'Y'
    # _energy_x_axis = 'Y'  # 1 means plot x-axis as energy in eV
    # _trans_y_axis = 'N'  # 1 means plot y-axis as transmission
    # _plot_each_ele_contribution = 'Y'  # 1 means plot each element's contribution
    # _plot_each_iso_contribution = 'N'  # 1 means plot each isotope's contribution
    # _plot_mixed = 'Y'  # 1 means plot mixed resonance
    # '''Export to clipboard for Excel or DataGraph?'''
    # _export_to_clipboard_boo = 'N'

    ''' Parse input formula str and return:
    (1) elements list, elemental ratio list
    (2) isotopes dict in the form of {element1: [iso11, iso12, iso13, ...], 
                                      element2: [iso21, iso22, iso23, ...], 
                                      element3: [iso31, iso32, iso33, ...], 
                                      ...}
    (3) isotopic ratio dict in the form of {element1: {iso11: iso_ratio11, iso12: iso_ratio12, iso13: iso_ratio13, ...},
                                            element2: {iso21: iso_ratio21, iso22: iso_ratio12, iso23: iso_ratio23, ...},
                                            element3: {iso31: iso_ratio31, iso32: iso_ratio12, iso33: iso_ratio33, ...},
                                            ...}
    '''
    formula_dict = _functions.input2formula(_input_ele_str)
    elements = _functions.dict_key_list(formula_dict)
    ratios = _functions.dict_value_list(formula_dict)
    sum_ratios = sum(ratios)
    isotope_dict = _functions.get_isotope_dicts(_database, elements)

    # DICT 1: Thickness dict with option for modification
    thick_cm_dict = _functions.repeat_value_dict(elements, _input_thick_cm)
    # if compound_boo == 'N':
    #     if special_thick_boo == 'Y':
    #         thick_cm_dict = _plot_functions.modify_thick_cm_dict_by_input(thick_cm_dict, special_thick_element_str,
    #                                                                       special_thick_cm_list)

    # DICT 2: Isotopic mass dict
    iso_mass_dicts = _functions.get_iso_mass_dicts_quick(elements, isotope_dict)

    # Dict 3: Molar mass dict
    molar_mass_dict = _functions.get_molar_mass_dict(elements)

    # DICT 4: Isotope at.% dict with option for modification
    iso_ratio_dicts = _functions.get_iso_ratio_dicts_quick(elements, isotope_dict)

    # DICT 5: Density dict
    density_gcm3_dict = _functions.get_density_dict(elements)

    # # Update DICT 3 & 4 & 5: isotopic ratio changes lead to |Density| & |Molar mass| changes
    # if enrichment_boo == 'Y':
    #     # Update isotope at.% ratio dict
    #     iso_ratio_dicts, enriched_element = _plot_functions.modify_iso_ratio_dicts(elements, isotope_dict,
    #                                                                                enriched_element_str,
    #                                                                                input_ratio_dict)
    #     # Update molar mass dict
    #     molar_mass_dict = _plot_functions.modify_molar_mass_dict_by_enrichment(molar_mass_dict, enriched_element,
    #                                                                            isotope_dict, iso_ratio_dicts,
    #                                                                            iso_mass_dicts)
    #     # Update density dict
    #     density_gcm3_dict = _plot_functions.modify_density_dict_by_enrichment(density_gcm3_dict, enriched_element,
    #                                                                           isotope_dict, iso_ratio_dicts)
    #
    # # Update DICT 5: Density dict, if special case encountered
    # if compound_boo == 'N':
    #     if special_density_boo == 'Y':
    #         # Stacked foils and would like to modify density for specific element
    #         density_gcm3_dict = _plot_functions.modify_density_dict_by_input(density_gcm3_dict,
    #                                                                          special_density_element_str,
    #                                                                          special_density_gcm3_list)
    # else:
    #     if special_density_boo == 'Y':
    #         # Not isolated elements or mixture or compound need density input currently
    #         input_tot_density = 0.7875

    print('Thickness (cm): ', thick_cm_dict)
    print('Density (g/cm^3): ', density_gcm3_dict)
    print('Isotopic ratio (at.%)', iso_ratio_dicts)
    print('Molar weight (g/mol): ', molar_mass_dict)

    '''For plotting the database'''
    sigma_iso_ele_eleisodict = {}  # For transmission calculation at isotope level
    sigma_iso_ele_sum_eledict = {}  # For transmission calculation at element level
    # sigma_iso_ele_sum_l_eledict = {}
    # sigma_iso_ele_l_eleisodict = {}
    df_raw_dict = {}  # Raw sigma data for elements and isotopes

    for el in elements:
        # isotopes_list = list(dict.keys(iso_ratio_dicts[el]))
        iso_ratio_list = list(dict.values(iso_ratio_dicts[el]))
        # iso_ratio_array = np.array(iso_ratio_list)
        # iso_mass_list = list(dict.values(iso_mass_dicts[el]))
        # iso_mass_array = np.array(iso_mass_list)
        ele_at_ratio = formula_dict[el] / sum_ratios

        # Get sigma related terms
        file_names = _functions.get_file_path(_database, el)
        x_energy, sigma_iso_ele_isodict, sigma_iso_ele_sum, df_raw_dict[el] \
            = get_xy_from_database(iso_ratio_dicts[el],
                                   file_names,
                                   energy_min,
                                   energy_max,
                                   iso_ratio_list,
                                   sub_x,
                                   ele_at_ratio)
        # Two level dict of isotopic array of (L * sigma * iso_ratio * ele_ratio)
        # sigma_iso_ele_l_eleisodict[el] = sigma_iso_ele_l_isodict
        # One level dict of elemental array of (L * sigma * iso_ratio * ele_ratio)
        # sigma_iso_ele_sum_l_eledict[el] = sigma_iso_ele_sum * thick_cm_dict[el]

        # Two level dict of isotopic array of (sigma * iso_ratio * ele_ratio)
        sigma_iso_ele_eleisodict[el] = sigma_iso_ele_isodict
        # One level dict of elemental array of (sigma * iso_ratio * ele_ratio)
        sigma_iso_ele_sum_eledict[el] = sigma_iso_ele_sum

    # Get Thickness * number of atoms per cm^3
    # if compound_boo == 'N':
        # Stacked foils or single foil
    mixed_l_n_avo = l_x_n_multi_ele_stack(elements,
                                          thick_cm_dict,
                                          density_gcm3_dict,
                                          molar_mass_dict)
    # else:
    #     # For compound
    #     thick_cm_list = list(dict.values(thick_cm_dict))
    #     thick_cm = thick_cm_list[0]
    #     compound_density = input_tot_density
    #     mixed_l_n_avo = _plot_functions.l_x_n_compound(elements,
    #                                                    thick_cm,
    #                                                    compound_density,
    #                                                    molar_mass_dict,
    #                                                    formula_dict,
    #                                                    sum_ratios)

    # Get the tot transmission for all
    # yi_values_l = list(dict.values(sigma_iso_ele_sum_l_eledict))
    # yi_values_l_sum = sum(yi_values_l)
    # # sum of (sigma * ele_ratio * iso_ratio * l)
    yi_values = list(dict.values(sigma_iso_ele_sum_eledict))
    yi_values_sum = sum(yi_values)
    # sum of (sigma * ele_ratio * iso_ratio)
    # print(yi_values)
    y_trans_tot = _functions.sig_l_2trans_quick(mixed_l_n_avo, yi_values_sum)

    return x_energy, y_trans_tot
