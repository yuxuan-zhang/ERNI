import numpy as np
import periodictable as pt
from periodictable import constants
import re
import os
import glob
import pandas as pd
from scipy import interpolate

'''
Energy, wavelength and time conversions
'''
def ev2lamda(energy_ev):  # function to convert energy in eV to angstrom
    energy_miliev = energy_ev * 1000
    lamda = np.sqrt(81.787/energy_miliev)
    return lamda


def time2lamda(time_tot_s, delay_us, source_to_detector_cm):  # function to convert time in us to angstrom
    time_tot_us = 1e6 * time_tot_s + delay_us
    lamda = 0.3956 * time_tot_us/source_to_detector_cm
    return lamda


def lamda2ev(lamda):  # function to convert angstrom to eV
    energy_miliev = 81.787/(lamda ** 2)
    energy_ev = energy_miliev/1000
    return energy_ev


def time2ev(time_tot_s, delay_us, source_to_detector_cm):  # function to convert time in us to energy in eV
    time_tot_us = 1e6 * time_tot_s + delay_us
    energy_miliev = 81.787/(0.3956 * time_tot_us/source_to_detector_cm) ** 2
    energy_ev = energy_miliev/1000
    return energy_ev


def atoms_per_cm3(density, mass):
    n_atoms = density * pt.constants.avogadro_number/mass
    print('Number of atoms per unit volume (#/cm^3): {}'.format(n_atoms))
    return n_atoms


def sig2trans(_thick_cm, _atoms_per_cm3, _ele_atomic_ratio, _sigma_b, _iso_atomic_ratio):
    neutron_transmission = np.exp(-1 * _thick_cm * _atoms_per_cm3 *
                                  _ele_atomic_ratio * _sigma_b * 1e-24 * _iso_atomic_ratio)
    return neutron_transmission


def sig2trans_quick(_thick_cm, _atoms_per_cm3, _sigma_portion_sum):
    neutron_transmission = np.exp(-1 * _thick_cm * _atoms_per_cm3 * 1e-24 * _sigma_portion_sum)
    return neutron_transmission


def sig_l_2trans_quick(l_n_avo, sigma_portion_sum):
    neutron_transmission = np.exp(-1 * l_n_avo * 1e-24 * sigma_portion_sum)
    return neutron_transmission


def get_isotope_dicts(_database, _element):
    # main_dir = os.path.dirname(os.path.abspath(__file__))
    isotope_dict = {}
    for _each in _element:
        path = 'data_web/' + _database + '/' + str(_each) + '*.csv'
        # path = main_dir + '/data_web/' + _database + '/' + _each + '*.csv'
        file_names = glob.glob(path)
        isotope_dict_mirror = {}
        for _i, file in enumerate(file_names):
            # Obtain element, z number from the basename
            _basename = os.path.basename(file)
            _name_number_csv = _basename.split('.')
            _name_number = _name_number_csv[0]
            _name = _name_number.split('-')
            _symbol = _name[1] + '-' + _name[0]
            isotope = str(_symbol)
            isotope_dict_mirror[isotope] = isotope
        isotopes = list(dict.values(isotope_dict_mirror))
        isotope_dict[_each] = isotopes
    return isotope_dict


def input2formula(_input):
    _input_parsed = re.findall(r'([A-Z][a-z]*)(\d*)', _input)
    _formula = {}
    # _natural_ele_boo_dict = {}
    # _natural_mix = {}
    # _ratio_array = {}
    for _element in _input_parsed:
        _element_list = list(_element)
        if _element_list[1] == '':
            _element_list[1] = 1
        _element_list[1] = int(_element_list[1])
        _formula[_element_list[0]] = _element_list[1]
        # _natural_ele_boo_dict[_element_list[0]] = 'Y'
    print('Parsed chemical formula: {}'.format(_formula))
    return _formula


def dict_key_list(_formula_dict):
    _keys = list(dict.keys(_formula_dict))
    return _keys


def dict_value_list(_formula_dict):
    _values = list(dict.values(_formula_dict))
    return _values


def boo_dict(_key_list, _y_or_n):
    _boo_dict = {}
    for key in _key_list:
        _boo_dict[key] = _y_or_n
    return _boo_dict

#########
def get_thick_dict(_key_list, _thick_mm):
    _thick_dict = {}
    for key in _key_list:
        _thick_dict[key] = _thick_mm
    return _thick_dict


def get_density_dict(_key_list):
    _density_dict = {}
    for key in _key_list:
        _density_dict[key] = pt.elements.isotope(key).density
        # key can be element ('Co') or isotopes ('59-Co')
        # Unit: g/cm3
    return _density_dict


####
def get_molar_mass_dict(elements):
    molar_mass_dict = {}
    for ele in elements:
        molar_mass_dict[ele] = pt.elements.isotope(ele).mass
    return molar_mass_dict


def get_iso_ratio_dict(isotopes):
    # natural_density = pt.elements.isotope(_element).density
    iso_ratio_dict = {}
    for iso in isotopes:
        iso_ratio_dict[iso] = pt.elements.isotope(iso).abundance
    return iso_ratio_dict


def get_iso_ratio_dicts(elements, iso_ratio_dict):
    # natural_density = pt.elements.isotope(_element).density
    iso_ratio_dicts = {}
    for el in elements:
        iso_ratio_dicts[el] = iso_ratio_dict
    return iso_ratio_dicts


def get_iso_ratio_dicts_quick(elements, isotope_dict):
    # natural_density = pt.elements.isotope(_element).density
    iso_ratio_dicts = {}
    for el in elements:
        iso_ratio_dict = {}
        for iso in isotope_dict[el]:
            iso_ratio_dict[iso] = pt.elements.isotope(iso).abundance/100
        iso_ratio_dicts[el] = iso_ratio_dict
    return iso_ratio_dicts


def get_iso_mass_dicts_quick(elements, isotope_dict):
    # natural_density = pt.elements.isotope(_element).density
    iso_mass_dicts = {}
    for el in elements:
        iso_mass_dict = {}
        for iso in isotope_dict[el]:
            iso_mass_dict[iso] = pt.elements.isotope(iso).mass
        iso_mass_dicts[el] = iso_mass_dict
    return iso_mass_dicts


def empty_2d_dict(top_level_name, top_base_dict):
    whole_dict = {}
    base_dict = {}
    for top in top_level_name:
        for base in top_base_dict[top]:
            base_dict[base] = 1
        whole_dict[top] = base_dict
    return whole_dict


def create_2d_dict(top_level_name, top_base_dict, value_dict):
    whole_dict = {}
    base_dict = {}
    for top in top_level_name:
        for base in top_base_dict[top]:
            base_dict[base] = value_dict[top][base]
        whole_dict[top] = base_dict
    return whole_dict

####
#########


def get_file_path(_database, _element):
    path = 'data_web/' + _database + '/' + _element + '*.csv'
    file_names = glob.glob(path)
    return file_names


def empty_dict(_key_list):
    _empty_dicts = {}
    _empty_dict = {}
    for key in _key_list:
        _empty_dicts[key] = _empty_dict
    return _empty_dicts


def repeat_value_dict(_key_list, value):
    _thick_dict = {}
    for key in _key_list:
        _thick_dict[key] = value
    return _thick_dict


def dict_value_by_key(_key_list, _value_list):
    p = 0
    _dict = {}
    for key in _key_list:
        _dict[key] = _value_list[p]
        p = p + 1
    return _dict


def ele_ratio_dict(element_list, thick_cm_dict, density_gcm3_dict, molar_mass_dict):
    mol_dict = {}
    ele_at_ratio_dict = {}
    mol_sum = 0.
    for ele in element_list:
        mol_dict[ele] = density_gcm3_dict[ele] * thick_cm_dict[ele] / molar_mass_dict[ele]
        mol_sum = mol_sum + mol_dict[ele]
    for ele in element_list:
        ele_at_ratio_dict[ele] = mol_dict[ele] / mol_sum
    return ele_at_ratio_dict


def boo_dict_invert_by_key(_key_list, _boo_dict):
    for key in _key_list:
        if _boo_dict[key] == 'Y':
            _boo_dict[key] = 'N'
        else:
            _boo_dict[key] = 'Y'
    return _boo_dict


def dict_replace_value_by_key(_dict, _key_list, _value_list):
    p = 0
    for key in _key_list:
        _dict[key] = _value_list[p]
        p = p + 1
    return _dict


def formula_ratio_array(_input, _all_ele_boo_dict, ratios_dict):
    _natural_ele = {}
    _ratio_array = {}
    for _element in _input:
        _natural_ele[_element] = _all_ele_boo_dict[_element]
        if _all_ele_boo_dict[_element] == 'Y':
            _ratio_array[_element] = []
        else:
            _ratio_array[_element] = ratios_dict[_element]
    print('Natual elements? ', _natural_ele)
    print('Isotope ratio array: ', _ratio_array)
    return _ratio_array


def get_normalized_data(_filename):
    df = pd.read_csv(_filename, header=None, skiprows=1)
    data_array = np.array(df[1])
    data = data_array[:int(len(data_array)/2)]
    ob = data_array[int(len(data_array)/2):]
    normalized_array = data/ob
    # OB at the end of 2773
    return normalized_array


def get_normalized_data_slice(_filename, _ignore):
    df = pd.read_csv(_filename, header=None, skiprows=1)
    data_array = np.array(df[1])
    data = data_array[:int(len(data_array)/2)]
    ob = data_array[int(len(data_array)/2):]
    normalized_array = data/ob
    normalized_array = normalized_array[_ignore:]
    # OB at the end of 2773
    return normalized_array


def get_normalized_data_range(_filename, range_min, range_max):
    df = pd.read_csv(_filename, header=None, skiprows=1)
    data_array = np.array(df[1])
    data = data_array[:int(len(data_array)/2)]
    ob = data_array[int(len(data_array)/2):]
    normalized_array = data/ob
    normalized_array = normalized_array[range_min:range_max]
    normalized_array = normalized_array[::-1]  # Flip array from descending to normal
    # OB at the end of 2773
    return normalized_array


def get_ob_range(_filename, range_min, range_max):
    df = pd.read_csv(_filename, header=None, skiprows=1)
    data_array = np.array(df[1])
    data = data_array[:int(len(data_array)/2)]
    ob = data_array[int(len(data_array)/2):]
    ob = ob[range_min:range_max]
    ob = ob[::-1]
    # normalized_array = data/ob
    # normalized_array = normalized_array[range_min:range_max]
    # normalized_array = normalized_array[::-1]  # Flip array from descending to normal
    # OB at the end of 2773
    return ob


def get_spectra_range(_filename, delay_us, source_to_detector_cm, range_min, range_max, time_lamda_ev_axis='eV'):
    """
    Get spectra file and convert time to wavelength or energy.
    :param _filename:
    :param delay_us:
    :param source_to_detector_cm:
    :param range_min:
    :param range_max:
    :param time_lamda_ev_axis:
    :return:
    """
    df_spectra = pd.read_csv(_filename, sep='\t', header=None)
    time_array = (np.array(df_spectra[0]))
    # flux_array = (np.array(df_spectra[1]))
    if time_lamda_ev_axis == 'lamda':
        lamda_array = time2lamda(time_array, delay_us, source_to_detector_cm)
        return lamda_array
    if time_lamda_ev_axis == 'eV':
        ev_array = time2ev(time_array, delay_us, source_to_detector_cm)
        ev_array = ev_array[range_min:range_max]
        ev_array = ev_array[::-1]  # Flip array from descending to normal
        return ev_array
    if time_lamda_ev_axis == 'time':
        time_array = time_array[range_min:range_max]
        return time_array


def get_spectra(_filename, delay_us, source_to_detector_cm, time_lamda_ev_axis='eV'):
    df_spectra = pd.read_csv(_filename, sep='\t', header=None)
    time_array = (np.array(df_spectra[0]))
    # flux_array = (np.array(df_spectra[1]))
    if time_lamda_ev_axis == 'lamda':
        lamda_array = time2lamda(time_array, delay_us, source_to_detector_cm)
        return lamda_array
    if time_lamda_ev_axis == 'eV':
        ev_array = time2ev(time_array, delay_us, source_to_detector_cm)
        return ev_array
    if time_lamda_ev_axis == 'time':
        return time_array


def get_spectra_slice(_filename, time_lamda_ev_axis, delay_us, source_to_detector_cm, _slice):
    df_spectra = pd.read_csv(_filename, sep='\t', header=None)
    time_array = (np.array(df_spectra[0]))
    # flux_array = (np.array(df_spectra[1]))
    if time_lamda_ev_axis == 'lamda':
        lamda_array = time2lamda(time_array, delay_us, source_to_detector_cm)
        return lamda_array
    if time_lamda_ev_axis == 'eV':
        ev_array = time2ev(time_array, delay_us, source_to_detector_cm)
        ev_array = ev_array[_slice:]
        return ev_array
    if time_lamda_ev_axis == 'lamda':
        return time_array


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

