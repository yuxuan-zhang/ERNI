import numpy as np
import periodictable as pt
import re
import os
import glob


def ev2lamda(energy):  # function to convert energy in eV to angstrom
    lamda = np.sqrt(81.787/(1000 * energy))
    return lamda


def time2lamda(time_tot):  # function to convert time in us to angstrom
    lamda = 0.3956 * time_tot/source_to_detector
    return lamda


def lamda2ev(lamda):  # function to convert angstrom to eV
    energy = 81.787/(1000 * lamda ** 2)
    return energy


def time2ev(time_tot, source_to_detector):  # function to convert time in us to energy in eV
    energy = 81.787/(1000 * (0.3956 * time_tot/source_to_detector) ** 2)
    return energy


def atoms_per_cm3(density, mass):
    n_atoms = density * pt.constants.avogadro_number/mass
    print('Number of atoms per unit volume (#/cm^3): {}'.format(n_atoms))
    return n_atoms


def sig2trans(_thick_cm, _atoms_per_cm3, _ele_atomic_ratio, _sigma_b, _iso_atomic_ratio):
    neutron_transmission = np.exp(-1 * _thick_cm * _atoms_per_cm3 *
                                  _ele_atomic_ratio * _sigma_b * 1e-24 * _iso_atomic_ratio)
    return neutron_transmission


def sig2trans_quick(_thick_mm, _atoms_per_cm3, _sigma_portion_sum):
    _thick_cm = _thick_mm/10
    neutron_transmission = np.exp(-1 * _thick_cm * _atoms_per_cm3 * 1e-24 * _sigma_portion_sum)
    return neutron_transmission


def get_isotope_dict(_database, _element):
    main_dir = os.path.dirname(os.path.abspath(__file__))
    isotope_dicts = {}
    for _each in _element:
        path = main_dir + '/data_web/' + _database + '/' + _each + '*.csv'
        file_names = glob.glob(path)
        isotope_dict = {}
        for _i, file in enumerate(file_names):
            # Obtain element, z number from the basename
            _basename = os.path.basename(file)
            _name_number_csv = _basename.split('.')
            _name_number = _name_number_csv[0]
            _name = _name_number.split('-')
            _symbol = _name[1] + '-' + _name[0]
            isotope = str(_symbol)
            isotope_dict[isotope] = isotope
        isotopes = list(dict.values(isotope_dict))
        isotope_dicts[_each] = isotopes
    return isotope_dicts


# def get_abundance_dicts(_isotope_dicts, _element):
#     abundance_dict = {}
#     abundance_dicts = {}
#     for _each in _element:
#         isotopes = list(dict.values(_isotope_dicts[_each])
#         for _iso in isotopes:
#             abundance_dict[_iso] = pt.elements.isotope(_iso).abundance / 100
#         abundance_dicts[_each] = abundance_dict
#     return abundance_dicts


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
    return _formula #, _natural_ele_boo_dict


def dict_key_list(_formula_dict):
    _elements = list(dict.keys(_formula_dict))
    return _elements


def dict_value_list(_formula_dict):
    _ratios = list(dict.values(_formula_dict))
    return _ratios


def boo_dict(_key_list):
    _boo_dict = {}
    for key in _key_list:
        _boo_dict[key] = 'Y'
    return _boo_dict


def thick_dict(_key_list, _thick_mm):
    _thick_dict = {}
    for key in _key_list:
        _thick_dict[key] = _thick_mm
    return _thick_dict


def empty_dict(_key_list):
    _empty_dicts = {}
    _empty_dict = {}
    for key in _key_list:
        _empty_dicts[key] = _empty_dict
    return _empty_dicts


def boo_dict_invert_by_key(_key_list, _boo_dict):
    for key in _key_list:
        if _boo_dict[key] == 'Y':
            _boo_dict[key] = 'N'
        else:
            _boo_dict[key] = 'Y'
    return _boo_dict


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

# def get_

        # def deter_xy(_energy_x_axis, ):
#
# if _energy_x_axis == 'Y':
#     _x_axis = x_energy
#     _x_words = 'Energy (eV)'
# else:
#     _x_axis = ev2lamda(x_energy)
#     _x_words = 'Wavelength (Å)'
#
# if _trans_y_axis == 'Y':
#     _y_axis = y_trans_tot
#     _y_words = 'Neutron transmission'
# else:
#     _y_axis = 1 - y_trans_tot
#     _y_words = 'Neutron attenuation'