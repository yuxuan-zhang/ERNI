import numpy as np
import periodictable as pt
import re

# Functions for energy, lamda and time conversions


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


def input2formula(_input, _natural_ele):
    _input_parsed = re.findall(r'([A-Z][a-z]*)(\d*)', _input)
    _formula = {}
    _natural_mix = {}
    _ratio_array = {}
    for _element in _input_parsed:
        _element_list = list(_element)
        if _element_list[1] == '':
            _element_list[1] = 1
        _element_list[1] = int(_element_list[1])
        _formula[_element_list[0]] = _element_list[1]
        _natural_mix[_element_list[0]] = _natural_ele
        _ratio_array[_element_list[0]] = []
    print('Parsed chemical formula: {}'.format(_formula))
    return _formula, _natural_mix, _ratio_array
