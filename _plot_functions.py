import _functions
import numpy as np
import pandas as pd
import periodictable as pt
import matplotlib.pyplot as plt
from matplotlib import style
import glob
import os
from scipy import interpolate
style.use('ggplot')


def get_pre_data(_database, _element):
    main_dir = os.path.dirname(os.path.abspath(__file__))
    path = main_dir + '/data_web/' + _database + '/' + _element + '*.csv'
    file_names = glob.glob(path)
    abundance_dict = {}
    density_dict = {}
    mass_dict = {}
    # z_number = {}
    for _i, file in enumerate(file_names):
        # Obtain element, z number from the basename
        _basename = os.path.basename(file)
        _name_number_csv = _basename.split('.')
        _name_number = _name_number_csv[0]
        _name = _name_number.split('-')
        _symbol = _name[1] + '-' + _name[0]
        abundance_dict[str(_symbol)] = pt.elements.isotope(_symbol).abundance / 100
        density_dict[str(_symbol)] = pt.elements.isotope(_symbol).density
        mass_dict[str(_symbol)] = pt.elements.isotope(_symbol).mass
    isotopes = list(dict.keys(abundance_dict))  # List of isotopes such as '238-U', ''235-U
    iso_abundance = list(dict.values(abundance_dict))  # List of isotopic abundance
    iso_density = list(dict.values(density_dict))  # List of isotopic density
    iso_mass = list(dict.values(mass_dict))  # List of isotopic molar mass
    return isotopes, abundance_dict, iso_abundance, iso_density, iso_mass, abundance_dict, density_dict, mass_dict, file_names


def get_mass_iso_ele(iso_abundance, iso_mass, ele_at_ratio, _natural_mix, _unnatural_ratio_dict):
    # Calculate the number of atoms per unit volume
    abundance_array = np.array(iso_abundance)
    mass_array = np.array(iso_mass)
    if _natural_mix == 'Y':
        abundance_array = abundance_array
    else:
        ratio_list = list(dict.values(_unnatural_ratio_dict))
        ratio_array = np.array(ratio_list)
        abundance_array = ratio_array
    mass_abundance_multiplied = mass_array * abundance_array
    mass_iso_ele = sum(mass_abundance_multiplied) * ele_at_ratio
    return mass_iso_ele


def get_xy(isotopes, file_names, energy_min, energy_max, iso_abundance, sub_x, ele_at_ratio, _natural_mix, _unnatural_ratio_dict):
    # Transmission calculation of summed and separated contributions by each isotopes
    df = pd.DataFrame()
    df_raw = pd.DataFrame()
    y_i_iso_ele_dict = {}
    # thick_cm = thick_mm/10
    y_i_iso_ele_sum = 0.

    if _natural_mix == 'Y':
        iso_at_ratio = iso_abundance
    else:
        ratio_list = list(dict.values(_unnatural_ratio_dict))
        ratio_array = np.array(ratio_list)
        iso_at_ratio = ratio_array
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
        spline = interpolate.interp1d(x=df['E_eV'], y=df['Sig_b'], kind='linear')
        y_i = spline(x_energy)
        # y_i_sum = y_i_sum + y_i * iso_abundance[i] * ele_at_ratio
        y_i_iso_ele_dict[_isotope] = y_i * iso_at_ratio[i] * ele_at_ratio
        y_i_iso_ele_sum = y_i_iso_ele_sum + y_i * iso_at_ratio[i] * ele_at_ratio

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

    return x_energy, y_i_iso_ele_dict, y_i_iso_ele_sum, df_raw


def set_xy(_all, thick_mm, mixed_atoms_per_cm3, sig_dict, _x_axis):
    for _i in _all:
        _y_axis_i = _functions.sig2trans_quick(thick_mm, mixed_atoms_per_cm3, sig_dict[_i])
        plt.plot(_x_axis, _y_axis_i, label=_i)
    return


def plot_xy(_all, _energy_x_axis, _trans_y_axis, _plot_each_contribution, _plot_mixed,
            x_energy, y_trans_tot, thick_mm, mixed_atoms_per_cm3, sig_iso_ele_dict):
    if _energy_x_axis == 'Y':
        _x_axis = x_energy
        _x_words = 'Energy (eV)'
    else:
        _x_axis = _functions.ev2lamda(x_energy)
        _x_words = 'Wavelength (Å)'

    if _trans_y_axis == 'Y':
        _y_axis = y_trans_tot
        _y_words = 'Neutron transmission'
    else:
        _y_axis = 1 - y_trans_tot
        _y_words = 'Neutron attenuation'

    if _plot_each_contribution == 'Y':
        if _trans_y_axis == 'Y':
            for _i in _all:
                _y_axis_i = _functions.sig2trans_quick(thick_mm, mixed_atoms_per_cm3, sig_iso_ele_dict[_i])
                plt.plot(_x_axis, _y_axis_i, label=_i)
        else:
            for _i in _all:
                _y_axis_i = _functions.sig2trans_quick(thick_mm, mixed_atoms_per_cm3, sig_iso_ele_dict[_i])
                plt.plot(_x_axis, 1 - _y_axis_i, label=_i)

    if _plot_mixed == 'Y':
        plt.plot(_x_axis, _y_axis, label='Mixture')

    plt.ylim(-0.01, 1.01)
    plt.xlabel(_x_words)
    plt.ylabel(_y_words)
    plt.legend(loc='best')
    plt.show()

