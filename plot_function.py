import _functions
import numpy as np
import pandas as pd
import periodictable as pt
import matplotlib.pyplot as plt
from matplotlib import style
from periodictable import constants
import glob
import os
from scipy import interpolate
style.use('ggplot')


def _get_plot(thick_mm, element, density_sample, ele_at_ratio, mass_abund_other, _database,
              energy_max, energy_min, energy_sub, _type_x_axis, _type_y_axis, _plot_each_contribution, _plot_mixed):
    # # Parameters
    # thick_mm = 0.26  # mm
    # element = 'U'
    # density_sample = 2  # 0.7875  # pt.elements.isotope(element).density  # g/cm3  https://en.wikipedia.org/wiki/Cadmium
    # ele_at_ratio = 0.25  # for single element, will be implanted for multiple elements compound
    # # element2 = '16-O'
    # mass_abund_other = 0  # pt.elements.isotope(element2).mass * (pt.elements.isotope(element2).abundance/100)
    # _database = 'ENDF_VIII'
    # energy_max = 300  # max incident energy in eV
    # energy_min = 1  # min incident energy in eV
    # energy_sub = 100
    # _type_x_axis = 'energy'  # 1 means plot x-axis as energy in eV
    # _type_y_axis = 'absorb'  # 1 means plot y-axis as transmission
    # _plot_each_contribution = 'Y'  # 1 means plot each isotope contribution
    # _plot_mixed = 'N'  # 1 means plot mixed resonance

    # _lamda_axis = 1 - _type_x_axis
    # _absorb_axis = 1 - _type_y_axis
    thick_cm = thick_mm / 10  # Thickness in cm
    sub_x = energy_sub * (energy_max - energy_min)  # subdivided new x-axis

    ### From .csv file basename obtain isotope name, isotope number, z number as dict.
    ### Convert dict to list which is callable in following steps.

    main_dir = os.path.dirname(os.path.abspath(__file__))
    path = main_dir + '/data_web/' + _database + '/' + element + '*.csv'
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
        _number_csv = _basename.split('-')
        _number = _number_csv[1].split('.')
        # z_number[_i] = z = int(_number[0])  # int is required in order to be used to call mass/density in periodictable
        abundance_dict[str(_symbol)] = pt.elements.isotope(_symbol).abundance/100
        density_dict[str(_symbol)] = pt.elements.isotope(_symbol).density
        mass_dict[str(_symbol)] = pt.elements.isotope(_symbol).mass
    isotopes = list(dict.keys(abundance_dict))  # List of isotopes such as '238-U', ''235-U
    iso_abundance = list(dict.values(abundance_dict))  # List of isotopic abundance
    iso_density = list(dict.values(density_dict))  # List of isotopic density
    iso_mass = list(dict.values(mass_dict))  # List of isotopic molar mass
    iso_abundance = [0, 0, 0.15, 0.85]

    ### Calculate the number of atoms per unit volume

    abundance_array = np.array(iso_abundance)
    mass_array = np.array(iso_mass)
    mass_abund_multiplied = mass_array * abundance_array
    sum_density = sum(mass_abund_multiplied) * ele_at_ratio + mass_abund_other * (1 - ele_at_ratio)
    mixed_atoms_per_cm3 = density_sample * pt.constants.avogadro_number/sum_density

    print('Number of atoms per unit volume (#/cm^3): {}'.format(mixed_atoms_per_cm3))

    # Transmission calculation of summed and separated contributions by each isotopes
    df = pd.DataFrame()
    df_raw = pd.DataFrame()
    trans_dict = {}
    absorb_dict = {}
    y_i_sum = 0.
    y_trans_tot = 1
    for i, _isotope in enumerate(isotopes):
        # Read database .csv file
        df = pd.read_csv(file_names[i], header=1) #names=['energy_inc', 'xs_tot'], float_precision='round_trip')
        # Drop rows beyond range
        df = df.drop(df[df.E_eV > energy_max].index)  # drop rows beyond range
        df = df.drop(df[df.E_eV < energy_min].index)  # drop rows beyond range
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
        y_i = spline(x_energy)  # 100 at.% for each isotope
        y_i_sum = y_i_sum + y_i * iso_abundance[i] * ele_at_ratio  # ele_at_ratio = 1 in this case

        ## For getting transmission contribution of each isotope, use the following
        if _plot_each_contribution == 'Y':
            y_trans_i = _functions.sig2trans(thick_cm, mixed_atoms_per_cm3, ele_at_ratio, y_i, iso_abundance[i])
            y_trans_tot = y_trans_tot * y_trans_i
            trans_dict[_isotope] = y_trans_i
            absorb_dict[_isotope] = 1 - y_trans_i
        # # plt.plot(df['E_eV'], df['Sig_b'], label=_isotope)
        # # plt.plot(x_energy, y_i, 'b+', label='Linear')
        # plt.plot(x_energy, y_trans_i, label=_isotope)
        # # plt.ylim(0.99, 1)
        # plt.legend(loc='best')
        # plt.show()

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

    print(trans_dict)
    print(absorb_dict)
    print(df_raw.head())
    x_lamda = _functions.ev2lamda(x_energy)
    y_trans_tot = _functions.sig2trans_quick(thick_cm, mixed_atoms_per_cm3, y_i_sum)
    y_absorb_tot = 1 - y_trans_tot

    # plt.plot(x_energy, y_i_sum, label='sig_sum')
    # plt.legend(loc='best')
    # plt.show()

    if _type_x_axis == 'energy':
        _x_axis = x_energy
        _x_words = 'Energy (eV)'
    else:
        _x_axis = x_lamda
        _x_words = 'Wavelength (Å)'

    if _type_y_axis == 'trans':
        _y_axis = y_trans_tot
        _y_words = 'Neutron transmission'
    else:
        _y_axis = y_absorb_tot
        _y_words = 'Neutron attenuation'

    if _plot_each_contribution == 'Y':
        if _type_y_axis == 'trans':
            for _each in isotopes:
                plt.plot(_x_axis, trans_dict[_each], label=_each)
        else:
            for _each in isotopes:
                plt.plot(_x_axis, absorb_dict[_each], label=_each)
    print(isotopes)

    if _plot_mixed == 'Y':
        plt.plot(_x_axis, _y_axis, label=element+' natural mixture')

    plt.ylim(-0.01, 1.01)
    plt.xlabel(_x_words)
    plt.ylabel(_y_words)
    plt.legend(loc='best')
    plt.show()

