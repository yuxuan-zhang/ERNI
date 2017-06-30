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
def time2ev(time_tot):  # function to convert time in us to energy in eV
    energy = 81.787/(1000 * (0.3956 * time_tot/source_to_detector) ** 2)
    return energy
def atoms_per_cm3(density, mass):
    n_atoms = density * pt.constants.avogadro_number/mass
    print('Number of atoms per unit volume (#/cm^3): ', n_atoms)
    return n_atoms
def sig2trans(_thick_cm, _atoms_per_cm3, _ele_atomic_ratio, _sigma_b, _iso_atomic_ratio):
    neutron_transmission = np.exp(-1 * _thick_cm * _atoms_per_cm3 *
                                  _ele_atomic_ratio * _sigma_b * 1e-24 * _iso_atomic_ratio)
    return neutron_transmission
def sig2trans_quick(_thick_cm, _atoms_per_cm3, _sigma_portion_sum):
    neutron_transmission = np.exp(-1 * _thick_cm * _atoms_per_cm3 * 1e-24 * _sigma_portion_sum)
    return neutron_transmission

# Parameters
thick_mm = 0.3  # mm
thick_cm = thick_mm/10  # Thickness in cm
density_para = 8.65  # g/cm3  https://en.wikipedia.org/wiki/Cadmium
ele_at_ratio = 1  # for single element, will be implanted for multiple elements compound
element = 'Cd'
element2 = '16-O'
mass_abund_other = 0  # pt.elements.isotope(element2).mass * (pt.elements.isotope(element2).abundance/100)
_database = 'ENDF_VIII'
energy_max = 3000  # max incident energy in eV
energy_min = 1  # min incident energy in eV
sub_x = 100 * (energy_max - energy_min)  # subdivided new x-axis

### From .csv file basename obtain isotope name, isotope number, z number as dict.
### Convert dict to list which is callable in following steps.

main_dir = os.path.dirname(os.path.abspath(__file__))
path = main_dir + '/data_web/' + _database + '/' + element + '*.csv'
file_names = glob.glob(path)
abundance_dict = {}
density_dict = {}
mass_dict = {}
z_number = {}
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
    # abundance = pt.elements.Cd[z].abundance
    # print(type(abundance))
    # abundance = abundance.append(abundance)
isotopes = list(dict.keys(abundance_dict))
iso_abundance = list(dict.values(abundance_dict))
iso_density = list(dict.values(density_dict))
iso_mass = list(dict.values(mass_dict))
# print(isotopes)
# print(iso_abundance)
# print(iso_density)
# print(iso_mass)
# print(abundance_dict)
# print(density_dict)
# print(mass_dict)


### Calculate the number of atoms per unit volume

abundance_array = np.array(iso_abundance)
mass_array = np.array(iso_mass)
mass_abund = mass_array * abundance_array
sum_density = sum(mass_abund) * ele_at_ratio + mass_abund_other * (1 - ele_at_ratio)
print(mass_abund)
print(sum_density)
mixed_atoms_per_cm3 = density_para * pt.constants.avogadro_number/sum_density
print('Number of atoms per unit volume (#/cm^3): {}'.format(mixed_atoms_per_cm3))

# Transmission calculation of summed and separated contributions by each isotopes
df = pd.DataFrame()
df_raw = pd.DataFrame()
y_i_sum = 0.
y_trans_tot = 1
for i, _isotope in enumerate(isotopes):
    # Read database .csv file
    df = pd.read_csv(file_names[i], header=1) #names=['energy_inc', 'xs_tot'], float_precision='round_trip')
    # Drop rows beyond range
    df = df.drop(df[df.E_eV > energy_max].index)  # drop rows beyond range
    df = df.drop(df[df.E_eV < energy_min].index)  # drop rows beyond range
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

    # ## For getting transmission contribution of each isotope, use the following
    # y_trans_i = sig2trans(thick_cm, mixed_atoms_per_cm3, ele_at_ratio, y_i, iso_abundance[i])
    # y_trans_tot = y_trans_tot * y_trans_i
    # print(y_trans_i[:5])
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

x_lamda = ev2lamda(x_energy)
y_trans_tot = sig2trans_quick(thick_cm, mixed_atoms_per_cm3, y_i_sum)
y_absorb_tot = 1 - y_trans_tot
# plt.plot(x_energy, y_i_sum, label='sig_sum')
# plt.legend(loc='best')
# plt.show()

# plt.plot(x_energy, y_absorb_tot, label='Cd natural mixture')
plt.plot(x_energy, y_trans_tot, label='UO2 natural mixture')
plt.ylim(0, 1.01)
plt.xlabel('Energy (eV)')
plt.ylabel('Neutron transmission')
plt.legend(loc='best')
plt.show()
