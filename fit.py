import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import _functions
import glob
import _fit_funtions
from lmfit import minimize, Parameters
import os
import periodictable as pt
from periodictable import constants
import peakutils as pku

# Parameters
source_to_detector_cm = 1610.9  # cm
delay_ms = 4.5 - 16.6127  # ms
delay_us = delay_ms * 1000
_slice = 220
_database = 'ENDF_VIII'
_input_element = 'Co'
energy_min = 0
energy_max = 400
energy_sub = 100
sub_x = energy_sub * (energy_max - energy_min)  # subdivided new x-axis
_thick_mm = 0.025
_thick_cm = _thick_mm/10
time_lamda_ev_axis = 'eV'
_name = 'foil6'
data_path = 'data/' + _name + '.csv'
spectra_path = 'data/spectra.txt'
# _fit_funtions.get_multi_data(file_name_signature, time_lamda_ev_axis, delay_us, source_to_detector_cm, _slice)


x_data_array = _functions.get_spectra_slice(spectra_path, time_lamda_ev_axis, delay_us,
                                                source_to_detector_cm, _slice)
print(x_data_array)
y_data_array = 1 - _functions.get_normalized_data_slice('data/'+_name+'.csv', _slice)/4.2
print(y_data_array)
plt.plot(x_data_array, y_data_array, 'r.', label=_name)

indexes = pku.indexes(y_data_array, thres=0.6, min_dist=50)

print(indexes)
peaks_x = pku.interpolate(x_data_array, y_data_array, ind=indexes)
print(peaks_x)
# plt.figure(figsize=(10, 6))
plt.plot(x_data_array[indexes], y_data_array[indexes], 'bx', label='peak')
plt.title('First estimation')


# paramas = Parameters()
# paramas.add('source_to_detector_cm', value=1610.9)
# paramas.add('delay_ms', value=4.5-16.6)
# paramas.add('density_gcm3', value=pt.elements.isotope(_input_element).density)

formula_dict = _functions.input2formula(_input_element)
elements = _functions.dict_key_list(formula_dict)
ratios = _functions.dict_value_list(formula_dict)

# To check whether the input are foils stacked
foils_stacked = ratios.count(ratios[0]) == len(ratios)
for _each_ in elements:
    if foils_stacked is True:
        ele_at_ratio = 1
### Get the data from database specified
    isotopes, abundance_dict, iso_abundance, iso_mass, file_names = _fit_funtions.get_pre_data_to_fit(_database, _input_element)
# Get (mass * iso_ratio * ele_ratio). In this case, ele_at_ratio = 1
    mass_iso_ele = _fit_funtions.get_mass_iso_ele_to_fit(iso_abundance, iso_mass, ele_at_ratio)
# Get density and then N (number of atoms per unit volume cm3)
    density_gcm3 = pt.elements.isotope(_input_element).density
    mixed_atoms_per_cm3 = density_gcm3 * pt.constants.avogadro_number/mass_iso_ele
# Get transmission and/or attenuation
    x_energy, sigma_iso_ele_isodict, sigma_iso_ele_sum, df_raw = \
        _fit_funtions.get_xy_to_fit(isotopes,
                                    file_names,
                                    energy_min,
                                    energy_max,
                                    iso_abundance,
                                    sub_x,
                                    ele_at_ratio)
aasdad = 1


def get_residual():


    model = _functions.sig2trans_quick(_thick_cm, mixed_atoms_per_cm3, sigma_iso_ele_sum)
    resudual = data - model
    return

trans_sum = _functions.sig2trans_quick(_thick_cm, mixed_atoms_per_cm3, sigma_iso_ele_sum)
attenu_sum = 1 - trans_sum

plt.plot(x_energy, attenu_sum, 'b-', label=_input_element+'-'+_database)

plt.xlim(0, 300)
# plt.ylim(-0.01, 1.01)
# plt.xlabel(_x_words)
# plt.ylabel(_y_words)
plt.legend(loc='best')
plt.show()
# Transmission calculation of summed and separated contributions by each isotopes


# def get_distance_delay(paramas, ):
