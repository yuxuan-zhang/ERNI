import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import _functions
import glob
import _fit_funtions
from lmfit import minimize, Parameters

# Parameters
# Parameters
source_to_detector_cm = 1610.9  # cm
delay_ms = 4.5 - 16.6  # ms
delay_us = delay_ms * 1000
_slice = 220
time_lamda_ev_axis = 'eV'
_name = 'foil6'
path = 'data/' + _name + '.csv'
spectra_path = 'data/spectra.txt'
file_names = glob.glob(path)
# _fit_funtions.get_multi_data(file_name_signature, time_lamda_ev_axis, delay_us, source_to_detector_cm, _slice)


x_axis_array = _functions.get_spectra_slice(spectra_path, time_lamda_ev_axis, delay_us,
                                                source_to_detector_cm, _slice)
print(x_axis_array)
y_axis_array = _functions.get_normalized_data_slice('data/'+_name+'.csv', _slice)
print(y_axis_array)

paramas = Parameters()
paramas.add('source_to_detector_cm', value=1610.9)
paramas.add('delay_ms', value=4.5-16.6)


def get_ele_database(_database, _element):
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



# def get_distance_delay(paramas, ):
