import glob
import os
import numbers
import re
import numpy as np
import periodictable as pt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.constants import Avogadro


def is_element_in_database(element='', database='ENDF_VIII'):
    """will try to find the element in the folder (database) specified
    
    Parameters:
    ==========
    element: string. Name of the element. Not case sensitive
    database: string (default is 'ENDF_VIII'). Name of folder that has the list of elements
    
    Returns:
    =======
    bool: True if element was found in the database
          False if element could not be found
    """
    if element == '':
        return False
    
    list_entry_from_database = get_list_element_from_database(database=database)
    if element.lower() in list_entry_from_database:
        return True
    return False    
    
def get_list_element_from_database(database=''):
    """return a string array of all the element from the database
    
    Parameters:
    ==========
    database: string. Name of database
    
    Raises:
    ======
    ValueError if database can not be found
    
    """
    _file_path = os.path.abspath(os.path.dirname(__file__))
    _database_folder = os.path.join(_file_path, 'reference_data', database)

    if not os.path.exists(_database_folder):
        raise ValueError("Database {} does not exist!".format(database))
    
    _list_files = glob.glob(_database_folder + '/*.csv')
    _list_short_files = [os.path.basename(_file) for _file in _list_files]
    _list_element = set([_name.split('-')[0].lower() for _name in _list_short_files])
    return _list_element
    
def checking_stack(stack={}, database='ENDF_VIII'):
    """This method makes sure that all the elements from the various stacks are 
    in the database and that the thickness has the correct format (float)
    
    Parameters:
    ==========
    stack: dictionary that defines the various stacks
    database: string (default is 'ENDF_VIII') name of database
    
    Raises:
    ======
    ValueError if one of the element in one of the stack can not be found 
    ValueError if thickness is not a float
    ValueError if elements and stoichiometric ratio arrays do not have the same size
    
    Return:
    ======
    True: for testing purpose only
    """
    for _keys in stack:
        _elements = stack[_keys]['elements']
        for _element in _elements:
            if not is_element_in_database(element=_element):
                raise ValueError("Element {} can not be found in the database".format(_element))

        _thickness = stack[_keys]['thickness']['value']
        if not isinstance(_thickness, numbers.Number):
            raise ValueError("Thickness {} should be a number!".format(_thickness))
    
        _stoichiometric_ratio = stack[_keys]['stoichiometric_ratio']
        if len(_stoichiometric_ratio) != len(_elements):
            raise ValueError("stoichiometric Ratio and Elements should have the same size!")
    
    return True    

def formula_to_dictionary(formula='', thickness=np.NaN, density=np.NaN, database='ENDF_VIII'):
    """create dictionary based on formula given
    
    Parameters:
    ===========
    formula: string
       ex: 'AgCo2'
       ex: 'Ag'
    thickness: float (in mm) default is np.NaN
    density: float (in g/cm3) default is np.NaN
    database: string (default is ENDV_VIII). Database where to look for elements
    
    Raises:
    =======
    ValueError if one of the element is missing from the database
    
    Return:
    =======
    the dictionary of the elements passed
      ex: {'AgCo2': {'elements': ['Ag','Co'],
                     'stoichiometric_ratio': [1,2],
                     'thickness': {'value': thickness,
                                   'units': 'mm',
                                   },
                     'density': {'value': density,
                                 'units': 'g/cm3',
                                 },
                    }
    """
    _formula_parsed = re.findall(r'([A-Z][a-z]*)(\d*)', formula)

    _dictionary = {}
    _elements_array = []
    _atomic_ratio_array = []
    for _element in _formula_parsed:
        [_single_element, _atomic_ratio] = list(_element)
        if not is_element_in_database(element=_single_element, database=database):
            raise ValueError("element {} not found in database!".format(_single_element))

        if _atomic_ratio == '':
            _atomic_ratio = 1

        _atomic_ratio_array.append(int(_atomic_ratio))
        _elements_array.append(_single_element)

    return {formula: {'elements': _elements_array,
                      'stoichiometric_ratio': _atomic_ratio_array,
                      'thickness': {'value': thickness,
                                    'units': 'mm'},
                      'density': {'value': density,
                                  'units': 'g/cm3',
                                  }
                      },
            }

def get_isotope_dicts(element='', database='ENDF_VIII'):
    """return a dictionary with list of isotopes found in database and name of database files
    
    Parameters:
    ===========
    element: string. Name of the element
      ex: 'Ag'
    database: string (default is ENDF_VIII)
    
    Returns:
    ========
    dictionary with isotopes and files 
      ex: {'Ag': {'isotopes': ['107-Ag','109-Ag'],
                  'file_names': ['Ag-107.csv','Ag-109.csv']}}
    
    """
    _file_path = os.path.abspath(os.path.dirname(__file__))
    _database_folder = os.path.join(_file_path, 'reference_data', database)
    _element_search_path = os.path.join(_database_folder, element + '-*.csv')
    list_files = glob.glob(_element_search_path)
    isotope_dict = {'isotopes': {'list': [], 
                                 'file_names': [],
                                 'density': {'value': np.NaN,
                                             'units': 'g/cm3'},
                                 'mass': {'value': [],
                                          'units': 'g/mol',
                                 },
                                 'isotopic_ratio': [],},
                    'density': {'value': np.NaN,
                                'units': 'g/cm3'},
                    'molar_mass': {'value': np.NaN,
                                   'units': 'g/mol'},
                    } 

    isotope_dict_mirror = {}
    _isotopes_list = []
    _isotopes_list_files = []
    _isotopes_mass = []
    _isotopes_density = []
    _isotopes_atomic_ratio = []
    _density = np.NaN
    _molar_mass = np.NaN
    
    for file in list_files:

        # Obtain element, z number from the basename
        _basename = os.path.basename(file)
        [filename, file_extension] = os.path.splitext(_basename)
        [_name, _number] = filename.split('-')
        _symbol = _number + '-' + _name
        isotope = str(_symbol)

        _isotopes_list.append(isotope)
        _isotopes_list_files.append(_basename)
        _isotopes_mass.append(get_mass(isotope))
        _isotopes_atomic_ratio.append(get_abundance(isotope))
        _isotopes_density.append(get_density(isotope))
        _density = get_density(element)
        _molar_mass = get_mass(element)
                                        
    isotope_dict['isotopes']['list'] = _isotopes_list
    isotope_dict['isotopes']['file_names'] = _isotopes_list_files              
    isotope_dict['isotopes']['mass']['value'] = _isotopes_mass
    isotope_dict['isotopes']['isotopic_ratio'] = _isotopes_atomic_ratio
    isotope_dict['isotopes']['density']['value'] = _isotopes_density
    isotope_dict['density']['value'] = _density
    isotope_dict['molar_mass']['value'] = _molar_mass
                    
    return isotope_dict   

def get_abundance(element):
    """return the abundance [0.,1.] of the defined element
    
    Parameters:
    ===========
    element: (sting)
    
    Returns:
    ========
    the abundance of the elemenet (value between 0 and 1)
    """
    return pt.elements.isotope(element).abundance / 100.

def get_mass(element):
    """return the molar mass (SI units) of an given isotope, or element
    
    Parameters:
    ===========
    element: string. Element or isotopes to get the mass from
    
    Returns:
    ========
    the molar mass of the element/isotope in SI units
    """
    return pt.elements.isotope(element).mass

def get_density(element):
    """return the density (g.cm-3) of the element
    
    Paramters:
    ==========
    element: string. Name of element
    
    Returns:
    ========
    the density of the element in g.cm-3 units
    """
    return pt.elements.isotope(element).density

def get_compound_density(list_density=[], list_ratio=[]):
    """"""
    _ratio_density = zip(list_ratio, list_density)
    _density_compound = 0
    
    _sum_ratio = np.array(list_ratio).sum()
    
    for _ratio, _density in _ratio_density:
        _density_compound += (_ratio * _density) / _sum_ratio
    return _density_compound

def get_database_data(file_name=''):
    """return the energy (eV) and Sigma (barn) from the file_name
    
    Parameters:
    ===========
    file_name: string ('' by default) name of csv file
    
    Returns:
    ========
    pandas dataframe
    
    Raises:
    =======
    IOError if file does not exist
    """
    if not os.path.exists(file_name):
        raise IOError("File {} does not exist!".format(file_name))
    df = pd.read_csv(file_name, header=1)
    return df

def get_interpolated_data(df=[], E_min=np.NaN, E_max=np.NaN, E_step=np.NaN):
    """return the interpolated x and y axis for the given x range [E_min, E_max] with step defined
    
    Parameters:
    ===========
    df: data frame 
    E_min: left range of new interpolated data
    E_max: right range of new interpolated data
    E_step: step of energy to use in interpolated data
    
    Returns:
    ========
    x_axis and y_axis of interpolated data over specified range
    """
    nbr_point = (E_max - E_min) / E_step
    
    # remove data outside specified range [x_min, x_max]
    #df = df.drop(df[df.E_eV < E_min].index)
    #df = df.drop(df[df.E_eV > E_max].index)

    # reset index
    df = df.reset_index(drop=True)
    
    # energy x_axis
    #x_axis = np.linspace(df['E_eV'].min(), df['E_eV'].max(), nbr_point)
    x_axis = np.linspace(E_min, E_max, nbr_point)
    y_axis_function = interp1d(x=df['E_eV'], y=df['Sig_b'], kind='linear')
    
    y_axis = y_axis_function(x_axis)   
    
    return {'x_axis': x_axis, 'y_axis': y_axis}

def get_sigma(database_file_name='', E_min=np.NaN, E_max=np.NaN, E_step=np.NaN):
    """retrieve the Energy and sigma axis for the given isotope
    
    Paramters:
    ==========
    database_file_name: string
    E_min: left range of new interpolated data
    E_max: right range of new interpolated data
    E_step: step of energy to use in interpolated data
    
    Returns:
    ========
    {'energy': np.array(), 'sigma': np.array}
    """
    _df = get_database_data(file_name=database_file_name)
    _dict = get_interpolated_data(df=_df, E_min=E_min, E_max=E_max, 
                                 E_step=E_step)
    return {'energy_eV': _dict['x_axis'],
            'sigma_b': _dict['y_axis']}

def get_atoms_per_cm3_of_layer(compound_dict={}):
    """calculate the atoms per cm3 of the given compound (layer)
    
    Paramters:
    ==========
    compound_dict: {} 
    
    Returns:
    ========
    dictionary
    """
    atoms_per_cm3 = {}
    
    _list_of_elements = compound_dict['elements']
    _stoichiometric_list = compound_dict['stoichiometric_ratio']

    _element_stochio = zip(_list_of_elements, _stoichiometric_list)
    _molar_mass_sum = 0
    for _element, _stochio in _element_stochio:
        _molar_mass_sum += _stochio * compound_dict[_element]['molar_mass']['value']

    _element_stochio = zip(_list_of_elements, _stoichiometric_list)
    for _element, _stochio in _element_stochio:
        _step1 = (compound_dict['density']['value'] * _stochio) / _molar_mass_sum
        atoms_per_cm3[_element] = Avogadro * _step1
    
    return atoms_per_cm3

def calculate_transmission(thickness_cm=np.NaN, atoms_per_cm3=np.NaN, sigma_b=[]):
    """calculate the transmission signal using the formula
    
    transmission = exp( - thickness_cm * atoms_per_cm3 * 1e-24 * sigma_b)
    
    Parameters:
    ===========
    thickness: float (in cm)
    atoms_per_cm3: float (number of atoms per cm3 of element/isotope)
    sigma_b: np.array of sigma retrieved from database
    
    Returns:
    ========
    transmission array
    """
    transmission = np.exp( -thickness_cm * 1e-24 * sigma_b * atoms_per_cm3)
    return np.array(transmission)

def set_distance_units(value=np.NaN, from_units='mm', to_units='cm'):
    """convert distance into new units
    
    Parameters:
    ===========
    value: float. value to convert
    from_units: string. Must be 'mm', 'cm' or 'm'
    to_units: string. must be 'mm','cm' or 'm'
    
    Returns:
    ========
    converted value
    
    Raises:
    =======
    ValueError if from_units is not a valid unit (see above)
    ValueError if to_units is not a valud unit
    """
    if from_units == to_units:
        return value

    coeff = 1
    if from_units == 'cm':
        if to_units == 'mm':
            coeff = 10
        elif to_units == 'm':
            coeff = 0.01
        else:
            raise ValueError("to_units not supported ['cm','m','mm']!")
    elif from_units == 'mm':
        if to_units == 'cm':
            coeff = 0.1
        elif to_units == 'm':
            coeff = 0.001
        else:
            raise ValueError("to_units not supported ['cm','m','mm']!")
    elif from_units == 'm':
        if to_units == 'mm':
            coeff = 1000
        elif to_units == 'cm':
            coeff = 100
        else:
            raise ValueError("to_units not supported ['cm','m','mm']!")
    else:
        raise ValueError("to_units not supported ['cm','m','mm']!")

    return coeff * value

def energy_to_lambda(energy_ev=[]):
    """convert into lambda from the energy array

    Parameters:
    ===========
    energy: array (in eV)

    Returns:
    ========
    lambda: array (in Angstroms)
    """
    energy_mev = energy_ev * 1000
    lambda_array = np.sqrt(81.787 / energy_mev)
    return lambda_array

def energy_to_time(energy_ev=[], delay_us=np.NaN, source_to_detector_cm=np.NaN):
    # delay values is normal 2.99 us with NONE actual MCP delay settings
    """convert energy (eV) to time (us)

    Parameters:
    ===========
    energy: array (in eV)

    Returns:
    ========
    time: array in us (micro seconds)
    """
    energy_mev = energy_ev * 1000
    time_tot_us = np.sqrt(81.787 / energy_mev) * source_to_detector_cm / 0.3956
    time_record_us = time_tot_us - delay_us
    time_record_ns = time_record_us * 1000
    time_record_s = time_record_us / 1e6
    return time_record_us

def energy_to_image_number(energy_ev=[], delay_us=np.NaN, time_resolution_us=np.NaN, source_to_detector_cm=np.NaN):
    # delay values is normal 2.99 us with NONE actual MCP delay settings
    """convert energy (eV) to image numbers (#)

    Parameters:
    ===========
    energy: array (in eV)

    Returns:
    ========
    image numbers: array in (micro seconds)
    """
    energy_mev = energy_ev * 1000
    time_tot_us = np.sqrt(81.787 / energy_mev) * source_to_detector_cm / 0.3956
    time_record_us = (time_tot_us - delay_us)
    image_number = time_record_us / time_resolution_us
    return image_number


def time_to_energy(time_record_s, delay_us=np.NaN, source_to_detector_cm=np.NaN):
    """convert time (s) to energy (eV)
    Parameters:
    ===========
    time (in s)

    Returns:
    ========
    energy: (in eV)
    """
    time_tot_us = 1e6 * time_record_s + delay_us
    energy_mev = 81.787 / (0.3956 * time_tot_us / source_to_detector_cm) ** 2
    energy_ev = energy_mev / 1000
    return energy_ev



def ev2lamda(energy_ev):  # function to convert energy in eV to angstrom
    energy_miliev = energy_ev * 1000
    lamda = np.sqrt(81.787 / energy_miliev)
    return lamda

def time2lamda(time_tot_s, delay_us, source_to_detector_cm):  # function to convert time in us to angstrom
    time_tot_us = 1e6 * time_tot_s + delay_us
    lamda = 0.3956 * time_tot_us / source_to_detector_cm
    return lamda

def lamda2ev(lamda):  # function to convert angstrom to eV
    energy_miliev = 81.787 / (lamda ** 2)
    energy_ev = energy_miliev / 1000
    return energy_ev

def time2ev(time_tot_s, delay_ms=0.00299, source_to_detector_cm=1612.5):  # function to convert time in us to energy in eV
    delay_us = delay_ms * 1000
    time_tot_us = 1e6 * time_tot_s + delay_us
    energy_miliev = 81.787 / (0.3956 * time_tot_us / source_to_detector_cm) ** 2
    energy_ev = energy_miliev / 1000
    return energy_ev


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

def get_spectra2(_filename, delay_us, source_to_detector_cm, time_lamda_ev_axis='eV'):
    df_spectra = pd.read_csv(_filename, sep='\t', header=None)
    time_array = (np.array(df_spectra[0]))
    counts_array = (np.array(df_spectra[1]))
    if time_lamda_ev_axis == 'lamda':
        lamda_array = time2lamda(time_array, delay_us, source_to_detector_cm)
        return lamda_array, counts_array
    if time_lamda_ev_axis == 'eV':
        ev_array = time2ev(time_array, delay_us, source_to_detector_cm)
        return ev_array, counts_array
    if time_lamda_ev_axis == 'time':
        return time_array, counts_array


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

'''
Energy, wavelength and time conversions
'''

def ev2lamda(energy_ev):  # function to convert energy in eV to angstrom
    energy_miliev = energy_ev * 1000
    lamda = np.sqrt(81.787 / energy_miliev)
    return lamda

def time2lamda(time_tot_s, delay_us, source_to_detector_cm):  # function to convert time in us to angstrom
    time_tot_us = 1e6 * time_tot_s + delay_us
    lamda = 0.3956 * time_tot_us / source_to_detector_cm
    return lamda

def lamda2ev(lamda):  # function to convert angstrom to eV
    energy_miliev = 81.787 / (lamda ** 2)
    energy_ev = energy_miliev / 1000
    return energy_ev

def time2ev(time_tot_s, delay_us, source_to_detector_cm):  # function to convert time in us to energy in eV
    time_tot_us = 1e6 * time_tot_s + delay_us
    energy_miliev = 81.787 / (0.3956 * time_tot_us / source_to_detector_cm) ** 2
    energy_ev = energy_miliev / 1000
    return energy_ev