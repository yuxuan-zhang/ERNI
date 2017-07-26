import _plot_functions
import _functions
import numpy as np
import pandas as pd

'''Describe your sample: '''
# Input sample name or names as str, case sensitive
_input_formula = 'UOGd'  # input('Please input the chemicals? ')
_input_thick_mm = 0.26  # float(input('Please input the thickness or majority thickness of stacked foils in mm : '))
_input_thick_cm = _input_thick_mm/10
_database = 'ENDF_VIII'
energy_max = 300  # max incident energy in eV
energy_min = 0  # min incident energy in eV
energy_sub = 100  # steps used to interpolate database
sub_x = energy_sub * (energy_max - energy_min)  # steps used to interpolate database
compound_boo = 'Y'  # Compound or single/multi elements foil/stacked foils: Y/N?

'''Input for dict modification in certain cases: '''
# Thickness input:
special_thick_boo = 'N'
special_thick_element_str = str
special_thick_mm_list = []
special_thick_cm_list = np.array(special_thick_mm_list)/10
# Enriched isotope ratio input:
enrichment_boo = 'Y'  # Isotopic enriched or depleted: Y/N?
enriched_element_str = 'U'
input_ratio_dict = {'U': [0., 0., .15, .85]}
                    # 'O': [1., 0., 0.]}  #{'233-U': 0., '234-U': 0., '235-U': 0.15, '238-U': 0.85}}
# Special density input:
special_density_boo = 'Y'
special_density_element_str = str
special_density_gcm3_list = []

'''How you want the data to be plotted?'''
_plot_or_not = 'Y'
_energy_x_axis = 'Y'  # 1 means plot x-axis as energy in eV
_trans_y_axis = 'N'  # 1 means plot y-axis as transmission
_plot_each_ele_contribution = 'Y'  # 1 means plot each element's contribution
_plot_each_iso_contribution = 'N'  # 1 means plot each isotope's contribution
_plot_mixed = 'N'  # 1 means plot mixed resonance
'''Export to clipboard for Excel or DataGraph?'''
_export_to_clipboard_boo = 'N'


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
formula_dict = _functions.input2formula(_input_formula)
elements = _functions.dict_key_list(formula_dict)
ratios = _functions.dict_value_list(formula_dict)
sum_ratios = sum(ratios)
isotope_dict = _functions.get_isotope_dicts(_database, elements)

# DICT 1: Thickness dict with option for modification
thick_cm_dict = _functions.repeat_value_dict(elements, _input_thick_cm)

# DICT 2: Isotopic mass dict
iso_mass_dicts = _functions.get_iso_mass_dicts_quick(elements, isotope_dict)

# Dict 3: Molar mass dict
molar_mass_dict = _functions.get_molar_mass_dict(elements)

# DICT 4: Isotope at.% dict with option for modification
iso_ratio_dicts = _functions.get_iso_ratio_dicts_quick(elements, isotope_dict)

# DICT 5: Density dict
density_gcm3_dict = _functions.get_density_dict(elements)


# Update DICT 2: isotopic ratio changes lead to |Density| & |Molar mass| changes
if compound_boo == 'N':
    if special_thick_boo == 'Y':
        thick_cm_dict = _plot_functions.modify_thick_cm_dict_by_input(thick_cm_dict, special_thick_element_str, special_thick_cm_list)

# Update DICT 3 & 4 & 5: isotopic ratio changes lead to |Density| & |Molar mass| changes
if enrichment_boo == 'Y':
    # Update isotope at.% ratio dict
    iso_ratio_dicts, enriched_element = _plot_functions.modify_iso_ratio_dicts(elements, isotope_dict, enriched_element_str, input_ratio_dict)
    # Update molar mass dict
    molar_mass_dict = _plot_functions.modify_molar_mass_dict_by_enrichment(molar_mass_dict, enriched_element, isotope_dict, iso_ratio_dicts, iso_mass_dicts)
    # Update density dict
    density_gcm3_dict = _plot_functions.modify_density_dict_by_enrichment(density_gcm3_dict, enriched_element, isotope_dict, iso_ratio_dicts)

# Update DICT 5: Density dict, if special case encountered
if compound_boo == 'N':
    if special_density_boo == 'Y':
        # Stacked foils and would like to modify density for specific element
        density_gcm3_dict = _plot_functions.modify_density_dict_by_input(density_gcm3_dict, special_density_element_str, special_density_gcm3_list)
else:
    if special_density_boo == 'Y':
        # Not isolated elements or mixture or compound need density input currently
        input_tot_density = 0.7875

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
# atoms_per_cm3_dict = {}

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
        = _plot_functions.get_xy_from_database(iso_ratio_dicts[el],
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
if compound_boo == 'N':
    # Stacked foils or single foil
    mixed_l_n_avo = _plot_functions.l_x_n_multi_ele_stack(elements,
                                                          thick_cm_dict,
                                                          density_gcm3_dict,
                                                          molar_mass_dict)
else:
    # For compound
    thick_cm_list = list(dict.values(thick_cm_dict))
    thick_cm = thick_cm_list[0]
    compound_density = input_tot_density
    mixed_l_n_avo = _plot_functions.l_x_n_compound(elements,
                                                   thick_cm,
                                                   compound_density,
                                                   molar_mass_dict,
                                                   formula_dict,
                                                   sum_ratios)

# Get the tot transmission for all
# yi_values_l = list(dict.values(sigma_iso_ele_sum_l_eledict))
# yi_values_l_sum = sum(yi_values_l)
# # sum of (sigma * ele_ratio * iso_ratio * l)
yi_values = list(dict.values(sigma_iso_ele_sum_eledict))
yi_values_sum = sum(yi_values)
# sum of (sigma * ele_ratio * iso_ratio)
# print(yi_values)
trans_sum = _functions.sig_l_2trans_quick(mixed_l_n_avo, yi_values_sum)
y_trans_tot = trans_sum


# Create the trans or absorb dict of ele for plotting if needed
y_ele_dict = {}
if _plot_each_ele_contribution == 'Y':
    for _ele in elements:
        if _trans_y_axis == 'Y':
            y_ele_dict[_ele] = _functions.sig_l_2trans_quick(mixed_l_n_avo, sigma_iso_ele_sum_eledict[_ele])
        else:
            y_ele_dict[_ele] = 1 - _functions.sig_l_2trans_quick(mixed_l_n_avo, sigma_iso_ele_sum_eledict[_ele])

# Create the trans or absorb dict : y_iso_dicts of isotopes for plotting if needed
y_iso_dicts = {}
if _plot_each_iso_contribution == 'Y':
    for _ele in elements:
        y_iso_dict = {}
        for _iso in isotope_dict[_ele]:
            if _trans_y_axis == 'Y':
                y_iso_dict[_iso] = _functions.sig_l_2trans_quick(mixed_l_n_avo,
                                                                 sigma_iso_ele_eleisodict[_ele][_iso])
            else:
                y_iso_dict[_iso] = 1 - _functions.sig_l_2trans_quick(mixed_l_n_avo,
                                                                     sigma_iso_ele_eleisodict[_ele][_iso])
        y_iso_dicts[_ele] = y_iso_dict

# Plot the theoretical neutron resonance
if _plot_or_not == 'Y':
    _plot_functions.plot_database(_energy_x_axis,
                                  _trans_y_axis,
                                  _plot_mixed,
                                  _plot_each_ele_contribution,
                                  _plot_each_iso_contribution,
                                  elements,
                                  isotope_dict,
                                  x_energy,
                                  y_trans_tot,
                                  y_ele_dict,
                                  y_iso_dicts,
                                  _input_formula)

# Export to clipboard for density and thickness manipulations with Excel or DataGraph
if _export_to_clipboard_boo == 'Y':
    _name = _input_formula
    df_yi_tot = pd.DataFrame(data=x_energy, index=None)
    df_yi_tot.rename(columns={0: 'eV' + _name}, inplace=True)
    df_yi_tot['lamda-' + _name] = _functions.ev2lamda(x_energy)
    df_yi_tot['sample_density-' + _name] = compound_density
    df_yi_tot['avo_divided-' + _name] = mixed_l_n_avo
    df_yi_tot['sigma-' + _name] = yi_values_sum

    for ele in elements:
        _ele_str = str(ele)
        df_yi_tot['sigma-' + _ele_str] = sigma_iso_ele_sum_eledict[ele]
        df_test = pd.DataFrame(sigma_iso_ele_eleisodict[ele])
        df_yi_tot = pd.concat([df_yi_tot, df_test], axis=1)
    print(df_yi_tot.head())
    # Export to clipboard
    df_yi_tot.to_clipboard(excel=True)
