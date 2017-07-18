import _plot_functions
import _functions
import major_plot_function as mpf
import periodictable as pt
from periodictable import constants
import numpy as np


# Parameters
# _input_formula = elements_str
_input_formula = 'Co'  # input('Please input the chemicals? ')
_input_thick_mm = 0.025  # float(input('Please input the thickness or majority thickness of stacked foils in mm : '))
_input_thick_cm = _input_thick_mm/10
_database = 'ENDF_VIII'
energy_max = 300  # max incident energy in eV
energy_min = 0  # min incident energy in eV
energy_sub = 100
sub_x = energy_sub * (energy_max - energy_min)  # subdivided new x-axis
_plot_or_not = 'Y'
_energy_x_axis = 'Y'  # 1 means plot x-axis as energy in eV
_trans_y_axis = 'N'  # 1 means plot y-axis as transmission
_plot_each_ele_contribution = 'Y'  # 1 means plot each element's contribution
_plot_each_iso_contribution = 'N'  # 1 means plot each isotope's contribution
_plot_mixed = 'N'  # 1 means plot mixed resonance
stacked_foil_boo = 'Y'  # Stacked foils or not: Y/N


# Function below is called to parse input formula and return elements and ratios
formula_dict = _functions.input2formula(_input_formula)
elements = _functions.dict_key_list(formula_dict)
ratios = _functions.dict_value_list(formula_dict)
isotope_dict = _functions.get_isotope_dicts(_database, elements)

# DICT 1: Thickness dict with option for modification
special_thick_boo = 'N'
if special_thick_boo == 'Y':
    special_thick_element_str = str
    special_thick_cm_list = []
    thick_cm_dict = mpf.modify_thick_cm_dict(elements, _input_thick_cm, special_thick_element_str, special_thick_cm_list)
else:
    thick_cm_dict = _functions.get_thick_dict(elements, _input_thick_cm)

# DICT 2: Isotopic mass dict
iso_mass_dicts = _functions.get_iso_mass_dicts_quick(elements, isotope_dict)

# DICT 3: Isotope atomic ratio dict with option for modification
# DICT 4: Density dict is also constructed since the modification of isotopic ratio changes solid density
enrichment_boo = 'N'
if enrichment_boo == 'Y':
    enriched_element_str = str
    input_ratio_dict = {}
    iso_ratio_dicts = mpf.modify_iso_ratio_dicts(elements, isotope_dict, enriched_element_str, input_ratio_dict)
else:
    iso_ratio_dicts = _functions.get_iso_ratio_dicts_quick(elements, isotope_dict)
    density_gcm3_dict = _functions.get_density_dict(elements)

# DICT 4: Update density dict if special case involved
special_density_boo = 'N'
modify_individual_boo = 'N'
modify_whole_boo = 'N'
if special_density_boo == 'Y':
    if modify_individual_boo == 'Y':
        # Stacked foils and would like to modify density for specific element
        special_density_element_str = str
        special_density_gcm3_list = []
        density_gcm3_dict = mpf.modify_density_dict(density_gcm3_dict, special_density_boo, special_density_element_str, special_density_gcm3_list)
    if modify_whole_boo == 'Y':
        sample_density = input('Total sample density in g/cm3: ')

print('Thickness (cm): ', thick_cm_dict)
print('Density (g/cm^3): ', density_gcm3_dict)
print('Isotopic at.%', iso_ratio_dicts)
# To check whether the input are foils stacked
# foils_stacked = ratios.count(ratios[0] == len(ratios))


mass_iso_ele_dict = {}  # For number of atoms per cm3 calculation
sigma_iso_ele_eleisodict = {}  # For transmission calculation at isotope level
sigma_iso_ele_sum_eledict = {}  # For transmission calculation at element level
sigma_iso_ele_sum_l_eledict = {}
sigma_iso_ele_l_eleisodict = {}
df_raw_dict = {}  # Raw sigma data for elements and isotopes

if modify_whole_boo == 'N':
    sample_density = 0.

for el in elements:
    # isotopes_list = list(dict.keys(iso_ratio_dicts[el]))
    iso_ratio_list = list(dict.values(iso_ratio_dicts[el]))
    iso_ratio_array = np.array(iso_ratio_list)
    iso_mass_list = list(dict.values(iso_mass_dicts[el]))
    iso_mass_array = np.array(iso_ratio_list)
    # if foils_stacked == 'Y':
    #     ele_at_ratio = 1
    # else:
    #     ele_at_ratio = formula_dict[el] / sum(ratios)
    ele_at_ratio = formula_dict[el] / sum(ratios)
    # A part for getting atoms_per_cm3, this part is irrelevant to fitting parameters, and will be exported for fitting
    if stacked_foil_boo == 'Y':
        mass_iso_ele_dict[el] = sum(iso_ratio_array * iso_mass_array) * 1
    else:
        mass_iso_ele_dict[el] = sum(iso_ratio_array * iso_mass_array) * ele_at_ratio

    # Total density calculation of mixture mixed by ele_at_ratio
    if modify_whole_boo == 'N':
        sample_density = sample_density + density_gcm3_dict[el] * ele_at_ratio

    file_names = _functions.get_file_path(_database, el)
    x_energy, sigma_iso_ele_isodict, sigma_iso_ele_l_isodict, sigma_iso_ele_sum, df_raw_dict[el] \
        = _plot_functions.get_xy_new(iso_ratio_dicts[el],
                                     thick_cm_dict[el],
                                     file_names,
                                     energy_min,
                                     energy_max,
                                     iso_ratio_list,
                                     sub_x,
                                     ele_at_ratio)
    # Two level dict of isotopic array of (L * sigma * iso_ratio * ele_ratio)
    sigma_iso_ele_l_eleisodict[el] = sigma_iso_ele_l_isodict
    # One level dict of elemental array of (L * sigma * iso_ratio * ele_ratio)
    sigma_iso_ele_sum_l_eledict[el] = sigma_iso_ele_sum * thick_cm_dict[el]

    # Two level dict of isotopic array of (sigma * iso_ratio * ele_ratio)
    sigma_iso_ele_eleisodict[el] = sigma_iso_ele_isodict
    # One level dict of elemental array of (sigma * iso_ratio * ele_ratio)
    sigma_iso_ele_sum_eledict[el] = sigma_iso_ele_sum


print(sample_density)
mass_iso_ele_list = list(dict.values(mass_iso_ele_dict))
mass_iso_ele_sum = sum(np.array(mass_iso_ele_list))
avo_divided = pt.constants.avogadro_number/mass_iso_ele_sum
mixed_atoms_per_cm3 = sample_density * avo_divided
# Use function: mixed_atoms_per_cm3 = _functions.atoms_per_cm3(density=sample_density, mass=mass_iso_ele_sum)

# sum of (sigma * ele_ratio * iso_ratio * l)
yi_values_l = list(dict.values(sigma_iso_ele_sum_l_eledict))
yi_values_l_sum = sum(yi_values_l)
# sum of (sigma * ele_ratio * iso_ratio)
yi_values = list(dict.values(sigma_iso_ele_sum_eledict))
yi_values_sum = sum(yi_values)


trans_sum = _functions.sigl2trans_quick(mixed_atoms_per_cm3, yi_values_l_sum)
y_trans_tot = trans_sum

### Create the trans or absorb dict of ele for plotting if needed
if _plot_each_ele_contribution == 'Y':
    y_ele_dict = {}
    for _ele in elements:
        if _trans_y_axis == 'Y':
            y_ele_dict[_ele] = _functions.sigl2trans_quick(mixed_atoms_per_cm3, sigma_iso_ele_sum_l_eledict[_ele])
        else:
            y_ele_dict[_ele] = 1 - _functions.sigl2trans_quick(mixed_atoms_per_cm3, sigma_iso_ele_sum_l_eledict[_ele])

y_iso_dicts = {}
### Create the trans or absorb dict : y_iso_dicts of isotopes for plotting if needed
if _plot_each_iso_contribution == 'Y':

    y_iso_dict = {}
    for _ele in elements:
        for _iso in isotope_dict[_ele]:
            if _trans_y_axis == 'Y':
                y_iso_dict[_iso] = _functions.sigl2trans_quick(mixed_atoms_per_cm3,
                                                               sigma_iso_ele_l_eleisodict[_ele][_iso])
            else:
                y_iso_dict[_iso] = 1 - _functions.sigl2trans_quick(mixed_atoms_per_cm3,
                                                                   sigma_iso_ele_l_eleisodict[_ele][_iso])
        y_iso_dicts[_ele] = y_iso_dict
        y_iso_dict = {}  # Clear for following set of isotopes
    # print(y_iso_dicts)


# ### Export to clipboard for density and thickness manipulations with Excel or DataGraph
# _name = _input_formula
# df_yi_tot = pd.DataFrame(data=x_energy, index=None)
# df_yi_tot.rename(columns={0: 'eV'+_name}, inplace=True)
# df_yi_tot['lamda-'+_name] = _functions.ev2lamda(x_energy)
# df_yi_tot['sample_density-'+_name] = sample_density
# df_yi_tot['avo_divided-'+_name] = avo_divided
# df_yi_tot['sigma-'+_name] = yi_values_sum
#
#
# # print(y_i_iso_ele_sum_dict)
# for _each_ in elements:
#     _ele_str = str(_each_)
#     df_yi_tot['sigma-'+_ele_str] = sigma_iso_ele_sum_eledict[_each_]
#     df_test = pd.DataFrame(sigma_iso_ele_eleisodict[_each_])
#     df_yi_tot = pd.concat([df_yi_tot, df_test], axis=1)
#
# print(df_yi_tot.head())
# # # Export to clipboard
# # df_yi_tot.to_clipboard(excel=True)


### Plot the theoretical neutron resonance
_plot_functions.plot_multi(_energy_x_axis, _trans_y_axis, _plot_mixed, _plot_each_ele_contribution, _plot_each_iso_contribution,
            elements, isotope_dict, x_energy, y_trans_tot, y_ele_dict, y_iso_dicts, _input_formula)

# if _plot_or_not == 'Y':
#     ### Determine x y axis types and captions
#     if _energy_x_axis == 'Y':
#         _x_axis = x_energy
#         _x_words = 'Energy (eV)'
#     else:
#         _x_axis = _functions.ev2lamda(x_energy)
#         _x_words = 'Wavelength (Å)'
#
#     if _trans_y_axis == 'Y':
#         _y_words = 'Neutron transmission'
#     else:
#         _y_words = 'Neutron attenuation'
#
#     ### Determine x y axis values
#     if _plot_mixed == 'Y':
#         if _trans_y_axis == 'Y':
#             _y_axis = y_trans_tot
#         else:
#             _y_axis = 1 - y_trans_tot
#         plt.plot(_x_axis, _y_axis, label=_input_formula)
#
#     if _plot_each_ele_contribution == 'Y':
#         for _ele in elements:
#             _y_each_axis = y_ele_dict[_ele]
#             plt.plot(_x_axis, _y_each_axis, label=_ele)
#
#     if _plot_each_iso_contribution == 'Y':
#         for _ele in elements:
#             for _iso in isotope_dict[_ele]:
#                 _y_each_axis = y_iso_dicts[_ele][_iso]
#                 plt.plot(_x_axis, _y_each_axis, label=_iso)
#
#     plt.ylim(-0.01, 1.01)
#     plt.xlabel(_x_words)
#     plt.ylabel(_y_words)
#     plt.legend(loc='best')
#     plt.show()
