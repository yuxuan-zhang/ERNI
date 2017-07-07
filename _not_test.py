import _functions

_input = 'UO3'
thick_mm = .3  # mm
_input_density = 0.7875  # g/cm3  not needed if _input is single element
_database = 'ENDF_VIII'
energy_max = 300  # max incident energy in eV
energy_min = 0  # min incident energy in eV
energy_sub = 100
sub_x = energy_sub * (energy_max - energy_min)  # subdivided new x-axis
# _multi_element = 'N'
_energy_x_axis = 'Y'  # 1 means plot x-axis as energy in eV
_trans_y_axis = 'N'  # 1 means plot y-axis as transmission
_plot_each_ele_contribution = 'Y'  # 1 means plot each element's contribution
_plot_each_iso_contribution = 'N'  # 1 means plot each isotope's contribution
_plot_mixed = 'N'  # 1 means plot mixed resonance

formula_dict = _functions.input2formula(_input)  # Function called to parse input formula and return elements and ratios
elements = _functions.dict_key_list(formula_dict)
ratios = _functions.dict_value_list(formula_dict)
natural_ele_boo_dict = _functions.boo_dict(elements)  # Dict for natural mixture
thick_mm_dict = _functions.thick_dict(elements, thick_mm)

# For unnatural mixture elements:
_natural_ele_input = 'Y'  # input('Is there any unnatural mixture? ')
if _natural_ele_input == 'Y':
    unnatural_ratio_dicts = {}
    unnatural_element_str = 'U,O'  # input('Please list all separated by only ",": ')
    unnatural_element = unnatural_element_str.split(',')
    natural_ele_boo_dict = _functions.boo_dict_invert_by_key(unnatural_element, natural_ele_boo_dict)
    isotope_dict = _functions.get_isotope_dict(_database, unnatural_element)
    print(isotope_dict)
    for ele in unnatural_element:
        isotopes = isotope_dict[ele]
        unnatural_ratio_dict = {}
        for iso in isotopes:
            unnatural_ratio_dict[iso] = float(input('Atomic ratio of {} in mixture: '.format(iso)))
        unnatural_ratio_dicts[ele] = unnatural_ratio_dict
    print(unnatural_ratio_dicts)
