import _one

folis_thickmm_dict = {'Gd': 1,
                      'Cd': 1,
                      'U': 1,
                      'Ag': 1,
                      'Au': 1,
                      'Ta': 1,
                      'In': 1,
                      'W': 1,
                      'Co': 1,
                      'Hf': 1,
                      'Pb': 1,
                      'Eu': 1,
                      'B': 1
                      }
# Parameters
_input = 'U'
thick_mm = .26  # mm
_natural_ele = ['Y']
_input_density = 0.7875  # not needed if _input is single element (Unit: g/cm3)
_input_ratios_dict = {}
                      #'U': [0, 0, 0.15, 0.85]
                      # 'O': [1, 0, 0]
                      # at.%  not needed if _natural_ele is 'Y'
_database = 'ENDF_VIII'
energy_max = 300  # max incident energy in eV
energy_min = 0  # min incident energy in eV
energy_sub = 100
_energy_x_axis = 'Y'  # 1 means plot x-axis as energy in eV
_trans_y_axis = 'N'  # 1 means plot y-axis as transmission
_plot_each_ele_contribution = 'Y'  # Y means plot each element's contribution
_plot_each_iso_contribution = 'N'  # Y means plot each isotope's contribution
_plot_mixed = 'N'  # Y means plot mixed resonance

_one.plot_input(_input, _natural_ele, thick_mm, _input_density, _input_ratios_dict, _database,
                energy_max, energy_min, energy_sub, _energy_x_axis, _trans_y_axis, _plot_mixed,
                _plot_each_ele_contribution, _plot_each_iso_contribution)