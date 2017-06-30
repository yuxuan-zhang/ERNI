import plot_function

# Parameters
thick_mm = 0.26  # mm
element = 'U'
density_sample = 2  # 0.7875  # pt.elements.isotope(element).density  # g/cm3  https://en.wikipedia.org/wiki/Cadmium
ele_at_ratio = 0.25  # for single element, will be implanted for multiple elements compound
# element2 = '16-O'
mass_abund_other = 0  # pt.elements.isotope(element2).mass * (pt.elements.isotope(element2).abundance/100)
_database = 'ENDF_VIII'
energy_max = 300  # max incident energy in eV
energy_min = 1  # min incident energy in eV
energy_sub = 100
sub_x = energy_sub * (energy_max - energy_min)  # subdivided new x-axis
_type_x_axis = 'energy'  # 1 means plot x-axis as energy in eV
_type_y_axis = 'absorb'  # 1 means plot y-axis as transmission
_plot_each_contribution = 'Y'  # 1 means plot each isotope contribution
_plot_mixed = 'N'  # 1 means plot mixed resonance


plot_function._get_plot(thick_mm, element, density_sample, ele_at_ratio, mass_abund_other, _database,
                        energy_max, energy_min, energy_sub,
                        _type_x_axis, _type_y_axis, _plot_each_contribution, _plot_mixed)
