import markdown
from mdutils.mdutils import MdUtils
from mdutils import Html
import matplotlib.pyplot as plt
import json
import clease.plot_post_process as pp


from clease.settings import CECrystal
from clease.settings import Concentration

# 0) setup crystal structure
#define concentration range of all elements
conc = Concentration(basis_elements=[['Mn', 'Ni'], ['As']])
conc.set_conc_ranges(ranges=[[(0,1),(0,1)], [(1,1)]])

#define crystal structure
settings = CECrystal(concentration=conc, 
    spacegroup=194, 
    basis=[(0.00000, 0.00000, 0.00000), (0.33333333, 0.66666667, 0.25)], 
    cell=[3.64580405, 3.64580405,   5.04506600, 90, 90, 120],  
    supercell_factor=8, 
    db_name="clease_MnNiAs.db", 
    basis_func_type='binary_linear', 
    max_cluster_dia=(7,7,7))


mdAcme = MdUtils(file_name='info')
mdAcme.create_md_file()

from clease import Evaluate
eva = Evaluate(settings=settings, scoring_scheme='loocv')
# scan different values of alpha and return the value of alpha that yields
# the lowest CV score
eva.set_fitting_scheme(fitting_scheme='l2')
alpha = eva.plot_CV(alpha_min=1E-6, alpha_max=10.0, num_alpha=50)
eva.set_fitting_scheme(fitting_scheme='l2', alpha=alpha)
print("The chosen value of alpha is")
print(alpha)
eva.fit()

fig = pp.plot_fit(eva)
plt.show()

fig.savefig('alpha-fit.png')
mdAcme.new_paragraph('License and Installation Instructions', bold_italics_code='bi', align='center')
mdAcme.new_header(level=1, title='Overview')
mdAcme.write("Welcome to <font color='red'>Acme Spinners</font>!\n\n")
mdAcme.new_paragraph(Html.image(path='alpha-fit.png', size='250', align='center'))
# plot ECI values
fig = pp.plot_eci(eva)
plt.show()
fig.savefig('ECI-values.png')
mdAcme.new_paragraph(Html.image(path='ECI-values.png', size='250', align='center'))
mdAcme.create_md_file()