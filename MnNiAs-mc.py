from ase.calculators.espresso import Espresso
from ase.db import connect
from clease.tools import update_db
import os
import re
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



curr_directory = os.getcwd()
pseudo_dir = curr_directory + '/pseudos'
out_dir = curr_directory + '/out'

# 1) DFT using QUANTUM ESPRESSO

pseudopotentials = {'Mn': 'Mn.pbesol-spn-rrkjus_psl.0.3.1.UPF',
                    'Ni': 'Ni.pbesol-spn-rrkjus_psl.1.0.0.UPF',
                    'As': 'As.pbesol-n-rrkjus_psl.1.0.0.UPF'}

input_data = {
        'control': {
           'calculation': 'relax',
           'restart_mode': 'from_scratch',
           'pseudo_dir': pseudo_dir,
           'prefix': 'tutorial',        
            },
        'system': {
           'ecutwfc': 40.424222147,
           'occupations': 'smearing',
           'degauss': 0.0146997171,
           'lda_plus_u' : '.FALSE.',
          }
        }

calc = Espresso(pseudopotentials = pseudopotentials,
                input_data = input_data,
                kspacing=0.4)
cwd = os.getcwd()

#Connecting to the relevant database
db_name = "clease_MnNiAs.db"
db = connect(db_name)

"""Function that extracts the energy from the espresso.pwo file
 in case the get_potential_energy() method from ASE fails"""
def get_energies(filename):
    energies = []
    with open(filename,"r") as file:
        for line in file:
            pattern="!"
            if re.search(pattern, line):
                energies.append(re.findall(r'-?\d+', line))
        # Get the last calculated energy value before the calculation converges
        a = energies[-1][0]
        b = energies[-1][1]
        d = float(f'{a}.{b}') 
        #Important: we need to convert back from Rydberg to eV
        d = 13.6057039763*d
        print("The final energies is ", d)
        return d
        
"""Updates the energy column for the given Atoms object given the id of the Atoms row (first column)"""
def update_energy(energy, id):
    with db.managed_connection() as con:
            cur = con.cursor()
            cur.execute(
                'UPDATE systems SET energy=? WHERE id=?',
                (energy, id))
            
"""Deletes from the database the configuration for which DFT calculations did not converge in 100 iterations"""
def handle_non_convergence(id):
    with db.managed_connection() as con:
            cur = con.cursor()
            cur.execute(
                'DELETE FROM systems WHERE id=?',
                (id,))

"""Updates the energy column for the given Atoms object given initial and final struct id once DFT is done"""
def update_converged_configs(initial_id, final_id):
    # get the id of the final structure whose calculation converged
    with db.managed_connection() as con:
            cur = con.cursor()
            cur.execute("SELECT energy FROM systems WHERE id = ? AND energy IS NOT NULL", (initial_id,))
            # Fetch the result
            result = cur.fetchone()
            if result != None: 
                energy = result[0]
                cur.execute('UPDATE systems SET energy = (SELECT energy FROM systems WHERE id =? AND energy IS NOT NULL) WHERE id = ?', (initial_id, final_id))
                cur.execute("UPDATE systems SET energy = NULL WHERE id = ?", (initial_id,))


for row in db.select(converged=False):
    db.update(row.id, queued=True)
    print("We are analzying the following configuration")
    print(row.id, row.name)
    atoms = row.toatoms()
    atoms.calc = calc
    directory = 'config_{}'.format(row.id)
    parent_dir = cwd
    path = os.path.join(parent_dir, directory)
    #if the first test fails, this try/except ensures that the calculations do not interfere with each other
    try:
        if not os.path.exists(path):
            os.mkdir(path)   
        # this was here just for testing purposes         
        #shutil.copy("espresso.pwo",path)
    except OSError as error:
        print(error)
        continue
    os.chdir(path)
    try: 
         # For some reason, this method from ASE
        # will sometimes fail so we will get the energies computed by QE manually by reading the output file
        atoms.get_potential_energy()
    except:
        # For some reason, it will fail so we will just get the energy value manually from the created espresso,pwo file
        try:
            energy = get_energies("espresso.pwo")
            update_energy(energy, row.id)
        # we need to handle the case where the DFT calculations do not converge: discard the configuration and delete
        # it from the database
        except: 
            handle_non_convergence(row.id)
            continue

    os.chdir(cwd)
    update_db(uid_initial=row.id, final_struct=atoms, db_name=db_name)
    print("We have updated the database for the following row")
    print(row.id)

"""Fine tunement to database once DFT has run successfully"""
for row in db.select(converged=True):
      dictionary = row.key_value_pairs
      final_struct_id = row.key_value_pairs.get("final_struct_id")
      initial_struct_id = row.id
      update_converged_configs(initial_struct_id, final_struct_id)  


from mdutils.mdutils import MdUtils
from mdutils import Html
mdAcme = MdUtils(file_name='info')
mdAcme.create_md_file()
# 2) CLUSTER EXPANSION: 
from clease import Evaluate
eva = Evaluate(settings=settings, scoring_scheme='loocv')

# scan different values of alpha and return the value of alpha that yields
# the lowest CV score
eva.set_fitting_scheme(fitting_scheme='l2')
#eva.set_fitting_scheme(fitting_scheme='l2', alpha=0.0001)
#alpha, cv = eva.alpha_CV(alpha_min=1E-6, alpha_max=10.0, num_alpha=50)
alpha = eva.plot_CV(alpha_min=1E-6, alpha_max=100.0, num_alpha=50)
eva.set_fitting_scheme(fitting_scheme='l2', alpha=alpha)
print("The chosen value of alpha is")
print(alpha)
eva.fit()

# plot ECI values
#eva.plot_ECI()
import clease.plot_post_process as pp
import json 
import matplotlib.pyplot as plt
import json
fig = pp.plot_fit(eva)
fig.savefig('alpha-fit.png')
mdAcme.new_paragraph('CE fitting results:', bold_italics_code='bi', align='center')
mdAcme.new_header(level=1, title='Overview')
mdAcme.write("The chosen value of alpha is: ")
mdAcme.write(str(alpha))
mdAcme.new_paragraph(Html.image(path='alpha-fit.png', size='250', align='center'))


# plot ECI values
fig = pp.plot_eci(eva)
fig.savefig('ECI-values.png')
mdAcme.new_paragraph(Html.image(path='ECI-values.png', size='250', align='center'))
mdAcme.create_md_file()

# save a dictionary containing cluster names and their ECIs
cluster_expansion = curr_directory + '/MnNiAs'
eva.save_eci(fname=cluster_expansion)
eci_file=open('MnNiAs.json', 'r')
eci = json.load(eci_file)
print(eci)
from clease.calculator import attach_calculator
atoms = settings.atoms.copy()*(7, 7, 7)
atoms = attach_calculator(settings, atoms=atoms, eci=eci)

for i in range(0,len(atoms),1):
    if atoms[i].symbol=='Mn' and i%2 == 0:
        atoms[i].symbol='Ni'


from clease.montecarlo import Montecarlo
from clease.montecarlo.observers import EnergyEvolution
# Monte Carlo results will be createed
directory = 'MC_5000K_results'
parent_dir = cwd
path = os.path.join(parent_dir, directory)
os.mkdir(path)  
os.chdir(path)

T = 5000
nsteps=1372000
mc = Montecarlo(atoms, T)

from clease.montecarlo.constraints import FixedElement
cnst = FixedElement('O')
mc.generator.add_constraint(cnst)

from clease.montecarlo.observers import EnergyEvolution, Snapshot, CorrelationFunctionObserver, AcceptanceRate
obs = EnergyEvolution(mc)
mc.attach(obs, interval=1000)
corr = CorrelationFunctionObserver(atoms.calc)
mc.attach(corr)
accrate = AcceptanceRate()
mc.attach(accrate)
snap = Snapshot(atoms, fname='snapshot')
mc.attach(snap, interval=137200)
mc.run(steps=nsteps)
energies = obs.energies
thermo = mc.get_thermodynamic_quantities()

f = open('MC_5000','w')
print(thermo, file=f)
print(energies, sep='\n', file=f)
f.close()



for i in range(5000, 0, -50):
    T = i
    #atoms = vasp.read_vasp(POSCARstring)
    #atoms = attach_calculator(settings, atoms=atoms, eci=eci)
    mc = Montecarlo(atoms, T)
    mc.generator.add_constraint(cnst)
    obs = EnergyEvolution(mc)
    mc.attach(obs, interval=1000)
    snap = Snapshot(atoms, fname='snapshot')
    mc.attach(snap, interval=137200)
    corr = CorrelationFunctionObserver(atoms.calc)
    mc.attach(corr)
    accrate = AcceptanceRate()
    mc.attach(accrate)
    
    mc.run(steps=137200)
    mc.run(steps=nsteps)

    energies = obs.energies
    rate=accrate.get_averages()
    thermo = mc.get_thermodynamic_quantities()
    corr_dict=corr.get_averages()
    
    from ase.build import sort
    from ase.io import vasp, db
    import pickle

    sorted_atoms=sort(atoms)
    vasp.write_vasp('POSCAR{}.vasp'.format(i), sorted_atoms, direct=False, wrap=False)
    POSCARstring='POSCAR{}.vasp'.format(i)
    
    f = open('MC_{}'.format(i),'w')
    print(thermo, file=f)
    print(energies, sep='\n', file=f)
    f.close()


os.chdir(cwd)



