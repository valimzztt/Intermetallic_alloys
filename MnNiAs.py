from ase.calculators.espresso import Espresso
from ase.db import connect
from clease.tools import update_db
import os
import re
import shutil

curr_directory = os.getcwd()
pseudo_dir = curr_directory + '/pseudos'
out_dir = curr_directory + '/out'

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






