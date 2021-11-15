import argparse
import os.path
import re
import Bio.PDB

class Superimposer():
    """ This function aligns two PDB structures using the Bio.PDB module in Biopython.
    Based on Anders Steen Christensen function: https://gist.github.com/andersx/6354971
    TODO check:
    https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_superposition/
    https://github.com/bioexcel/biobb_common/blob/master/biobb_common/tools/test_fixtures.py (!!!)
    TODO hydrogens: test with atom.element param
    """

    def __init__(self, reference, sample, type, output, **kwargs):
        # Input/Output files
        self.reference = reference
        self.sample = sample
        self.output = output
        self.type = type

        self.check_io()

    def is_valid_type(self, t):
        """ Checks if type is correct """
        types = ['structure', 'ligand']
        return t in types

    def is_valid_extension(self, file, label, ext):
        """ Checks if extension is correct """
        filename, file_extension = os.path.splitext(file)
        if file_extension[1:] != ext:
            raise SystemExit('Format %s in %s is not compatible, it must be %s. Exit 1.' % (file_extension[1:], label, ext))
        return True

    def check_io(self):
        """ Checks the reference, sample and output paths and the type """
        # Check reference
        if not os.path.exists(self.reference):
            raise SystemExit('Unexisting reference file. Exit 2.')
        self.is_valid_extension(self.reference, 'reference', 'pdb')
        # Check sample
        if not os.path.exists(self.sample):
            raise SystemExit('Unexisting sample file. Exit 3.')
        self.is_valid_extension(self.sample, 'sample', 'pdb')
        # Check output
        if os.path.dirname(self.output) and not os.path.exists(os.path.dirname(self.output)):
            raise SystemExit('Unexisting output directory. Exit 4.')
        self.is_valid_extension(self.output, 'output', 'pdb')
        # Check type
        if not self.is_valid_type(self.type):
            raise SystemExit('Type %s is not correct. It must be one of: structure, ligand. Exit 5.' % self.type)


    def superimpose_ligands(self):
        """ Aligning two ligands """
        # Start the parser
        pdb_parser = Bio.PDB.PDBParser(QUIET = True)

        # Get the ligands
        ref_structure = pdb_parser.get_structure("reference", self.reference)
        sample_structure = pdb_parser.get_structure("sample", self.sample)

        # Use the first model in the pdb-files for alignment
        ref_model    = ref_structure[0]
        sample_model = sample_structure[0]

        # Get all atoms from both ligands
        all_atoms_ref = ref_model.get_atoms()
        all_atoms_samp = sample_model.get_atoms()

        # Make a list of all the atoms (in the ligands)
        ref_atoms = []
        sample_atoms = []

        # Fill array with atoms of reference ligand
        for atom in all_atoms_ref:
            if not re.search("^H[0-9]{0,}$", atom.get_name()):
                ref_atoms.append(atom)

        # Fill array with atoms of sample ligand
        for atom in all_atoms_samp:
            if not re.search("^H[0-9]{0,}$", atom.get_name()):
                sample_atoms.append(atom)

        # If both ligands don't have the same length, exit
        if not len(ref_atoms) == len(sample_atoms):
            print('Both ligands must have the same number of atoms. Exit 6.')
            raise SystemExit('')

        # Now we initiate the superimposer:
        super_imposer = Bio.PDB.Superimposer()
        super_imposer.set_atoms(ref_atoms, sample_atoms)
        super_imposer.apply(sample_model.get_atoms())

        # Print RMSD:
        print('RMSd: ' + str(super_imposer.rms))

        # Save the aligned version sample
        io = Bio.PDB.PDBIO()
        io.set_structure(sample_structure) 
        io.save(self.output)
        print('Created ' + self.output + ' ligand')

    def superimpose_structures(self):
        """ Aligning two structures """
        # Select what residues numbers you wish to align
        # and put them in a list
        start_id = 1
        end_id   = 70
        atoms_to_be_aligned = range(start_id, end_id + 1)

        # Start the parser
        pdb_parser = Bio.PDB.PDBParser(QUIET = True)

        # Get the structures
        ref_structure = pdb_parser.get_structure("reference", self.reference)
        sample_structure = pdb_parser.get_structure("sample", self.sample)

        # Use the first model in the pdb-files for alignment
        ref_model    = ref_structure[0]
        sample_model = sample_structure[0]

        # Make a list of the atoms (in the structures) you wish to align.
        # In this case we use CA atoms whose index is in the specified range
        ref_atoms = []
        sample_atoms = []

        # Iterate of all chains in the model in order to find all residues
        for ref_chain in ref_model:
          # Iterate of all residues in each model in order to find proper atoms
          for ref_res in ref_chain:
            # Check if residue number ( .get_id() ) is in the list
            if ref_res.get_id()[1] in atoms_to_be_aligned:
              # Append CA atom to list
              ref_atoms.append(ref_res['CA'])

        # Do the same for the sample structure
        for sample_chain in sample_model:
          for sample_res in sample_chain:
            if sample_res.get_id()[1] in atoms_to_be_aligned:
              sample_atoms.append(sample_res['CA'])

        # Now we initiate the superimposer:
        super_imposer = Bio.PDB.Superimposer()
        super_imposer.set_atoms(ref_atoms, sample_atoms)
        super_imposer.apply(sample_model.get_atoms())

        # Print RMSD:
        print('RMSd: ' + str(super_imposer.rms))

        # Save the aligned version sample
        io = Bio.PDB.PDBIO()
        io.set_structure(sample_structure) 
        io.save(self.output)
        print('Created ' + self.output + ' structure')


    def superimpose_structures2(self, rmsd_cutoff=1, remove_hetatm=False, remove_hydrogen=False):
        pdb_a = self.reference
        pdb_b = self.sample
        print("Checking RMSD between:")
        print("     PDB_A: "+pdb_a)
        print("     PDB_B: "+pdb_b)
        pdb_parser = Bio.PDB.PDBParser(QUIET=True)
        st_a = pdb_parser.get_structure("st_a", pdb_a)[0]
        st_b = pdb_parser.get_structure("st_b", pdb_b)[0]

        # Get all atoms from both ligands
        """all_atoms_ref = st_a.get_atoms()
        all_atoms_samp = st_b.get_atoms()"""

        print("     Ignoring HETAMT in RMSD")
        residues_a = [list(res.get_atoms()) for res in st_a.get_residues() if not res.id[0].startswith('H_')]
        residues_b = [list(res.get_atoms()) for res in st_b.get_residues() if not res.id[0].startswith('H_')]
        all_atoms_ref = [atom for residue in residues_a for atom in residue]
        all_atoms_samp = [atom for residue in residues_b for atom in residue]

        # Make a list of all the atoms (in the ligands)
        ref_atoms = []
        sample_atoms = []

        # Fill array with atoms of reference ligand
        for atom in all_atoms_ref:
            if not re.search("^H[0-9]{0,}$", atom.get_name()):
                ref_atoms.append(atom)

        # Fill array with atoms of sample ligand
        for atom in all_atoms_samp:
            if not re.search("^H[0-9]{0,}$", atom.get_name()):
                sample_atoms.append(atom)

        print("     Atoms ALIGNED in PDB_A: "+str(len(ref_atoms)))
        print("     Atoms ALIGNED in PDB_B: "+str(len(sample_atoms)))

        super_imposer = Bio.PDB.Superimposer()
        super_imposer.set_atoms(ref_atoms, sample_atoms)
        super_imposer.apply(sample_atoms)
        print('     RMS: '+str(super_imposer.rms))
        print('     RMS_CUTOFF: '+str(rmsd_cutoff))

        raise SystemExit('Superposition successful. Exit 0.')

        ###########################################

        if remove_hetatm:
            print("     Ignoring HETAMT in RMSD")
            residues_a = [list(res.get_atoms()) for res in st_a.get_residues() if not res.id[0].startswith('H_')]
            residues_b = [list(res.get_atoms()) for res in st_b.get_residues() if not res.id[0].startswith('H_')]
            atoms_a = [atom for residue in residues_a for atom in residue]
            atoms_b = [atom for residue in residues_b for atom in residue]
        else:
            atoms_a = st_a.get_atoms()
            atoms_b = st_b.get_atoms()

        if remove_hydrogen:
            print("     Ignoring Hydrogen atoms in RMSD")
            atoms_a = [atom for atom in atoms_a if not atom.get_name().startswith('H')]
            atoms_b = [atom for atom in atoms_b if not atom.get_name().startswith('H')]

        print("     Atoms ALIGNED in PDB_A: "+str(len(atoms_a)))
        print("     Atoms ALIGNED in PDB_B: "+str(len(atoms_b)))
        super_imposer = Bio.PDB.Superimposer()
        super_imposer.set_atoms(atoms_a, atoms_b)
        super_imposer.apply(atoms_b)
        print('     RMS: '+str(super_imposer.rms))
        print('     RMS_CUTOFF: '+str(rmsd_cutoff))
        return super_imposer.rms < rmsd_cutoff


    def launch(self):
        if(self.type == 'ligand'):
            print('Superimposing ligands...')
            self.superimpose_ligands()
        elif(self.type == 'structure'):
            print('Superimposing structures...')
            self.superimpose_structures()
            #self.superimpose_structures2()

def main():
    parser = argparse.ArgumentParser(description = "Superimposer for two pdb structures.", formatter_class = lambda prog: argparse.RawTextHelpFormatter(prog, width = 99999))
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('-r', '--reference', required = True, help = 'Path for the reference pdb file')
    required_args.add_argument('-s', '--sample', required = True, help = "Path for the sample pdb file")
    required_args.add_argument('-t', '--type', required = True, help = "Type of superimposition. Accepted values: structure, ligand")
    required_args.add_argument('-o', '--output', required = True, help = "Path for the output pdb file")

    args = parser.parse_args()

    Superimposer(reference = args.reference, sample = args.sample, type = args.type, output = args.output).launch()

if __name__ == '__main__':
    main()