from pymol import cmd


# Function to do the point mutations in a pdb file using PyMol
def pdb_creator(file_name, chain_name, atom_number, mutation):
    # pymol.finish_launching()
    atom_selected = chain_name + '/' + atom_number + '/'

    # Load molecule
    cmd.fetch(file_name)
    cmd.select(file_name)
    cmd.show("licorice")
    cmd.hide("cartoon")
    cmd.deselect()
    cmd.save("imported_pdb.png")

    # Image of un-mutated residue
    cmd.select(file_name)
    cmd.show("licorice")
    cmd.hide("cartoon")
    cmd.select("mutant", atom_selected)
    cmd.zoom("mutant", "4")
    cmd.show("licorice")
    cmd.deselect()
    cmd.save("not_mutated.png")

    # Saving pdb file
    cmd.select("hetera", "hetatm")
    cmd.alter("hetera", "type='ATOM'")
    cmd.alter("mutant", "type='HETATM'")
    cmd.save("altered.pdb")

    # Clearing Screen
    cmd.reinitialize()
    cmd.fetch(file_name)
    cmd.select("mutant", atom_selected)

    # Turn on the wizard
    cmd.wizard("mutagenesis")
    cmd.do("refresh_wizard")

    # Apply mutation
    cmd.get_wizard().set_mode(mutation)
    cmd.get_wizard().do_select(atom_selected)
    cmd.get_wizard().apply()

    # Saving Mutated pdb
    cmd.select("hatm", "hetatm")
    cmd.alter("hatm", "type='ATOM'")
    cmd.select("new", atom_selected)
    cmd.alter("new", "type='HETATM'")
    cmd.save("mutated.pdb")

    # Image of Mutated Residue
    cmd.select("mutant", atom_selected)
    cmd.zoom("mutant", "4")
    cmd.show("licorice")
    cmd.hide("cartoon")
    cmd.deselect()
    cmd.save("mutated.png")

    # Quitting PyMol
    cmd.quit()


FILE_NAME = input("Enter Existing PDB : ")
CHAIN_NAME = input("Enter Name of the chain to be modified : ")
ATOM_NUMBER = input("Enter the atom number : ")
MUTATION = input("Enter the required Mutation (AA Name): ")

pdb_creator(FILE_NAME, CHAIN_NAME, ATOM_NUMBER, MUTATION)
