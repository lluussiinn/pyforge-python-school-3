from rdkit import Chem

def substructure_search(mols, mol):
    
    sub_mol = Chem.MolFromSmiles(mol) # defining molecule and substructure to search for
    
    if sub_mol is None: 
        raise ValueError("No substructure")
    
    matching_molecules = []

    for smiles in mols:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is not None and molecule.HasSubstructMatch(sub_mol):
            matching_molecules.append(smiles)


    return matching_molecules


if __name__ == "__main__":
    molecules = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    substructure = "c1ccccc1"
    matching_molecules = substructure_search(molecules, substructure)
    print(matching_molecules)
