from rdkit import Chem
from fastapi import FastAPI, HTTPException, status, UploadFile
from pydantic import BaseModel
from io import StringIO
import pandas as pd

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


# creating empty list to store dictionaries with identifier and smiles later
molecules = []

class Molecule(BaseModel):
    identifier: str
    smiles: str

app = FastAPI() 

@app.post("/molecules", status_code = status.HTTP_201_CREATED)
def add_molecule(new_molecule: Molecule):
    '''
    Add molecule (smiles) and its identifier
    '''
    for mol in molecules:
        if mol['identifier'] == new_molecule.identifier:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail = "Molecule already exists")
    molecules.append(new_molecule.model_dump()) #dict() didn't work here, so I used model_dump() instead
    return {"message": "Molecule added successfully", "molecules": molecules}


@app.get("/molecules/{identifier}")
def retrieve_molecules(identifier: str):
    '''
    Get molecule by identifier
    '''
    for mol in molecules:
        if mol['identifier'] == identifier:
            return {"message": "Molecule retrieved successfully", "molecule": mol}
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail = "Molecule not found")


@app.put("/molecules/{identifier}")
def update_molecule(identifier: str, updated_mol: Molecule):
    '''
    Updating a molecule by identifier
    '''
    for mol in molecules:
        if mol['identifier'] == identifier:
            mol.update(updated_mol.model_dump())
            return {"message": "Molecule updated successfully", "molecules": molecules}
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail = "Molecule not found")


@app.delete("/molecules/{identifier}")
def delete_molecule(identifier: str):
    '''
    Delete a molecule by identifier
    '''
    global molecules
    for mol in molecules:
        if mol['identifier'] == identifier:
            molecules.remove(mol)
            return {"message": "Molecule deleted successfully", "molecules": molecules}
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Molecule not found")

@app.get("/molecules")
def get_all_molecules():
    '''
    List all molecules
    '''
    return {"message": "Molecules retrieved successfully", "molecules": molecules}

from rdkit import Chem


@app.post("/molecules/search")
def search_substructure(substructure: str):
    '''
    Substructure search for all added molecules
    '''
    smiles_list = [smiles['smiles'] for smiles in molecules]

    try:
        matching_molecules = substructure_search(smiles_list, substructure)
        return {"message": "Substructure search completed successfully", "matching_molecules": matching_molecules}
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@app.post("/molecules/upload")
def upload_file(file: UploadFile):
    ''' 
    Upload csv file with molecules .
    '''
    contents = file.file.read()
    df = pd.read_csv(StringIO(contents.decode('utf-8')))

    if 'identifier' not in df.columns or 'smiles' not in df.columns:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="CSV must contain 'identifier' and 'smiles' columns")
   
    molecules.extend(df.to_dict(orient='records'))
    return {"message": "File uploaded successfully", "filename": file.filename, "molecules": molecules}



