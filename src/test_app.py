import pytest
from fastapi.testclient import TestClient
from main import app  
import pandas as pd
from io import StringIO

client = TestClient(app)

# Test for adding a molecule

def test_add_molecule():
    response = client.post("/molecules", json={"identifier": "test1", "smiles": "CCO"})
    assert response.status_code == 201
    assert response.json() == {"message": "Molecule added successfully", "molecules": [{"identifier": "test1", "smiles": "CCO"}]}

def test_add_molecule_duplicate():
    client.post("/molecules", json={"identifier": "test2", "smiles": "CCO"})
    response = client.post("/molecules", json={"identifier": "test2", "smiles": "CCO"})
    assert response.status_code == 400
    assert response.json() == {"detail": "Molecule already exists"}

# Test for retrieving a molecule
def test_retrieve_molecule():
    client.post("/molecules", json={"identifier": "test3", "smiles": "CCC"})
    response = client.get("/molecules/test3")
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule retrieved successfully", "molecule": {"identifier": "test3", "smiles": "CCC"}}

def test_retrieve_molecule_not_found():
    response = client.get("/molecules/nonexistent")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}

# Test for updating a molecule
def test_update_molecule():
    client.post("/molecules", json={"identifier": "test4", "smiles": "CCCC"})
    response = client.put("/molecules/test4", json={"identifier": "test4", "smiles": "CCO"})
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule updated successfully", "molecules": [{"identifier": "test4", "smiles": "CCO"}]}

def test_update_molecule_not_found():
    response = client.put("/molecules/nonexistent", json={"identifier": "nonexistent", "smiles": "CCO"})
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}

# Test for deleting a molecule
def test_delete_molecule():
    client.post("/molecules", json={"identifier": "test5", "smiles": "CCO"})
    response = client.delete("/molecules/test5")
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule deleted successfully", "molecules": []}

def test_delete_molecule_not_found():
    response = client.delete("/molecules/nonexistent")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}

# Test for getting all molecules
def test_get_all_molecules():
    client.post("/molecules", json={"identifier": "test6", "smiles": "C1CCCCC1"})
    response = client.get("/molecules")
    assert response.status_code == 200
    assert response.json() == {"message": "Molecules retrieved successfully", "molecules": [{"identifier": "test6", "smiles": "C1CCCCC1"}]}

# Test for substructure search
def test_search_substructure():
    client.post("/molecules", json={"identifier": "test7", "smiles": "C1CCOCC1"})
    response = client.post("/molecules/search", json={"substructure": "O"})
    assert response.status_code == 200
    assert response.json() == {"message": "Substructure search completed successfully", "matching_molecules": ["C1CCOCC1"]}

def test_search_substructure_no_match():
    response = client.post("/molecules/search", json={"substructure": "N"})
    assert response.status_code == 200
    assert response.json() == {"message": "Substructure search completed successfully", "matching_molecules": []}

# Test for uploading a file
def test_upload_file():
    csv_content = "identifier,smiles\nfile1,C1CCCCC1\n"
    file = ("molecules.csv", csv_content, "text/csv")
    response = client.post("/molecules/upload", files={"file": file})
    assert response.status_code == 200
    assert response.json() == {"message": "File uploaded successfully", "filename": "molecules.csv", "molecules": [{"identifier": "file1", "smiles": "C1CCCCC1"}]}

def test_upload_file_invalid():
    invalid_csv_content = "invalid_column,smiles\nfile2,C1CCCCC1\n"
    file = ("invalid.csv", invalid_csv_content, "text/csv")
    response = client.post("/molecules/upload", files={"file": file})
    assert response.status_code == 400
    assert response.json() == {"detail": "CSV must contain 'identifier' and 'smiles' columns"}