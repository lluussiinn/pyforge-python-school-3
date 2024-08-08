import pytest
from fastapi.testclient import TestClient
from main import app, molecules

@pytest.fixture(scope="function", autouse=True)
def clear_molecules():
    """
    Automatically clearing the molecules list before each test.
    """
    molecules.clear()

@pytest.fixture(scope="module")
def test_client():
    return TestClient(app)

# test adding molecules 
@pytest.mark.parametrize("identifier, smiles, expected_status_code, expected_response",[
        ("mol1", "CCO", 201, {"message": "Molecule added successfully", "molecules": [{"identifier": "mol1", "smiles": "CCO"}]}),
        ("mol1", "CCO", 400, {"detail": "Molecule already exists"}),
        ("mol2", "c1ccccc1", 201, {"message": "Molecule added successfully", "molecules": [{"identifier": "mol2", "smiles": "c1ccccc1"}]}),
        ("mol2", "c1ccccc1", 400, {"detail": "Molecule already exists"}),
        ("mol3", "CC(=O)O", 201, {"message": "Molecule added successfully", "molecules": [{"identifier": "mol3", "smiles": "CC(=O)O"}]})])
def test_add_molecule(identifier, smiles, expected_status_code, expected_response, test_client):
    test_client.post("/molecules", json={"identifier": identifier, "smiles": smiles})
    response = test_client.post("/molecules", json={"identifier": identifier, "smiles": smiles})
    assert response.status_code == expected_status_code
    assert response.json() == expected_response

# test retrieving molecules 
@pytest.mark.parametrize("identifier, expected_status_code, expected_response",[
        ("mol1", 200, {"message": "Molecule retrieved successfully", "molecule": {"identifier": "mol1", "smiles": "CCO"}}),
        ("mol2", 200, {"message": "Molecule retrieved successfully", "molecule": {"identifier": "mol2", "smiles": "c1ccccc1"}}),
        ("mol3", 200, {"message": "Molecule retrieved successfully", "molecule": {"identifier": "mol3", "smiles": "CC(=O)O"}}),
        ("nonexistent", 404, {"detail": "Molecule not found"})])

def test_retrieve_molecule(identifier, expected_status_code, expected_response, test_client):
    if identifier in ["mol1","mol2","mol3"]:
        test_client.post("/molecules", json={"identifier": identifier, "smiles": {"mol1": "CCO", "mol2": "c1ccccc1","mol3": "CC(=O)O"}[identifier]})
    response = test_client.get(f"/molecules/{identifier}")
    assert response.status_code == expected_status_code
    assert response.json() == expected_response

# test updating molecules
@pytest.mark.parametrize(
    "identifier, smiles, expected_status_code, expected_response",[
        ("mol1", "CCC", 200, {"message": "Molecule updated successfully", "molecules": [{"identifier": "mol1", "smiles": "CCC"}]}),
        ("mol2", "CCO", 200, {"message": "Molecule updated successfully", "molecules": [{"identifier": "mol2", "smiles": "CCO"}]}),
        ("mol3", "C1=CC=CC=C1", 200, {"message": "Molecule updated successfully", "molecules": [{"identifier": "mol3", "smiles": "C1=CC=CC=C1"}]}),
        ("wrong_input", "CCC", 404, {"detail": "Molecule not found"})])
    
def test_update_molecule(identifier, smiles, expected_status_code, expected_response, test_client):
    if identifier in ["mol1","mol2", "mol3"]:
        test_client.post("/molecules", json={"identifier": identifier, "smiles": {
            "mol1": "CCO","mol2": "c1ccccc1","mol3": "CC(=O)O"}[identifier]})
    response = test_client.put(f"/molecules/{identifier}", json={"identifier": identifier, "smiles": smiles})
    assert response.status_code == expected_status_code
    assert response.json() == expected_response

# test deleting molecules
@pytest.mark.parametrize(
    "identifier, expected_status_code, expected_response",[
        ("mol1", 200, {"message": "Molecule deleted successfully", "molecules": [{"identifier": "mol2", "smiles": "c1ccccc1"}, {"identifier": "mol3", "smiles": "CC(=O)O"}]}),
        ("mol2", 200, {"message": "Molecule deleted successfully", "molecules": [{"identifier": "mol1", "smiles": "CCO"}, {"identifier": "mol3", "smiles": "CC(=O)O"}]}),
        ("mol3", 200, {"message": "Molecule deleted successfully", "molecules": [{"identifier": "mol1", "smiles": "CCO"}, {"identifier": "mol2", "smiles": "c1ccccc1"}]}),
        ("nonexistent", 404, {"detail": "Molecule not found"})])
    
def test_delete_molecule(identifier, expected_status_code, expected_response, test_client):
    if identifier in ["mol1", "mol2", "mol3"]:
        test_client.post("/molecules", json={
            "identifier": identifier, "smiles": {"mol1": "CCO", "mol2": "c1ccccc1", "mol3": "CC(=O)O"}[identifier]})
    response = test_client.delete(f"/molecules/{identifier}")
    assert response.status_code == expected_status_code
    assert response.json() == expected_response

#test getting all molecules
def test_get_all_molecules(test_client):
    for identifier, smiles in {"mol1": "CCO", "mol2": "c1ccccc1", "mol3": "CC(=O)O"}.items():
        test_client.post("/molecules", json={"identifier": identifier, "smiles": smiles})
    response = test_client.get("/molecules")
    assert response.status_code == 200
    assert response.json() == {"message": "Molecules retrieved successfully",
        "molecules": [{"identifier": "mol1", "smiles": "CCO"},
                      {"identifier": "mol2", "smiles": "c1ccccc1"},
                      {"identifier": "mol3", "smiles": "CC(=O)O"}]}

@pytest.mark.parametrize(
    "substructure, expected_status_code, expected_response", [
        ("CC", 200, {"message": "Substructure search completed successfully", "matching_molecules": ["CCO", "CC(=O)O"]}),
        ("c1ccccc1", 200, {"message": "Substructure search completed successfully", "matching_molecules": ["c1ccccc1"]}),
        ("CCCCCC", 200, {"message": "Substructure search completed successfully", "matching_molecules": []}),
        ("invalid_smiles", 400, {"detail": "No substructure"})
    ])
def test_search_substructure(substructure, expected_status_code, expected_response, test_client):
    # Adding test molecules
    for identifier, smiles in {"mol1": "CCO", "mol2": "c1ccccc1", "mol3": "CC(=O)O"}.items():
        test_client.post("/molecules", json={"identifier": identifier, "smiles": smiles})
    
    response = test_client.post("/molecules/search", json={"substructure": substructure})
    assert response.status_code == expected_status_code
    assert response.json() == expected_response