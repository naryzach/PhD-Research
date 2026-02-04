import sys
from Bio.PDB import Structure, Model, Chain, Residue, Atom, PDBIO

def create_mock_structure(chain_id):
    s = Structure.Structure("test")
    m = Model.Model(0)
    c = Chain.Chain(chain_id)
    r = Residue.Residue((" ", 1, " "), "ALA", " ")
    a = Atom.Atom("CA", [0,0,0], 20.0, 1.0, " ", "CA", 0, "C")
    r.add(a)
    c.add(r)
    m.add(c)
    s.add(m)
    return s

def try_save(chain_id):
    print(f"Testing save with chain ID: '{chain_id}'")
    s = create_mock_structure(chain_id)
    io = PDBIO()
    io.set_structure(s)
    try:
        io.save("test_out.pdb")
        print("Save successful")
        return True
    except Exception as e:
        print(f"Save failed as expected: {e}")
        return False

if __name__ == "__main__":
    # Test valid
    try_save("A")
    # Test invalid
    if not try_save("AA"):
        sys.exit(0) # Success if it fails
    else:
        print("Error: 'AA' chain should have failed but didn't.")
        sys.exit(1)
