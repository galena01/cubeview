# Visualize MOs from a Molcas calculation.
# Use `export MOLCAS_MOLDEN=ON` before running Molcas to generate a Molden file.

# Here's an example input for Molcas:
''' c2h4_guessorb.inp
&GATEWAY
Title = C2H4 molecule as an example for generating Molden file
Coord
6

C     0.      0.      0.658 
C     0.      0.     -0.658 
H     0.      0.924   1.252 
H     0.     -0.924   1.252 
H     0.     -0.924  -1.252 
H     0.      0.924  -1.252 
Basis = 6-31G*
Group = X Y Z

&SEWARD
'''
# Here we only do seward, which computes the integrals and generates a initial guess without doing SCF calculation, we can visualize the MOs from this step.



from pyscf import gto, scf, cubeview
import re
import numpy as np


def parse_molden(filename):
    atoms = []
    mos = []
    with open(filename, "r") as f:
        lines = f.readlines()

    # Parse atoms
    in_atoms = False
    for line in lines:
        if line.strip().startswith("[Atoms]"):
            in_atoms = True
            continue
        if in_atoms:
            if line.strip().startswith("["):  # End of atoms section
                in_atoms = False
                continue
            parts = line.split()
            if len(parts) >= 6:
                label = parts[0]
                idx = int(parts[1])
                Z = int(parts[2])
                x, y, z = map(float, parts[3:6])
                atoms.append((label, idx, Z, x, y, z))

    # Parse MOs
    in_mo = False
    mo_data = {}
    coeffs = []
    for line in lines:
        if line.strip().startswith("[MO]"):
            in_mo = True
            continue
        if in_mo:
            if line.strip().startswith("Sym="):
                if mo_data:  # Start of a new MO, save the previous one
                    mo_data["coeffs"] = coeffs
                    mos.append(mo_data)
                    mo_data = {}
                    coeffs = []
                mo_data["sym"] = line.split("=")[1].strip()
            elif line.strip().startswith("Ene="):
                mo_data["ene"] = float(line.split("=")[1].strip())
            elif line.strip().startswith("Spin="):
                mo_data["spin"] = line.split("=")[1].strip()
            elif line.strip().startswith("Occup="):
                mo_data["occup"] = float(line.split("=")[1].strip())
            elif re.match(r"^\s*\d+", line):  # Coefficient lines
                parts = line.split()
                if len(parts) == 2:
                    idx = int(parts[0])
                    val = float(parts[1])
                    coeffs.append(val)
    # Add the last MO if exists
    if mo_data:
        mo_data["coeffs"] = coeffs
        mos.append(mo_data)

    # Sort by mo energy
    mos.sort(key=lambda x: x["ene"])

    # Convert to numpy arrays
    n_basis = len(mos[0]["coeffs"])
    n_mo = len(mos)

    mo_coeff = np.zeros((n_basis, n_mo))
    mo_energy = np.zeros(n_mo)
    mo_occ = np.zeros(n_mo)
    mo_spin = []
    mo_sym = []

    for i, mo in enumerate(mos):
        mo_coeff[:, i] = mo["coeffs"]
        mo_energy[i] = mo["ene"]
        mo_occ[i] = mo["occup"]
        mo_spin.append(mo["spin"])
        mo_sym.append(mo["sym"])

    return atoms, mo_coeff, mo_energy, mo_occ, mo_spin, mo_sym

if __name__=="__main__":

    # Rebuild the molecule in PySCF to generate the cube files.
    # MAKE SURE THE COORDINATES AND BASIS SET ARE THE SAME AS IN MOLCAS
    # Here we set cart=True because Molcas uses Cartesian basis functions.
    mol = gto.M(
    atom = """
 C     0.      0.      0.658 
 C     0.      0.     -0.658 
 H     0.      0.924   1.252 
 H     0.     -0.924   1.252 
 H     0.     -0.924  -1.252 
 H     0.      0.924  -1.252 
""",
    basis = "6-31g*", cart = True,
    spin = 0, charge = 0, symmetry = "D2h",
    verbose = 3)

    atoms, mo_coeff, mo_energy, mo_occ, mo_spin, mo_sym = parse_molden('04-molcas_files/c2h4_guessorb.guessorb.molden')

    # This example assumes a restricted system, for unrestricted cases, you need to handle alpha and beta orbitals separately.

    cv = cubeview.viewer(
        mf_or_mol = mol,
        mo_coeff = mo_coeff,
        mo_energy = mo_energy,
        mo_occ = mo_occ,
        irreps = mo_sym,
    )
    cv.prepare(ao_component=2)  

    cv.serve()