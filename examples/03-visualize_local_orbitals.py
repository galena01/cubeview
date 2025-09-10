'''
This example shows how to visualize local orbitals using PySCF's cubeview module.
'''

from pyscf import gto, scf, lo, cubeview

if __name__=="__main__":
    mol = gto.M(
    atom = """
 C     0.      0.      0.658 
 C     0.      0.     -0.658 
 H     0.      0.924   1.252 
 H     0.     -0.924   1.252 
 H     0.     -0.924  -1.252 
 H     0.      0.924  -1.252 
""",
    basis = "6-31g*",
    spin = 0, charge = 0, symmetry = "D2h",
    verbose = 3)

    mf = scf.RHF(mol)
    mf.kernel()

    
    c_lo = lo.orth_ao(mf, 'nao')
    # print(c_lo.shape)

    cv = cubeview.viewer(mol, mo_coeff=c_lo)
    cv.prepare(ao_component=2)  # Render all orbitals, and show the two largest components of the local orbitals
    
    # Note: ao_component is calculated directly from the square of MO coefficients, which is not a awfully reliable method for analyzing molecular orbitals. PR is welcomed if you want to implement methods like Hirshfeld, AIM, etc.

    # Occupation numbers, Irreps, and MO energies are NOT available for local orbitals.

    cv.serve()
