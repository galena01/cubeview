'''
This example shows how to visualize molecular orbitals using PySCF's cubeview module.
'''


from pyscf import gto, scf, cubeview

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

    mf = scf.UHF(mol)   # Support for UHF and UKS methods
    mf.kernel()

    cv = cubeview.viewer(mf)
    cv.prepare(mo_list=([1,2,5,6], [3,4,7,8]))  # Specify the MOs to visualize for alpha and beta spins, or leave empty to visualize all MOs
    cv.serve()
